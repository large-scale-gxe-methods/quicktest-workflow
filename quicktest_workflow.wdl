task process_phenos {
	
	File phenofile
	File samplefile
	String sample_id_header
	String outcome
	String exposure
	String covar_names
	String delimiter
	String missing
	Int ppmem

	command {
		python3 /format_quicktest_phenos.py ${phenofile} ${sample_id_header} ${outcome} "${exposure}" "${covar_names}" "${delimiter}" ${missing} "${samplefile}"
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/quicktest-workflow"
		memory: ppmem + "GB"
	}

        output {
                File pheno_fmt = "quicktest_phenotypes.csv"
		File covar_file = "covar_string.txt"
	}
}

task run_interaction {

	File genofile
	Boolean is_bgen
	File phenofile
	String outcome
	Boolean binary_outcome
	File covar_file
	String? missing = "NA"
	Boolean robust
	Int? memory = 10
	Int? disk = 20
	Int monitoring_freq

	String covar_string = read_string(covar_file)

	command {
		dstat -c -d -m --nocolor 1 > system_resource_usage.log &
		atop -x -P PRM 1 | grep '(quicktest)' > process_resource_usage.log &

		/quicktest-1.1_bgen_v1.2/quicktest \
			--geno ${genofile} \
			${true="--bgen" false="" is_bgen} \
			--pheno ${phenofile} \
			--npheno ${outcome} \
			${covar_string} \
			--method-mean \
			${true="--method-binary 0.5" false="" binary_outcome} \
			--method-interaction \
			--missing-code ${missing} \
			${true="--method-robust" false="" robust} \
			--out quicktest_res \
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/quicktest-workflow"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
		dx_timeout: "7D0H00M"
	}

	output {
        File res = "quicktest_res"
	File system_resource_usage = "system_resource_usage.log"
	File process_resource_usage = "process_resource_usage.log"
    }
}

task standardize_output {

	File resfile
    Boolean robust

	String outfile_base = basename(resfile)
	String outfile = "${outfile_base}.fmt"

	command {
		python3 /format_quicktest_output.py ${resfile} ${outfile} ${robust}
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/quicktest-workflow"
		memory: "2 GB"
	}

	output {
		File res_fmt = "${outfile}"
	}
}

task cat_results {

	Array[File] results_array

	command {
		head -1 ${results_array[0]} > all_results.txt && \
			for res in ${sep=" " results_array}; do tail -n +2 $res >> all_results.txt; done
	}
	
	runtime {
		docker: "quay.io/large-scale-gxe-methods/probabel-workflow"
		disks: "local-disk 5 HDD"
	}
	output {
		File all_results = "all_results.txt"
	}
}


workflow run_quicktest {

	Array[File] genofiles
	Boolean is_bgen = false
	Float? maf
	Array[File] samplefiles
	File phenofile
	String? sample_id_header = "sampleid"
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names = ""
	String? delimiter = ","
	String? missing = "NA"
	Boolean robust
	Int? memory
	Int? disk
	Int? monitoring_freq = 1

	Int ppmem = 2 * ceil(size(phenofile, "GB")) + 1

	scatter (samplefile in samplefiles) {
		call process_phenos {
			input:
				phenofile = phenofile,
				samplefile = samplefile,
				sample_id_header = sample_id_header,
				outcome = outcome,
				exposure = exposure_names,
				covar_names = covar_names,
				delimiter = delimiter,
				missing = missing,
				ppmem = ppmem
		}	
	}

	scatter (i in range(length(genofiles))) {
                call run_interaction {
                        input:
                                genofile = genofiles[i],
				is_bgen = is_bgen,
                                phenofile = process_phenos.pheno_fmt[i],
                                covar_file = process_phenos.covar_file[i],
				outcome = outcome,
                                binary_outcome = binary_outcome,
				missing = missing,
                                robust = robust,
                                memory = memory,        
                                disk = disk,
				monitoring_freq = monitoring_freq
                }
        }

	scatter (resfile in run_interaction.res) {
		call standardize_output {
			input:
				resfile = resfile,
                robust = robust
		}
	}	

	call cat_results {
		input:
			results_array = standardize_output.res_fmt
	}

        output {
                File results = cat_results.all_results
		Array[File] system_resource_usage = run_interaction.system_resource_usage
		Array[File] process_resource_usage = run_interaction.process_resource_usage
	}

	parameter_meta {
		genofiles: "Array of genotype filepaths in Oxford (.gen.gz) format (may be gzipped)."
		is_bgen: "Optional boolean to specify whether the input genotype file is in .bgen format (default is False, i.e. .gen format)."
		maf: "Minor allele frequency threshold for pre-filtering variants as a fraction (default is 0.001)."
		samplefile: "Optional .sample file accompanying the .gen/.bgen file. Required for proper function if .gen/.bgen does not store sample identifiers."
		phenofile: "Phenotype filepath."	
		sample_id_header: "Optional column header name of sample ID in phenotype file."
		outcome: "Column header name of phenotype data in phenotype file."
                binary_outcome: "Boolean: is the outcome binary? Otherwise, quantitative is assumed."
		exposure_names: "Column header name(s) of the exposures for genotype interaction testing (space-delimited). Only one exposure is currently allowed."
		covar_names: "Column header name(s) of any covariates for which only main effects should be included selected covariates in the pheno data file (space-delimited). This set should not overlap with exposure_names."
		delimiter: "Delimiter used in the phenotype file."
		missing: "Missing value key of phenotype file."
                robust: "Boolean: should robust (a.k.a. sandwich/Huber-White) standard errors be used?"
		memory: "Requested memory (in GB)."
		disk: "Requested disk space (in GB)."
		monitoring_freq: "Delay between each output for process monitoring (in seconds). Default is 1 second."
	}

        meta {
                author: "Ye Chen"
                email: "ychen115@mgh.harvard.edu"
                description: "Run gene-environment interaction tests using the QuickTest program."
        }
}
