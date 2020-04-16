task process_phenos {
	
	File phenofile
	File samplefile
	String sample_id_header
	String outcome
	String exposure
	String? covar_names = ""
	String? delimiter = ","
	String? missing = "NA"

	command {
		python3 /format_quicktest_phenos.py ${phenofile} ${sample_id_header} ${outcome} "${exposure}" "${covar_names}" "${delimiter}" ${missing} "${samplefile}"
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/quicktest-workflow"
		memory: "2 GB"
	}

        output {
                File pheno_fmt = "quicktest_phenotypes.csv"
		File covar_file = "covar_string.txt"
	}
}

task run_interaction {

	File genofile
	Boolean? is_bgen = false
	File phenofile
	String outcome
	Boolean binary_outcome
	File covar_file
	String? missing = "NaN"
	Boolean robust
	Int? memory = 10
	Int? disk = 20

	String covar_string = read_string(covar_file)

	command {
		echo "" > resource_usage.log
		dstat -c -d -m --nocolor 10 1>>resource_usage.log &
		/quicktest-1.1_bgen_v1.2/quicktest \
			--geno ${genofile} \
			${default="" true="--bgen" false="" is_bgen} \
			--pheno ${phenofile} \
			--npheno ${outcome} \
			${covar_string} \
			--method-mean \
			${default="" true="--method-binary 0.5" false="" binary_outcome} \
			--method-interaction \
			--missing-code ${missing} \
			${default="" true="--method-robust" false="" robust} \
			--out quicktest_res \
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/quicktest-workflow"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
	}

	output {
        File res = "quicktest_res"
	File resource_usage = "resource_usage.log"
    }
}

task standardize_output {

	File resfile
	String outfile_base = basename(resfile)
	String outfile = "${outfile_base}.fmt"

	command {
		python3 /format_quicktest_output.py ${resfile} ${outfile}
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
	Boolean? is_bgen
	Float? maf
	Array[File] samplefiles
	File phenofile
	String? sample_id_header = "sampleid"
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names
	String? delimiter
	String? missing
	Boolean robust
	Int? memory
	Int? disk

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
				missing = missing
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
                                disk = disk
                }
        }

	scatter (resfile in run_interaction.res) {
		call standardize_output {
			input:
				resfile = resfile
		}
	}	

	call cat_results {
		input:
			results_array = standardize_output.res_fmt
	}

        output {
                File results = cat_results.all_results
		Array[File] resource_usage = run_interaction.resource_usage
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
	}

        meta {
                author: "Ye Chen"
                email: "ychen115@mgh.harvard.edu"
                description: "Run gene-environment interaction tests using the QuickTest program."
        }
}
