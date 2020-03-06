task process_phenos {
	
	File phenofile
	String sample_id_header
	String outcome
	String covar_headers
	String exposure
	String? delimiter = ","
	String? missing = "NA"

	command {
		python3 /format_quicktest_phenos.py ${phenofile} ${sample_id_header} ${outcome} "${covar_headers}" ${exposure} "${delimiter}" ${missing}
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
	File? samplefile
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
			--bgen \
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
	String exposure
	String outfile_base = basename(resfile)
	String outfile = "${outfile_base}.fmt"

	command {
		python3 /format_quicktest_output.py ${resfile} ${exposure} ${outfile}
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
	Float? maf
	File? samplefile
	File phenofile
	String? sample_id_header
	String outcome
	Boolean binary_outcome
	String covar_headers
	String exposures
	String? delimiter
	String? missing
	Boolean robust
	Int? memory
	Int? disk
	Int? threads

	call process_phenos {
		input:
			phenofile = phenofile,
			sample_id_header = sample_id_header,
			outcome = outcome,
			covar_headers = covar_headers,
			exposure = exposures,
			delimiter = delimiter,
			missing = missing
	}	

        scatter (i in range(length(genofiles))) {
                call run_interaction {
                        input:
                                genofile = genofiles[i],
                                phenofile = process_phenos.pheno_fmt,
                                covar_file = process_phenos.covar_file,
				outcome = outcome,
                                binary_outcome = binary_outcome,
				missing = missing,
                                robust = robust,
                                memory = memory,        
                                disk = disk
                }
        }

        
        parameter_meta {
                genofiles: "Imputed genotypes in bgen format"
                phenofile: "Tab-delimited phenotype file with Family IDs in the first column, Individual IDs in the second column, missing in the third column and the outcome of interest (quantitative or binary) in the following columns. Phenotyple file with multiple outcomes is supported. Phenotype file used in analysis can be specified by npheno option with the name of phenotype variable"
                npheno: "Speciﬁes which column in the phenotype ﬁle to use for analysis. If value matches any token in the header, the corresponding column is used. Otherwise, value is interpreted as a number, and the (value+3)−th column will be used. The default value is 1, meaning to analyse the phenotype in the 4−th column of the phenotype ﬁle."
                ncovar: "Speciﬁes that a column in the phenotype ﬁle is to be used as a covariate in the analysis. The argument value is handled in the same way as for −−npheno. The −−ncovar option can be used multiple times to specify multiple covariates, but at present only the ﬁrst covariate is used for interaction analysis."
                binary_outcome: "Boolean: is the outcome binary? Otherwise, quantitative is assumed. Please code the binary outcome as 0 & 1, the software set the split point to be 0.5 to code levels of the outcome. "
                interaction: "Boolean: should an interaction term be included? The first covariate in the phenotype file will be used. Defaults to true. If multiple ncovar option exist, only the ﬁrst covariate is used for interaction analysis."
                robust: "Boolean: should robust/sandwich/Huber-White standard errors be used? Default to be false."
                Compute_option: "Optional computational results provided by the software, such as −−compute−alphaHat, computing a simple method-of-moments estimator for alpha, −−compute−MAF, computing the (expected) minor allele frequency for each SNP/locus, −−compute−rSqHat, computing r−squared-hat for each SNP, which is the (estimated) fraction of variance in unobserved 0/1/2 genotype explained by the the individual mean genotypes."
                out_names: "Names to be included for distinguishing output files."
                memory: "Memory required for the modeling step (in GB)."
                disk: "Disk space required for the modeling step (in GB)."
        }

        meta {
                author: "Ye Chen"
                email: "ychen115@mgh.harvard.edu"
                description: "Run gene-environment interaction tests using the QuickTest program."
        }
}
