task run_interaction {
	File genofile
	File phenofile
	String npheno
	String ncovar
	Boolean? binary_outcome
	Boolean? interaction
	Boolean? robust
	String? compute_option
	String out_name
	Int? memory = 10
	Int? disk = 20
    Boolean? bgen


	command {
		/quicktest-1.1_bgen_v1.2/quicktest \
		--geno ${genofile} \
		--pheno ${phenofile} \
		--npheno ${npheno} \
		--ncovar ${ncovar} \
		--method-mean \
		${default="" true="--method-binary 0.5" false="" binary_outcome} \
		${default="--method-interaction" true="--method-interaction" false="" interaction} \
		${default="" true="--method-robust" false="" robust} \
		${compute_option} \
		--out quicktest_${out_name}.out \
        ${default="--bgen" true="--bgen" false="" bgen}
	}

	runtime {
		docker: "deadpan39/gxe_quicktest_12112019:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
	}

	output {
        File res = "quicktest_${out_name}.out"
    }
}


workflow run_quicktest {

        Array[File] genofiles
        File phenofile
        String npheno
        String ncovar
        Boolean? binary_outcome
        Boolean? interaction
        Boolean? robust
        String? compute_option
        Array[String] out_names
        Int? memory
        Int? disk
        Boolean? bgen

        
        scatter (i in range(length(genofiles))) {
                call run_interaction {
                        input:
                                genofile = genofiles[i],
                                phenofile = phenofile,
                                npheno = npheno,
                                ncovar = ncovar,
                                binary_outcome = binary_outcome,
                                interaction = interaction,
                                robust = robust,
                                compute_option = compute_option,
                                out_name = out_names[i],
                                memory = memory,        
                                disk = disk,
                                bgen = bgen
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
                description: "Run interaction tests using the Quicktest."
        }
}
