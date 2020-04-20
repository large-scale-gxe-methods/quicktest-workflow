import sys
import pandas as pd


phenofile, sample_id_header, outcome, exposure, covar_names, delimiter, missing, samplefile = sys.argv[1:9]

phenos = pd.read_csv(phenofile, sep=delimiter, na_values=missing)

if samplefile != "":  # Optional sample file ensures correct ordering of IDs
    sample = pd.read_csv(samplefile, sep=" ", skiprows=2, header=None, usecols=[1],
                         names=["id"])
    phenos = pd.merge(sample, phenos, how="left", left_on="id",
                      right_on=sample_id_header)


covars = [] if covar_names == "" else covar_names.split(" ")
output_cols = ["id", "id", outcome, exposure] + covars
phenos = phenos.loc[:, output_cols]
phenos.columns.values[:2] = ["FID", "IID"]

phenos.to_csv("quicktest_phenotypes.csv", sep=" ", index=False, na_rep="NA")

covar_string = "--ncovar " + exposure
for covar in covars:
    covar_string += " --ncovar " + covar
with open("covar_string.txt", "w") as textfile:
    textfile.write(covar_string)
