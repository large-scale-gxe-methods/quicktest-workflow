import sys
import pandas as pd
import numpy as np
from scipy.stats import chi2


resfile, outfile, robust = sys.argv[1:4]

names_dict = {'id1': 'SNPID', 'AlleleA': 'Allele1', 'AlleleB': 'Allele2', 
              'snp.beta': 'Beta_Main', 'snp.se': 'SE_Beta_Main', 
              'interaction.beta': 'Beta_Interaction', 
              'interaction.se': 'SE_Beta_Interaction',
              'cov.snp.interaction': 'Cov_Main_Interaction'}

if robust == "false":
    pd.read_csv(resfile, sep=" ").to_csv(outfile, sep=" ", index=False,
                                         na_rep="NaN")
    sys.exit(0)

def calc_joint_test(beta_main, var_beta_main, beta_int, var_beta_int, cov_main_int, df=1):
     # Given a set of summary statistics from interaction testing, calculate
     # a p-value for the joint test of genetic effect (main + interaction)
     det_sigma = var_beta_main * var_beta_int - cov_main_int**2  # Determinant of covariance matrix
     T_stat = (beta_main**2 * var_beta_int -  # Full expression for 2-DF Chi-squared statistic
         2 * beta_main * beta_int * cov_main_int + beta_int**2 * var_beta_main) / det_sigma
     T_stat[det_sigma <= 0] = np.nan  # In case of possible numerical rounding errors
     return 1 - chi2.cdf(T_stat, df=df + 1)

res = (pd.read_csv(resfile, sep=" ")
       .rename(columns=names_dict)
       .filter(list(names_dict.values()))
       .assign(Var_Beta_Main = lambda x: x.SE_Beta_Main**2,
               Var_Beta_Interaction = lambda x: x.SE_Beta_Interaction**2)
       .assign(P_Value_Main = lambda x: 1 - chi2.cdf(x.Beta_Main ** 2 / x.Var_Beta_Main, df=1),
               P_Value_Interaction = lambda x: 1 - chi2.cdf(x.Beta_Interaction**2 / x.Var_Beta_Interaction, df=1),
               P_Value_Joint = lambda x: calc_joint_test(x.Beta_Main, x.Var_Beta_Main, x.Beta_Interaction, x.Var_Beta_Interaction,
                                                         x.Cov_Main_Interaction))
       .drop(["SE_Beta_Main", "SE_Beta_Interaction"], axis="columns"))
res.to_csv(outfile, sep=" ", index=False, na_rep="NaN")
