import sys
import pandas as pd
from scipy.stats import chi2


resfile, outfile, robust = sys.argv[1:4]

names_dict = {'id1': 'SNPID', 'AlleleA': 'Allele1', 'AlleleB': 'Allele2', 
	      'snp.beta': 'Beta_Main', 'snp.se': 'SE_Beta_Main', 
	      'interaction.beta': 'Beta_Interaction_1', 
	      'interaction.se': 'SE_Beta_Interaction_1_1'}

if robust == "false":
    pd.read_csv(resfile, sep=" ").to_csv(outfile, sep=" ", index=False,
                                         na_rep="NaN")
    sys.exit(0)

res = (pd.read_csv(resfile, sep=" ")
       .rename(columns=names_dict)
       .filter(list(names_dict.values()))
       .assign(Var_Beta_Main = lambda x: x.SE_Beta_Main ** 2,
	       Var_Beta_Interaction_1_1 = lambda x: x.SE_Beta_Interaction_1_1 ** 2)
       .assign(P_Value_Main = lambda x: 1 - chi2.cdf(x.Beta_Main ** 2 / x.Var_Beta_Main, df=1),
	       P_Value_Interaction = lambda x: 1 - chi2.cdf(x.Beta_Interaction_1 ** 2 / x.Var_Beta_Interaction_1_1, df=1),
	       P_Value_Joint = lambda x: 1 - chi2.cdf(x.Beta_Main ** 2 / x.Var_Beta_Main + x.Beta_Interaction_1 ** 2 / x.Var_Beta_Interaction_1_1, df=2))
       .drop(["SE_Beta_Main", "SE_Beta_Interaction_1_1"], axis="columns"))
res.to_csv(outfile, sep=" ", index=False, na_rep="NaN")
