import pandas as pd

CLINICAL_DIR = "/mnt/d/Healthcare_AI_Challenge/Healthcare_AI_Challenge_2026/CLINICAL_DATA"
clinical = pd.read_csv(f"{CLINICAL_DIR}/BASE_DATA1.txt")[["T_ID", "T_SEX", "T_AGE"]]
disease = pd.read_csv(f"{CLINICAL_DIR}/BASE_DATA2.txt")[["T_ID", "T_HTN", "T_DM", "T_LIP"]]

covar = clinical.rename(columns={"T_ID": "IID", "T_SEX": "SEX", "T_AGE": "AGE"})
covar.insert(0, "FID", covar["IID"])
covar.to_csv("covar.txt", sep="\t", index=False)

# T_HTN/T_DM/T_LIP: 1=no, 2=yes, 99999=unknown -> NA.
# pheno_logistic.txt keeps plink's 1=control/2=case coding (logistic GWAS, chr1 validation only).
# pheno.txt shifts to 0.5/1.5 so plink2 fits a linear model (it treats {0,1,2} as binary); the
# shift doesn't change the per-SNP beta, which is the risk difference per allele (linear probability model).
pheno = disease.rename(columns={"T_ID": "IID", "T_HTN": "HTN", "T_DM": "DM", "T_LIP": "LIP"})
pheno[["HTN", "DM", "LIP"]] = pheno[["HTN", "DM", "LIP"]].replace(99999, pd.NA).astype("Int64")
pheno.insert(0, "FID", pheno["IID"])
pheno.to_csv("pheno_logistic.txt", sep="\t", index=False, na_rep="NA")
pheno[["HTN", "DM", "LIP"]] = pheno[["HTN", "DM", "LIP"]].astype(float) - 0.5
pheno.to_csv("pheno.txt", sep="\t", index=False, na_rep="NA")
