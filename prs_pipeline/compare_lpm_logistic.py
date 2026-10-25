import numpy as np
import pandas as pd

# chr1, fold0 discovery set: linear probability model vs logistic GWAS, per disease
rows = []
for ph in ["HTN", "DM", "LIP"]:
    lo = pd.read_csv(f"chk_logistic.{ph}.glm.logistic.hybrid", sep="\t").set_index("ID")
    li = pd.read_csv(f"chk_lpm.{ph}.glm.linear", sep="\t").set_index("ID").reindex(lo.index)
    ok = lo["P"].notna() & li["P"].notna()
    lo, li = lo[ok], li[ok]
    rows.append({
        "pheno": ph, "n_snps": int(ok.sum()),
        "r_neglog10p": np.corrcoef(-np.log10(li["P"]), -np.log10(lo["P"]))[0, 1],
        "r_beta": np.corrcoef(li["BETA"], lo["BETA"])[0, 1],
        "sign_agreement": (np.sign(li["BETA"]) == np.sign(lo["BETA"])).mean(),
        "top100_overlap": len(set(li["P"].nsmallest(100).index) & set(lo["P"].nsmallest(100).index)),
        "n_p1e-3_lpm": int((li["P"] < 1e-3).sum()), "n_p1e-3_logistic": int((lo["P"] < 1e-3).sum()),
        "top_snp_lpm": li["P"].idxmin(), "top_snp_logistic": lo["P"].idxmin(),
    })
pd.DataFrame(rows).to_csv("lpm_vs_logistic_chr1.csv", index=False)
print(pd.DataFrame(rows).round(4).to_string(index=False))
