import sys
import pandas as pd

PHENO, OUT = sys.argv[1], sys.argv[2]
prefix = f"PRS_{PHENO}"


def read_scores(s):
    sc = pd.read_csv(f"prs_{PHENO}_{s}.sscore", sep="\t").set_index("IID")
    cols = [c for c in sc.columns if c.startswith("BETA_") and c.endswith("_SUM")]
    return sc[cols].rename(columns=lambda c: f"{prefix}_{c[len('BETA_'):-len('_SUM')]}")


for name, sets in [("prs_all_methods.csv", [f"fold{f}" for f in range(5)]), ("prs_insample_all_methods.csv", ["all"])]:
    prs = pd.concat([read_scores(s) for s in sets]).rename_axis("T_ID").reset_index()
    assert len(prs) == 10000 and prs["T_ID"].is_unique
    prs.to_csv(f"{OUT}/{name}", index=False)
