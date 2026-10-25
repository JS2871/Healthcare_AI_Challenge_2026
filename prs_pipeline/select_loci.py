import sys
import pandas as pd
from statsmodels.stats.multitest import multipletests

PHENO, OUT = sys.argv[1], sys.argv[2]
SETS = ["fold0", "fold1", "fold2", "fold3", "fold4", "all"]
P_THRESHOLDS = {"p1e-5": 1e-5, "p1e-4": 1e-4, "p1e-3": 1e-3}
SPLITHALF_P = 1e-3
METHODS = list(P_THRESHOLDS) + ["fdr05", "splithalf"]
APOE_LOCI = {"rs429358", "rs7412"}


def read_gwas(part, s):
    df = pd.read_csv(f"gwas_{part}_{s}.{PHENO}.glm.linear", sep="\t").dropna(subset=["P"])
    return df.set_index("ID")


summary = []
for s in SETS:
    disc = read_gwas("discovery", s)
    # dedupe correlated LD-tag SNPs (rs1530440_ld1..20) to the lowest-p representative per base locus
    disc["base_id"] = disc.index.str.replace(r"_ld\d+$", "", regex=True)
    # the GWAS excludes rs429358/rs7412 themselves but not their LD tags (rs7412_ld1..): drop those too,
    # otherwise the PRS re-encodes apoe_geno_group
    disc = disc[~disc["base_id"].isin(APOE_LOCI)]
    lead = disc.loc[disc.groupby("base_id")["P"].idxmin()].copy()

    selected = {m: lead["P"] < t for m, t in P_THRESHOLDS.items()}
    selected["fdr05"] = pd.Series(multipletests(lead["P"], alpha=0.05, method="fdr_bh")[0], index=lead.index)
    a, b = read_gwas("halfA", s).reindex(lead.index), read_gwas("halfB", s).reindex(lead.index)
    selected["splithalf"] = (a["P"] < SPLITHALF_P) & (b["P"] < SPLITHALF_P) & (a["BETA"] * b["BETA"] > 0)

    score = lead[["A1"]].copy()
    for m in METHODS:
        score[f"BETA_{m}"] = lead["BETA"].where(selected[m], 0.0)
        sel = lead[selected[m]].sort_values("P")
        summary.append({
            "fold": s, "method": m, "n_selected": len(sel),
            "n_named_loci": int((~sel["base_id"].str.startswith("bg_")).sum()),
            "top_locus": sel["base_id"].iloc[0] if len(sel) else None,
            "top_p": sel["P"].iloc[0] if len(sel) else None,
            "named_loci": ";".join(sel.loc[~sel["base_id"].str.startswith("bg_"), "base_id"]),
        })
    score = score[(score.filter(like="BETA_") != 0).any(axis=1)]
    score.rename_axis("SNP").reset_index().to_csv(f"{OUT}/score_files/score_{s}.txt", sep="\t", index=False)

pd.DataFrame(summary).to_csv(f"{OUT}/selection_summary.csv", index=False)
