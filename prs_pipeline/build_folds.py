import sys
import numpy as np
import pandas as pd

# fold_assignment.csv is copied from 20260927_Step1/prs_pipeline (seed=42, same folds as every earlier Step 1 run)
HERE = sys.argv[1]
folds = pd.read_csv(f"{HERE}/fold_assignment.csv")


def write_keep(ids, path):
    pd.DataFrame({"FID": ids, "IID": ids}).to_csv(path, sep="\t", index=False)


def write_set(disc, name, seed):
    write_keep(disc, f"folds/discovery_{name}.txt")
    # random halves of the discovery set, for split-half replication
    half = np.random.RandomState(seed).permutation(len(disc)) % 2
    write_keep(disc[half == 0], f"folds/halfA_{name}.txt")
    write_keep(disc[half == 1], f"folds/halfB_{name}.txt")


for f in range(5):
    write_set(folds.loc[folds.fold != f, "T_ID"].values, f"fold{f}", 100 + f)
    write_keep(folds.loc[folds.fold == f, "T_ID"].values, f"folds/heldout_fold{f}.txt")

# in-sample comparison: discover on all 10,000 and score the same 10,000
write_set(folds["T_ID"].values, "all", 200)
write_keep(folds["T_ID"].values, "folds/heldout_all.txt")
