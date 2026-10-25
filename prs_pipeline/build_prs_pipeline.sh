#!/bin/bash
# 20260930 revision of 20260927_Step1/prs_pipeline/build_prs_pipeline.sh: PRS for the three
# cardiometabolic diagnoses (T_HTN, T_DM, T_LIP) instead of for the atrophy targets. The atrophy
# targets are never used here, so no PRS can leak the target.
#
# Two schemes per disease:
#   out-of-fold (used as the model feature) - 5-fold, same fold assignment as 20260927; each fold's
#       held-out individuals are scored with variants/weights from the other 4 folds only.
#   in-sample (comparison only) - GWAS on all 10,000, score the same 10,000. Its PRS partly re-encodes
#       each person's own observed diagnosis (overfit), see 20260930_Step1_Features.ipynb.
# Selection methods (within each discovery set): p1e-5, p1e-4, p1e-3, fdr05, splithalf - as in 20260927.
#
# GWAS model: linear probability model (disease coded 0.5/1.5 -> plink2 --glm linear), age + sex
# covariates, APOE SNPs (rs429358, rs7412) excluded from the GWAS and their LD tags excluded in select_loci.py. Logistic --glm is ~100x slower here (~13h for the full pipeline);
# step 3 checks on chr1 that the linear model gives the same SNP ranking/effect signs as logistic.
#
# Outputs, one directory per disease (../genomics_htn, ../genomics_dm, ../genomics_lip):
#   prs_all_methods.csv           - out-of-fold PRS per T_ID for every method (10,000 rows)
#   prs_insample_all_methods.csv  - in-sample PRS per T_ID for every method (10,000 rows)
#   selection_summary.csv         - per fold (fold0-4, all) x method: n loci selected, n named loci, top locus
#   score_files/score_*.txt       - SNP/A1/beta per method (0 = not selected)
# Shared output here: lpm_vs_logistic_chr1.csv. Requires plink2 and the merged combined.pgen
# (built by the 20260927 pipeline in /tmp/step1_prs_work_20260927, rebuilt below if missing).
set -euo pipefail
SNP=/mnt/d/Healthcare_AI_Challenge/Healthcare_AI_Challenge_2026/SNP_ARRAY
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${1:-/tmp/step1_prs_work_20260930}
mkdir -p "$WORK/folds"
cd "$WORK"

# 1. Genome-wide pfile (reuse 20260927's merge if present).
if [ ! -f combined.pgen ]; then
  if [ -f /tmp/step1_prs_work_20260927/combined.pgen ]; then
    for ext in pgen pvar psam; do ln -sf /tmp/step1_prs_work_20260927/combined.$ext .; done
  else
    rm -f merge_list.txt
    for i in $(seq 1 22); do echo "$SNP/synthetic_kchip_chr${i}" >> merge_list.txt; done
    plink2 --pmerge-list merge_list.txt bfile --make-pgen --out combined --threads 8
  fi
fi
echo -e "rs429358\nrs7412" > apoe_snps.txt

# 2. Folds (from fold_assignment.csv) + in-sample set, covariates, phenotypes.
python3 "$HERE/build_folds.py" "$HERE"
python3 "$HERE/build_covar_pheno.py"

# 3. Linear-probability vs logistic check (chr1, fold0 discovery set).
plink2 --pfile combined --keep folds/discovery_fold0.txt --exclude apoe_snps.txt --chr 1 \
  --pheno pheno_logistic.txt --pheno-name HTN,DM,LIP --covar covar.txt --covar-name AGE,SEX \
  --glm hide-covar firth-fallback cols=+beta,+a1freq --out chk_logistic --threads 8
plink2 --pfile combined --keep folds/discovery_fold0.txt --exclude apoe_snps.txt --chr 1 \
  --pheno pheno.txt --pheno-name HTN,DM,LIP --covar covar.txt --covar-name AGE,SEX \
  --glm hide-covar cols=+beta,+a1freq --out chk_lpm --threads 8
python3 "$HERE/compare_lpm_logistic.py"
cp lpm_vs_logistic_chr1.csv "$HERE/"

# 4. GWAS per set on the full discovery set and on each half (all three diseases in one pass).
for s in fold0 fold1 fold2 fold3 fold4 all; do
  for part in discovery halfA halfB; do
    plink2 --pfile combined --keep folds/${part}_${s}.txt --exclude apoe_snps.txt \
      --pheno pheno.txt --pheno-name HTN,DM,LIP --covar covar.txt --covar-name AGE,SEX \
      --glm hide-covar cols=+beta,+a1freq --out gwas_${part}_${s} --threads 8
  done
done

# 5-7. Per disease: select loci, score held-out sets, concatenate.
for pair in HTN:genomics_htn DM:genomics_dm LIP:genomics_lip; do
  PHENO=${pair%%:*}
  OUT="$HERE/../${pair##*:}"
  mkdir -p "$OUT/score_files"
  python3 "$HERE/select_loci.py" "$PHENO" "$OUT"
  for s in fold0 fold1 fold2 fold3 fold4 all; do
    plink2 --pfile combined --keep folds/heldout_${s}.txt \
      --score "$OUT/score_files/score_${s}.txt" 1 2 header-read cols=+scoresums \
      --score-col-nums 3-7 --out prs_${PHENO}_${s} --threads 8
  done
  python3 "$HERE/concat_prs.py" "$PHENO" "$OUT"
done
