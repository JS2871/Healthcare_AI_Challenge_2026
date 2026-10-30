# Multi-Omics-Based Classification of Dementia Risk

This repository contains the code developed for the **KNIH–KSBi Healthcare AI Challenge 2026** (국립보건연구원–한국생명정보학회 공동주관), under the *Cohort-Based Healthcare AI Model Development and Validation* (코호트 기반 헬스케어 인공지능 모델 개발·검증) track.

## Overview

Early identification of people at high risk of dementia may delay progression and extend independent living,
but assessing that risk early is still hard in clinical practice. This project builds a multimodal model called **K-BAP (Korean Brain Atrophy prediction)** that combines **genomic** (SNP array) and **clinical-epidemiological** data. It classifies participants by the
severity of **brain atrophy**, which is measured from brain-MRI ROI volumes.

Brain imaging is used **only to define the targets**. The predictors are genomic and clinical features.

## Model: K-BAP (Korean Brain Atrophy Prediction)

### K-BAP Objectives

**Objective 1: Quantification of Brain Atrophy Severity**
Integrate atrophy measurements across multiple brain regions into a unified composite score to quantify multidimensional brain atrophy severity and define high-risk groups based on this score.

**Objective 2: Prediction of High-Risk Brain Atrophy**
Develop predictive models using genomic, clinical, and epidemiological data to identify individuals at high risk of brain atrophy prior to MRI examination.


## Pipeline

| Step | File | What it does |
|---|---|---|
| 1 | `Step1_Targets.ipynb` | Builds 5 atrophy-severity targets (composite, hippocampus, whole brain, ventricles, amygdala) from ROI volumes: ICV adjustment (OLS; two-limit Tobit for the double-censored ventricles), sign correction, and a PCA composite. These are turned into binary top-20/10/5% labels and ordered `risk5`/`risk4` classes. |
| 1b | `prs_pipeline/build_prs_pipeline.sh` | Cardiometabolic PRS (hypertension, diabetes, dyslipidemia) via 5-fold **out-of-fold** GWAS (plink2, linear probability model, age+sex covariates, APOE SNPs excluded), with five variant-selection rules (p<1e-5/1e-4/1e-3, BH-FDR<0.05, split-half replication). It also extracts APOE (rs429358/rs7412) dosages. |
| 2 | `Step2_Features.ipynb` | Assembles the predictor matrix: curated clinical variables (with sentinel-code handling and disease durations), the APOE genotype group, and the out-of-fold PRS. The PRS selection method is chosen on the training split only. |
| 3 | `Step3_Modeling.ipynb` | Models: an 80/20 split, then a preprocessing pipeline with an L1 feature selector and `RandomizedSearchCV`. Regression, binary classification (LR, linear/RBF SVM, RF, LightGBM; Youden threshold chosen on training CV) and ordinal multiclass (including a custom proportional-odds logistic estimator). It also covers permutation importance, odds ratios, ROC curves and confusion matrices. |

### `prs_pipeline/`
- `build_prs_pipeline.sh`: driver (merge per-chromosome PLINK files → folds → GWAS → locus selection → scoring)
- `build_folds.py`: discovery/held-out keep files and split-half sets per fold
- `build_covar_pheno.py`: plink covariate and phenotype files
- `compare_lpm_logistic.py`: chr1 check that compares the linear-probability and logistic GWAS
- `select_loci.py`: LD-tag deduplication, APOE exclusion, the five selection rules, and score files
- `concat_prs.py`: concatenates per-fold scores into one PRS table

## Requirements

- Python 3 with `pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-learn` (0.24-compatible code), `lightgbm`,
  `xgboost`, `matplotlib`, `seaborn`, `joblib`, `openpyxl`
- [`plink2`](https://www.cog-genomics.org/plink/2.0/) for the PRS pipeline

## Running

1. Edit the path constants at the top of each script or notebook to point to your local copy of the challenge data and to an
   output directory: `CLINICAL_DIR`, `SNP`/`SNP_ARRAY_DIR`, `STEP_DIR`, `TABLE_DIR`, `FIG_DIR`, `NOTEBOOK_DIR`,
   `SUPPLE_FIG_DIR`, `SUPPLE_TABLE_DIR`.
2. Run `Step1_Targets.ipynb`.
3. Run `bash prs_pipeline/build_prs_pipeline.sh [work_dir]` (this needs plink2 and takes several hours).
4. Run `Step2_Features.ipynb`, then `Step3_Modeling.ipynb`.

## Data availability and omitted content

The data are synthetic (K-CHIP SNP array, baseline/follow-up survey data, and brain-MRI ROI volumes). They were
provided **only to challenge participants** under the competition's synthetic-data use and security agreement
(합성데이터 이용 및 보안 서약서). That agreement prohibits redistributing the data and publicly releasing outputs or
derived results. For that reason, this repository contains **code only**:

- no input data, and no participant-level derived files (PRS tables, genotype extracts, fold assignments)
- no selected-variant lists or GWAS summary results
- notebook outputs are cleared, and result figures, tables and numeric findings are removed from the text

If you want to run the code, you need your own authorized access to the challenge data.
