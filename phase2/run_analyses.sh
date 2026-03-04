#!/bin/bash

echo "$(date): $(uname -n)"
echo "Executing: $0 $@"

# format the initial knowns covariates by modality
uv run phase2/analysis/format_covariates.py

# pseudobulk convert the modalities and save per cell-type
uv run phase2/analyses/pseudobulk_convert.py

# run the per cell-type pseudobulk data prep and generate non-target variance components
phase2/run_prep_pb_jobs.sh prep_pb_data.py rna
phase2/run_prep_pb_jobs.sh prep_pb_data.py atac

# format the covariates to include the non-target variance components

# run the variance partition analyis per cell-type
phase2/run_prep_pb_jobs.sh run_variance_partition.py rna
phase2/run_prep_pb_jobs.sh run_variance_partition.py atac

# run the age regression analysis per cell-type
phase2/run_regression_jobs.sh ols rna
phase2/run_regression_jobs.sh ols atac
phase2/run_regression_jobs.sh rlm rna
phase2/run_regression_jobs.sh rlm atac

# post process the age regressions
uv run phase2/analyses/post_pseudobulk_regression.py --modality rna --regression-type ols --no-volcano-per-celltype
uv run phase2/analyses/post_pseudobulk_regression.py --modality rna --regression-type rlm --no-volcano-per-celltype
