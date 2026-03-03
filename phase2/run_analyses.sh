#!/bin/bash

echo "$(date): $(uname -n)"
echo "Executing: $0 $@"

# run the per cell-type pseudobulk data prep
phase2/run_prep_pb_jobs.sh prep_pb_data.py rna
phase2/run_prep_pb_jobs.sh prep_pb_data.py atac

# run the variance partition analyis per cell-type
phase2/run_prep_pb_jobs.sh run_variance_partition.py rna
phase2/run_prep_pb_jobs.sh run_variance_partition.py atac

# run the age regression analysis per cell-type
phase2/run_regression_jobs.sh ols rna
phase2/run_regression_jobs.sh ols atac
phase2/run_regression_jobs.sh rlm rna
phase2/run_regression_jobs.sh rlm atac
