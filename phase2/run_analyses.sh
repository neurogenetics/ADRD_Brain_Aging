#!/bin/bash

echo "$(date): $(uname -n)"
echo "Executing: $0 $@"

# format the initial knowns covariates by modality
uv run phase2/analysis/format_covariates.py

# pseudobulk convert the modalities and save per cell-type
uv run phase2/analyses/pseudobulk_convert.py --aggregate-type sum

# run the per cell-type pseudobulk data prep and generate non-target variance components
MODALITIES="rna atac"
for MODALITY in ${MODALITIES[@]}; do
  phase2/run_prep_pb_jobs.sh prep_pb_data.py ${MODALITY}
done

# format the covariates to include the non-target variance components

# run the variance partition analyis per cell-type
for MODALITY in ${MODALITIES[@]}; do
  phase2/run_prep_pb_jobs.sh run_variance_partition.py ${MODALITY}
done

# run the age regression analysis per cell-type
REGRESSTYPES="wls vwrlm ols rlm"
for MODALITY in ${MODALITIES[@]}; do
  for REGRESSTYPE in ${REGRESSTYPES[@]}; do
    phase2/run_regression_jobs.sh ${REGRESSTYPE} ${MODALITY}
  done
done

# post process the age regressions
for MODALITY in ${MODALITIES[@]}; do
  for REGRESSTYPE in ${REGRESSTYPES[@]}; do
    uv run phase2/analyses/post_pseudobulk_regression.py --modality ${MODALITY} --regression-type ${REGRESSTYPE} --no-volcano-per-celltype
  done
done

# check the general regression against the robust to check for outlier driven results
for MODALITY in ${MODALITIES[@]}; do
  uv run phase2/analyses/filter_regression_type_differences.py --modality ${MODALITY}
done

# for age associated features see if there is correlation between gene ~ atac for chromatin peaks cis-proximal to the gene
uv run phase2/analyses/cis_correlation.py --covariates specified --covariates-list PCA_0_endo PCA_1_endo PCA_2_endo PCA_3_endo PCA_0_exog PCA_1_exog PCA_2_exog PCA_3_exog

# for age associated features where the cis-proximal atac peaks are correlated with gene expression perform a mediation analysis of these pairs
nohup uv run phase2/analyses/run_mediation.py --endo-covariates specified --endo-covariates-list PCA_0_endo PCA_1_endo PCA_2_endo PCA_3_endo --exog-covariates specified --exog-covariates-list PCA_0_exog PCA_1_exog PCA_2_exog PCA_3_exog --debug &

# generate latent features per-celltype using cNMF
uv run phase2/analyses/run_cnmf.py --modality rna --covariates sample_id sex gex_pool --components 4 5 6 7 8 9 10 11 12 13 14 15 16
