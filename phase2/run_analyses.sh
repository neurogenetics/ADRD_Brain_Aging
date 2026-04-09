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

# visualize similarities of age effects between tissues
for MODALITY in ${MODALITIES[@]}; do
  uv run phase2/figures/celltype_age_effect_similarity.py --modality ${MODALITY}
done

# for age associated features see if there is correlation between gene ~ atac for chromatin peaks cis-proximal to the gene
uv run phase2/analyses/cis_correlation.py --covariates specified --covariates-list PCA_0_endo PCA_1_endo PCA_2_endo PCA_3_endo PCA_0_exog PCA_1_exog PCA_2_exog PCA_3_exog

# generate summary figure for the cis correlation analysis but for just the age associated features in both modalities
uv run phase2/figures/cis_correlation_summary.py

# for age associated features where the cis-proximal atac peaks are correlated with gene expression perform a mediation analysis of these pairs
tmux new -s brainage
uv run phase2/analyses/cis_correlation_mediation.py --endo-covariates specified --endo-covariates-list PCA_0_endo PCA_1_endo PCA_2_endo PCA_3_endo --exog-covariates specified --exog-covariates-list PCA_0_exog PCA_1_exog PCA_2_exog PCA_3_exog --debug
tmux attach-session -t brainage

# generate latent features per-celltype using cNMF
tmux new -s brainage
uv run phase2/analyses/cnmf_latent_generation.py --modality rna --covariates sample_id sex gex_pool --components 4 5 6 7 8 9 10 11 12 13 14 15 16
tmux attach-session -t brainage

# review cNMF stability figures and run latent based analysis using the determined K
# Astrocytes, K=11
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type Astrocytes --k 11
# Endothelial, K=14
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type Endothelial --k 14
# ExN_BCL11B, K=7
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "ExN BCL11B" --k 7
# ExN_CUX2, K=6
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "ExN CUX2" --k 6
# ExN_LAMP5, K=8
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "ExN LAMP5" --k 8
# ExN_RELN, K=5
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "ExN RELN" --k 5
# ExN_RMST, K=10
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "ExN RMST" --k 10
# ExN_RORB, K=8
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "ExN RORB" --k 8
# ExN_SEMA3E, K=6
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "ExN SEMA3E" --k 6
# ExN_THEMIS, K=6
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "ExN THEMIS" --k 6
# InN_LAMP5, K=10
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "InN LAMP5" --k 10
# InN_PAX6, K=7
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "InN PAX6" --k 7
# InN_PVALB, K=6
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "InN PVALB" --k 6
# InN_SST, K=8
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "InN SST" --k 8
# InN_VIP, K=10
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "InN VIP" --k 10
# Microglia, K=12
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type Microglia --k 12
# Microglia, K=12
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type Microglia --k 12
# Oligodendrocytes, K=7
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type Oligodendrocytes --k 7
# OPCs, K=10
uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type OPCs --k 10
