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

# compute enrichments for age associated features
uv run phase2/analyses/feature_enrichment.py --modality rna --cell-type-specificity
uv run phase2/analyses/feature_enrichment.py --modality atac --cell-type-specificity
uv run phase2/analyses/feature_enrichment.py --modality rna --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_har.csv --name har
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_har.csv --name har
uv run phase2/analyses/feature_enrichment.py --modality rna --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_haqer.csv --name haqer
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_haqer.csv --name haqer
uv run phase2/analyses/feature_enrichment.py --modality rna --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_telomere.csv --distance --name telomere
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_telomere.csv --distance --name telomere
uv run phase2/analyses/feature_enrichment.py --modality rna --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_centromere.csv --distance --name centromere
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_centromere.csv --distance --name centromere

# for age associated features see if there is correlation between gene ~ atac for chromatin peaks cis-proximal to the gene, using rna covariates
uv run phase2/analyses/cis_correlation.py --covariates specified --covariates-list PCA_0 PCA_1 PCA_2 PCA_3

# generate summary figure for the cis correlation analysis but for just the age associated features in both modalities
uv run phase2/figures/cis_correlation_summary.py

# for age associated features where the cis-proximal atac peaks are correlated with gene expression perform a conditioned analysis of these pairs, using rna covariates
uv run phase2/analyses/cis_conditioned_regression.py --endo-covariates specified --endo-covariates-list PCA_0 PCA_1 PCA_2 PCA_3
# summarize and visualize the cis conditioned regression analysis results
uv run phase2/figures/cis_conditioned_regression_summary.py

# for age associated features where the cis-proximal atac peaks are correlated with gene expression perform a mediation analysis of these pairs
phase2/run_mediation_jobs.sh
uv run python phase2/analyses/post_cis_correlation_mediation.py

# generate latent features per-celltype using cNMF
tmux new -s brainage
CELLTYPES=("Astrocytes" "Endothelial" "ExN BCL11B" "ExN CUX2" "ExN LAMP5" "ExN RELN" "ExN RMST" "ExN RORB" "ExN SEMA3E" "ExN THEMIS" "InN LAMP5" "InN PAX6" "InN PVALB" "InN SST" "InN VIP" "Microglia" "OPCs" "Oligodendrocytes")
for CELLTYPE in "${CELLTYPES[@]}"; do
  # uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality rna --covariates gex_pool ancestry sample_id sex
  # uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality atac --covariates atac_pool ancestry sample_id sex
  uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality rna --workers 48
  uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality atac --workers 48
done
tmux attach-session -t brainage

# review cNMF stability figures and run latent based analysis using the determined K
for CELLTYPE in "${CELLTYPES[@]}"; do
  uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "$CELLTYPE" --k auto --covariates gex_pool ancestry sex
  uv run phase2/analyses/cnmf_latent_regressions.py --modality atac --cell-type "$CELLTYPE" --k auto --covariates atac_pool ancestry sex
done

# combine the cNMF latent regression output and compute FDRs
uv run phase2/analyses/post_cnmf_latent_regressions.py --modality rna
uv run phase2/analyses/post_cnmf_latent_regressions.py --modality atac
