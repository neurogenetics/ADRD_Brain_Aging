#!/bin/bash

echo "$(date): $(uname -n)"
echo "Executing: $0 $@"

# format the initial knowns covariates by modality
uv run phase2/analyses/format_covariates.py --exclude-ids Aging134

# pseudobulk convert the modalities and save per cell-type
uv run phase2/analyses/pseudobulk_convert.py --aggregate-type sum --exclude-ids Aging134

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

# since using WLS check cell-types and modalities for correlations between cell counts and age
uv run phase2/analyses/cell_counts_regression.py

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
    uv run phase2/analyses/post_pseudobulk_regression.py --modality ${MODALITY} --regression-type ${REGRESSTYPE}
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

# power analysis of this study
# min sample sizes used per modality
uv run phase2/analyses/regression_power_analysis.py --sizes 8,15
# max sample sized used per modality
uv run phase2/analyses/regression_power_analysis.py --sizes 35,34

# compute enrichments for age associated features
uv run phase2/analyses/feature_enrichment.py --modality rna --cell-type-specificity --name cell_specific
uv run phase2/analyses/feature_enrichment.py --modality atac --cell-type-specificity --name cell_specific
uv run phase2/analyses/feature_enrichment.py --modality rna --metrics-file /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_gene_metrics.csv --name gene_metrics
uv run phase2/analyses/feature_enrichment.py --modality rna --metrics-file /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/Gene_hallmarks_clean.csv --name gene_aging_hallmarks
uv run phase2/analyses/feature_enrichment.py --modality rna --metrics-file /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/Gene_hallmarks_clean1.csv --name gene_aging_hallmarks_level1
uv run phase2/analyses/feature_enrichment.py --modality rna --metrics-file /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/Gene_hallmarks_clean2.csv --name gene_aging_hallmarks_level2
uv run phase2/analyses/feature_enrichment.py --modality rna --metrics-file /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/Gene_hallmarks_clean3.csv --name gene_aging_hallmarks_level3
uv run phase2/analyses/feature_enrichment.py --modality rna --metrics-file /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/Gene_hallmarks_clean4.csv --name gene_aging_hallmarks_level4
uv run phase2/analyses/feature_enrichment.py --modality rna --metrics-file /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/Gene_hallmarks_clean5.csv --name gene_aging_hallmarks_level5
uv run phase2/analyses/feature_enrichment.py --modality atac --metrics-file /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/quants/aging_phase2_peak_attributes.csv --name peak_metrics
uv run phase2/analyses/feature_enrichment.py --modality rna --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_har.csv --name har &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_har.csv --name har &
uv run phase2/analyses/feature_enrichment.py --modality rna --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_haqer.csv --name haqer &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_haqer.csv --name haqer &
uv run phase2/analyses/feature_enrichment.py --modality rna --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_telomere.csv --distance --name telomere &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_telomere.csv --distance --name telomere &
uv run phase2/analyses/feature_enrichment.py --modality rna --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_centromere.csv --distance --name centromere &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_centromere.csv --distance --name centromere &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_Encode4_cCRE.Proximal_enhancer.csv --name Encode4_Proximal_enhancer &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_Encode4_cCRE.CA-CTCF.csv --name Encode4_CA-CTCF &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_Encode4_cCRE.CA-TF.csv --name Encode4_CA-TF &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_Encode4_cCRE.CA.csv --name Encode4_CA &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_Encode4_cCRE.Distal_enhancer.csv --name Encode4_Distal_enhancer &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_Encode4_cCRE.TF.csv --name Encode4_TF &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_Encode4_cCRE.CA-H3K4me3.csv --name Encode4_CA-H3K4me34 &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_Encode4_cCRE.Promoter.csv --name Encode4_Promoter &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_TEs.csv --name TEs &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/DNAm/clocksites_annotation.csv --name clocksites &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/DNAm/clocksites_annotation_Horvath2013.csv --name clocksites_Horvath2013 &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/DNAm/clocksites_annotation_Shireby2020.csv --name clocksites_Shireby2020 &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/DNAm/clocksites_annotation_Tong2024_BrainClock.csv --name clocksites_Tong2024_BrainClock &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/DNAm/clocksites_annotation_Tong2024_Glia-In.csv --name clocksites_Tong2024_Glia-In &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/DNAm/clocksites_annotation_Tong2024_Glia-Sin.csv --name clocksites_Tong2024_Glia-Sin &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/DNAm/clocksites_annotation_Tong2024_Neu-In.csv --name clocksites_Tong2024_Neu-In &
uv run phase2/analyses/feature_enrichment.py --modality atac --annotation-csv /mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/DNAm/clocksites_annotation_Tong2024_Neu-Sin.csv --name clocksites_Tong2024_Neu-Sin &

# summarize the enrichment results and visualize
uv run phase2/figures/summarize_metrics_enrichments.py
uv run phase2/figures/summarize_all_enrichments.py
uv run phase2/figures/generate_enrichment_plots.py

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
  uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality rna --exclude-ids Aging134 --workers 48
  uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality atac --exclude-ids Aging134 --workers 48
done
# for CELLTYPE in "${CELLTYPES[@]}"; do
# uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality rna --cnmf-dir-name cnmf_harmony --covariates gex_pool ancestry sample_id sex
# uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality atac --cnmf-dir-name cnmf_harmony --covariates atac_pool ancestry sample_id sex
# uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality rna --use-aaf --cnmf-dir-name cnmf_aaf --workers 48
# uv run phase2/analyses/cnmf_latent_generation.py --cell-type "$CELLTYPE" --modality atac --use-aaf --cnmf-dir-name cnmf_aaf --workers 48
# done
tmux attach-session -t brainage

# development test/comparison
uv run phase2/development/analyses/latent_generation.py --cell-type "ExN SEMA3E" --modality rna

# review cNMF stability figures and run latent based analysis using the determined K
for CELLTYPE in "${CELLTYPES[@]}"; do
  uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "$CELLTYPE" --k auto --covariates PCA_0 PCA_1 PCA_2 PCA_3
  uv run phase2/analyses/cnmf_latent_regressions.py --modality atac --cell-type "$CELLTYPE" --k auto --covariates PCA_0 PCA_1 PCA_2 PCA_3
done
# for CELLTYPE in "${CELLTYPES[@]}"; do
#   # uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "$CELLTYPE" --k auto
#   # uv run phase2/analyses/cnmf_latent_regressions.py --modality atac --cell-type "$CELLTYPE" --k auto
#   # uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "$CELLTYPE" --k auto --covariates gex_pool ancestry sex
#   # uv run phase2/analyses/cnmf_latent_regressions.py --modality atac --cell-type "$CELLTYPE" --k auto --covariates atac_pool ancestry sex
#   # uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "$CELLTYPE" --k auto --cnmf-dir-name cnmf_aaf
#   # uv run phase2/analyses/cnmf_latent_regressions.py --modality atac --cell-type "$CELLTYPE" --k auto --cnmf-dir-name cnmf_aaf
#   # uv run phase2/analyses/cnmf_latent_regressions.py --modality rna --cell-type "$CELLTYPE" --k auto --cnmf-dir-name cnmf_aaf --covariates gex_pool ancestry sex
#   # uv run phase2/analyses/cnmf_latent_regressions.py --modality atac --cell-type "$CELLTYPE" --k auto --cnmf-dir-name cnmf_aaf --covariates atac_pool ancestry sex
# done

# combine the cNMF latent regression output and compute FDRs
uv run phase2/analyses/post_cnmf_latent_regressions.py --modality rna
uv run phase2/analyses/post_cnmf_latent_regressions.py --modality atac
# uv run phase2/analyses/post_cnmf_latent_regressions.py --modality rna --cnmf-dir-name cnmf_aaf
# uv run phase2/analyses/post_cnmf_latent_regressions.py --modality atac --cnmf-dir-name cnmf_aaf

# compare the age association latent factor with others across cell-types and modalities
uv run phase2/analyses/compare_age_latent_factors.py

# network analyses using both Quantum walk and random walk
for CELLTYPE in "${CELLTYPES[@]}"; do
  uv run phase2/analyses/run_qml_pipeline.py --cell-type "$CELLTYPE" --modality rna --epochs 100
done
# for CELLTYPE in "${CELLTYPES[@]}"; do
#   for MODALITY in ${MODALITIES[@]}; do
#     uv run phase2/analyses/run_qml_pipeline.py --cell-type "$CELLTYPE" --modality "$MODALITY" --epochs 100
#   done
# done

# network analysis comparisons across cell-types
uv run phase2/analyses/summarize_qml_results.py --help

# network analysis visualizations
for CELLTYPE in "${CELLTYPES[@]}"; do
  uv run phase2/figures/export_walk_topologies.py --cell-type "$CELLTYPE" --modality rna
done
# for CELLTYPE in "${CELLTYPES[@]}"; do
#   for MODALITY in ${MODALITIES[@]}; do
#     uv run phase2/figures/export_walk_topologies.py --cell-type "$CELLTYPE" --modality "$MODALITY"
#   done
# done
