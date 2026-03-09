#!/bin/bash

echo "$(date): $(uname -n)"
echo "Executing: $0 $@"

REGRESSIONTYPE=$1 #"wls" "ols" "rlm"
MODALITY=$2       #"rna" "atac"

CELLTYPES="Astrocytes Endothelial ExN_BCL11B ExN_CUX2 ExN_LAMP5 ExN_RELN ExN_RMST ExN_RORB ExN_SEMA3E ExN_THEMIS InN_LAMP5 InN_PAX6 InN_PVALB InN_SST InN_VIP Microglia OPCs Oligodendrocytes"

for CELLTYPE in ${CELLTYPES[@]}; do
  echo "Running: uv run phase2/analyses/pseudobulk_regression.py --modality ${MODALITY} --cell-type ${CELLTYPE} --regression-type ${REGRESSIONTYPE}"
  uv run phase2/analyses/pseudobulk_regression.py --modality ${MODALITY} --cell-type ${CELLTYPE} --regression-type ${REGRESSIONTYPE} --covariates specified --covariates-list PCA_0 PCA_1 PCA_2 PCA_3 &
done

wait
