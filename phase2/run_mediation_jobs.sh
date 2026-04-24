#!/bin/bash

echo "$(date): $(uname -n)"
echo "Executing: $0 $@"

CELLTYPES="Astrocytes Endothelial ExN_BCL11B ExN_CUX2 ExN_LAMP5 ExN_RELN ExN_RMST ExN_RORB ExN_SEMA3E ExN_THEMIS InN_LAMP5 InN_PAX6 InN_PVALB InN_SST InN_VIP Microglia OPCs Oligodendrocytes"

for CELLTYPE in ${CELLTYPES[@]}; do
  echo "Running: phase2/analyses/cis_correlation_mediation.py --cell-type ${CELLTYPE} "
  uv run phase2/analyses/cis_correlation_mediation.py --cell-type ${CELLTYPE} --endo-covariates specified --endo-covariates-list PCA_0_endo PCA_1_endo PCA_2_endo PCA_3_endo --exog-covariates specified --exog-covariates-list PCA_0_exog PCA_1_exog PCA_2_exog PCA_3_exog
done
