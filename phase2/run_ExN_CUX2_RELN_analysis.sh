# Analyze ExN CUX2 RELN in the HBCC EC data for age and SEA-AD EC data for AD

DATADIR="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"

# since the SAHA analysis we have for the SEAAD EC data has the more granular cell-types already included can use existing utilites
# since the many of the utilities will run for all cell-types at granularity speciefied, which we don't need, just the one subtype
# subset the mudata object to the ExN CUX2 and then run utitlities specifying anno instead of anno_coarse
uv run python -c "
import mudata as md
data_dir = '/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2'
mdata = md.read_h5mu(f'{data_dir}/public/seaad/seaad_ec_multiome_labeled.h5mu')
cux2_reln = mdata[mdata.obs.anno == 'ExN CUX2 RELN']
print(cux2_reln)
cux2_reln.obs.anno = cux2_reln.obs.anno.cat.remove_unused_categories()
cux2_reln.write(f'{data_dir}/public/seaad/seaad_ec_multiome_exn_cux2_reln.h5mu')
"

# for the HBCC EC data will need to subset previous ExN CUX2 cells import the more granular cell-type labels
# to include ExN CUX2 RELN, etc. SAHA labels for the RNA cells and then train a classifier to predict the more
# granular ExN CUX2 subtype labels for the ATAC cells

# convert to pseudobulk per cell-type
uv run phase2/SEAAD/analyses/pseudobulk_convert.py \
  --input-file "$DATADIR"/public/seaad/seaad_ec_multiome_exn_cux2_reln.h5mu \
  --project seaad_ec_multiome \
  --work-dir "$DATADIR"/public/seaad \
  --cell-type-col anno \
  --sample-col donor_id \
  --aggregate-type sum \
  --exclude-ids H21.33.018,H21.33.001,H20.33.036 \
  --debug

# format the initial knowns covariates by modality
uv run phase2/SEAAD/analyses/format_covariates.py \
  --input-file "$DATADIR"/public/seaad/seaad_ec_multiome_exn_cux2_reln.h5mu \
  --cell-type-col anno \
  --sample-col donor_id \
  --project seaad_ec_multiome \
  --output-prefix seaad_ec_multiome.granular \
  --work-dir "$DATADIR"/public/seaad \
  --exclude-ids H21.33.018,H21.33.001,H20.33.036 \
  --debug

# run the per cell-type pseudobulk data prep and generate non-target variance components
CELLTYPE="ExN_CUX2_RELN"
MODALITIES="rna atac"
for MODALITY in ${MODALITIES[@]}; do
  uv run phase2/SEAAD/analyses/prep_pb_data.py \
    --modality ${MODALITY} \
    --cell-type ${CELLTYPE} \
    --covariates-file "$DATADIR"/public/seaad/sample_info/seaad_ec_multiome.granular.covariates.${MODALITY}.csv \
    --features-file "$DATADIR"/quants/aging_phase2.features.csv &
done

# run the variance partition analyis per cell-type
for MODALITY in ${MODALITIES[@]}; do
  uv run phase2/SEAAD/analyses/run_variance_partition.py --modality ${MODALITY} --cell-type ${CELLTYPE} &
done

# run the age regression analysis per cell-type
REGRESSTYPES="wls vwrlm"
for MODALITY in ${MODALITIES[@]}; do
  for REGRESSTYPE in ${REGRESSTYPES[@]}; do
    uv run phase2/SEAAD/analyses/pseudobulk_regression.py --modality ${MODALITY} --cell-type ${CELLTYPE} \
      --regression-type ${REGRESSTYPE} --covariates specified --covariates-list PCA_0 PCA_1 PCA_2 PCA_3 \
      --only-autosomal "$DATADIR"/quants/aging_phase2.features.csv &
  done
done

# check the general regression against the robust to check for outlier driven results
CELLTYPES="ExN_CUX2 ExN_CUX2_RELN"
for CELLTYPE in ${CELLTYPES[@]}; do
  for MODALITY in ${MODALITIES[@]}; do
    uv run phase2/analyses/filter_regression_type_differences.py \
      --modality ${MODALITY} \
      --cell-type ${CELLTYPE} \
      --target-variable dx \
      --work-dir "$DATADIR"/public/seaad \
      --project seaad_ec_multiome
  done
done

# visualize similarities of disease status effects between course and ganular cell-type
for MODALITY in ${MODALITIES[@]}; do
  uv run phase2/figures/pairwise_celltype_effect_comparison.py \
    --modality ${MODALITY} \
    --cell-type-1 ExN_CUX2 \
    --cell-type-2 ${CELLTYPE} \
    --target-variable dx \
    --work-dir "$DATADIR"/public/seaad \
    --project seaad_ec_multiome
done
