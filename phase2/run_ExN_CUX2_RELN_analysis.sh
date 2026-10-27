# Analyze ExN CUX2 RELN in the SEA-AD EC data
# since the SAHA analysis we have for the SEAAD EC data has the more granular cell-types already included can use existing utilites
# since the many of the utilities will run for all cell-types at granularity speciefied, which we don't need, just the one subtype
# subset the mudata object to the ExN CUX2 and then run utitlities specifying anno instead of anno_coarse
# convert to pseudobulk per cell-type
uv run python -c "
import mudata as md
data_dir = '/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2'
mdata = md.read_h5mu(f'{data_dir}/public/seaad/seaad_ec_multiome_labeled.h5mu')
cux2_reln = mdata[mdata.obs.anno == 'ExN CUX2 RELN']
print(cux2_reln)
cux2_reln.obs.anno = cux2_reln.obs.anno.cat.remove_unused_categories()
cux2_reln.write(f'{data_dir}/public/seaad/seaad_ec_multiome_exn_cux2_reln.h5mu')
"

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
PROJECTNAME="seaad_ec_multiome"
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
  phase2/SEAAD/run_var_part_jobs.sh run_variance_partition.py ${MODALITY}
done

# since using WLS check cell-types and modalities for correlations between cell counts and age
uv run phase2/analyses/cell_counts_regression.py --target-variable dx \
  --work-dir "$DATADIR"/public/seaad --project seaad_ec_multiome

# run the age regression analysis per cell-type
REGRESSTYPES="wls vwrlm"
for MODALITY in ${MODALITIES[@]}; do
  for REGRESSTYPE in ${REGRESSTYPES[@]}; do
    phase2/SEAAD/run_regression_jobs.sh ${REGRESSTYPE} ${MODALITY}
  done
done

# post process the age regressions
for MODALITY in ${MODALITIES[@]}; do
  for REGRESSTYPE in ${REGRESSTYPES[@]}; do
    uv run phase2/analyses/post_pseudobulk_regression.py \
      --modality ${MODALITY} \
      --regression-type ${REGRESSTYPE} \
      --target-variable dx \
      --no-volcano-per-celltype \
      --work-dir "$DATADIR"/public/seaad \
      --project seaad_ec_multiome
  done
done

# check the general regression against the robust to check for outlier driven results
for MODALITY in ${MODALITIES[@]}; do
  uv run phase2/analyses/filter_regression_type_differences.py \
    --modality ${MODALITY} \
    --target-variable dx \
    --work-dir "$DATADIR"/public/seaad \
    --project seaad_ec_multiome
done

# visualize similarities of disease status effects between tissues
## union of significant features as heatmap
for MODALITY in ${MODALITIES[@]}; do
  uv run phase2/figures/celltype_target_effect_similarity.py \
    --modality ${MODALITY} \
    --target-variable dx \
    --work-dir "$DATADIR"/public/seaad \
    --project seaad_ec_multiome
done
