# pull neccessary files, fragments and raw matrix, for SEA-AD from Synapse
uv run python phase2/development/SEAAD/pull_synapse_filetypes.py

# detect ambient RNA
# get the file IDs (suffix) to process
ls public/seaad/src_data/*_raw_feature_bc_matrix.h5 >public/seaad/suffix_ids.list
sed -i s/_raw_feature_bc_matrix\.h5//g public/seaad/suffix_ids.list
sed -i s"/public\/seaad\/src_data\///"g public/seaad/suffix_ids.list
# format the swarm file to run the cellbender script per ID (file) found
rm /data/ADRD/brain_aging/phase2/scripts/seaad_cellbender.swarm
while read SUFFIXID; do
  echo /data/ADRD/brain_aging/phase2/scripts/run_cellbender.sh $SUFFIXID >>/data/ADRD/brain_aging/phase2/scripts/seaad_cellbender.swarm
done <public/seaad/suffix_ids.list
# submit the cellbender swarm
swarm -f /data/ADRD/brain_aging/phase2/scripts/seaad_cellbender.swarm \
  -g 32 \
  -t 2 \
  --module cellbender \
  --logdir /data/ADRD/brain_aging/phase2/hpc_logs \
  --time=12:00:00 \
  --partition=gpu \
  --gres=gpu:v100x:1 \
  --job-name swarm_cellbender_seaad
# check for failed jobs correct problem and resubmit
# looks like for 4 samples ran out of memory and killed (not VRAM actually mem)
while read SUFFIXID; do
  SAMPLEFILE="public/seaad/cellbender/"$SUFFIXID"_cellbender_filtered.h5"
  if [ ! -f "$SAMPLEFILE" ]; then
    echo "Missing file: $SAMPLEFILE"
    echo /data/ADRD/brain_aging/phase2/scripts/run_cellbender.sh $SUFFIXID >>/data/ADRD/brain_aging/phase2/scripts/seaad_retry_cellbender.swarm
  fi
done <public/seaad/suffix_ids.list
# submit the cellbender swarm
swarm -f /data/ADRD/brain_aging/phase2/scripts/seaad_retry_cellbender.swarm \
  -g 64 \
  -t 2 \
  --module cellbender \
  --logdir /data/ADRD/brain_aging/phase2/hpc_logs \
  --time=12:00:00 \
  --partition=gpu \
  --gres=gpu:v100x:1 \
  --job-name swarm_cellbender_retry_seaad

# detect cell doublets for RNA
DATADIR="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"
while read SUFFIXID; do
  uv run python /mnt/labshare/raph/working/ADRD_Brain_Aging/phase2/development/SEAAD/doublet_detection.py \
    --input-file "$DATADIR"/public/seaad/cellbender/"$SUFFIXID"_cellbender_filtered.h5 \
    --output-file "$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_rna_filtered.h5ad \
    --full-doublet-obs-file "$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_rna_unfiltered_obs.csv \
    --log-file "$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_rna.log \
    --plot-file "$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_rna.png
done <"$DATADIR"/public/seaad/suffix_ids.list

# quantify the chromatin peaks based on Phase2 Aging consensus peaks sed and detect doublets
while read SUFFIXID; do
  uv run python phase2/development/SEAAD/chromatin_peak_prep.py \
    --fragment-file "$DATADIR"/public/seaad/src_data/"$SUFFIXID"_atac_fragments.tsv.gz \
    --cpeaks-bed "$DATADIR"/src_data/aging_phase2_consensus_atac_peaks.bed \
    --output-file "$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_atac_filtered.h5ad \
    --min-tsse 2.5 \
    --tsse-plot "$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_atac.html \
    --log-file "$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_atac.log \
    --temp-dir "$DATADIR"/tmp &
done <"$DATADIR"/public/seaad/suffix_ids.list
wait

# make sure we have all the post-doublet files
while read SUFFIXID; do
  for MODALITY in "rna" "atac"; do
    SAMPLEFILE="$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_"$MODALITY"_filtered.h5ad
    if [ ! -f "$SAMPLEFILE" ]; then
      echo "Missing file: $SAMPLEFILE"
    fi
  done
done <"$DATADIR"/public/seaad/suffix_ids.list

# per sample combine the RNA and ATAC data into a mudata object, also integrate donor info into the obbs
uv run python phase2/development/SEAAD/combine_modalities.py --output-dir "$DATADIR"/public/seaad

# combine the per sample mudata objects into a single cohort mudata object
uv run python phase2/development/SEAAD/concat_mudata.py --output-file "$DATADIR"/public/seaad/seaad_ec_multiome.h5mu

# do preliminary cell-type label assignments using Celltypist
uv run python phase2/development/SEAAD/prep_celltypist_input.py \
  --input-file "$DATADIR"/public/seaad/seaad_ec_multiome.h5mu \
  --output-file "$DATADIR"/public/seaad/seaad_ec_rna_celltypist.h5ad \
  --detect-hv-features

# run Celltypist
sinteractive --constraint=gpuv100x --gres=lscratch:10,gpu:v100x:1 --mem=96g
module load rapids-singlecell
module load singularity
cd /data/ADRD/brain_aging/phase2
singularity pull celltypist-latest.sif docker://quay.io/teichlab/celltypist:latest
singularity run \
  -B /data/ADRD/brain_aging/phase2:/data \
  celltypist-latest.sif \
  celltypist --indata /data/public/seaad/seaad_ec_rna_celltypist.h5ad \
  --model /data/celltypist/Adult_Human_PrefrontalCortex.pkl --majority-voting \
  --prefix seaad_ec_rna_celltypist_PFC_ --outdir /data/public/seaad

# check initial clustering based on only RNA
uv run python phase2/development/SEAAD/initial_clustering_rna.py \
  --input-file "$DATADIR"/public/seaad/seaad_ec_multiome.h5mu \
  --output-file "$DATADIR"/public/seaad/seaad_ec_rna_clustered.h5ad \
  --model-dir "$DATADIR"/public/seaad/seaad_ec_rna_scvi \
  --detect-hv-features

# integrate the cell-type predicted labels to the initial clustering
uv run python phase2/development/SEAAD/annotate_clusters.py \
  --input-h5ad "$DATADIR"/public/seaad/seaad_ec_rna_clustered.h5ad \
  --predictions-csv "$DATADIR"/public/seaad/seaad_ec_rna_celltypist_PFC_predicted_labels.csv \
  --output-h5ad "$DATADIR"/public/seaad/seaad_ec_rna_annotated.h5ad \
  --output-plot "$DATADIR"/public/seaad/umap_celltype_annotations.png

# perform subclustering by broad cell-type
# excitatory neurons
uv run python phase2/development/SEAAD/subcluster_rna.py \
  --raw-input "$DATADIR"/public/seaad/seaad_ec_multiome.h5mu \
  --processed-h5ad "$DATADIR"/public/seaad/seaad_ec_rna_annotated.h5ad \
  --model-dir "$DATADIR"/public/seaad/models/seaad_ec_ExN_scvi \
  --resolution-dir "$DATADIR"/public/seaad/resolution_selection \
  --name-prefix "seaad_ec_ExN" \
  --subset-col leiden_scvi \
  --subset-values "0,10,11,2,3,8,9" \
  --annotation-col celltype_major \
  --detect-hv-features
# inhibitory neurons
uv run python phase2/development/SEAAD/subcluster_rna.py \
  --raw-input "$DATADIR"/public/seaad/seaad_ec_multiome.h5mu \
  --processed-h5ad "$DATADIR"/public/seaad/seaad_ec_rna_annotated.h5ad \
  --model-dir "$DATADIR"/public/seaad/models/seaad_ec_InN_scvi \
  --resolution-dir "$DATADIR"/public/seaad/resolution_selection \
  --name-prefix "seaad_ec_InN" \
  --subset-col leiden_scvi \
  --subset-values "1,12,13,14,15,4,5,6" \
  --annotation-col celltype_major \
  --detect-hv-features
# non-neuronal
uv run python phase2/development/SEAAD/subcluster_rna.py \
  --raw-input "$DATADIR"/public/seaad/seaad_ec_multiome.h5mu \
  --processed-h5ad "$DATADIR"/public/seaad/seaad_ec_rna_annotated.h5ad \
  --model-dir "$DATADIR"/public/seaad/models/seaad_ec_NonNeuronal_scvi \
  --resolution-dir "$DATADIR"/public/seaad/resolution_selection \
  --name-prefix "seaad_ec_NonNeuronal" \
  --subset-col leiden_scvi \
  --subset-values "16,17,18,19,20,7" \
  --annotation-col celltype_major \
  --detect-hv-features
