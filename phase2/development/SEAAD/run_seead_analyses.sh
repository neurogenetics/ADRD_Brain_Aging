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
    --log-file "$DATADIR"/public/seaad/doublet_det/"$SUFFIXID"_atac.log \
    --temp-dir "$DATADIR"/tmp &
done <"$DATADIR"/public/seaad/suffix_ids.list
