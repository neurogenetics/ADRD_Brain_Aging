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
