#!/bin/bash

SAMPLEID=${1}

WRKDIR="/data/ADRD/brain_aging/phase2"
THISJOBDIR="$WRKDIR/hpc_logs/$SAMPLEID"
mkdir -p "$THISJOBDIR"
INPUT="$WRKDIR/public/seaad/src_data/"$SAMPLEID"_raw_feature_bc_matrix.h5"
OUTPUT="$WRKDIR/public/seaad/cellbender/"$SAMPLEID"_cellbender.h5"

module load cellbender

cd "$THISJOBDIR"
cellbender remove-background \
  --cuda \
  --input $INPUT \
  --output $OUTPUT
