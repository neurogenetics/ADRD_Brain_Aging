#!/bin/bash

echo "Executing: ${0} $@"

for POOL in {1..6}
do
    for LANE in {1..8}
    do
        SAMPLE="Aging_P00${POOL}_SCRN_${LANE}"
        echo ${SAMPLE}
        if test -f testing/${SAMPLE}_cellbender_filtered.h5
        then
            echo "${SAMPLE} already done skipping"        
        else
            gsutil cp gs://nihnialng-aging-brain/nisc/${SAMPLE}/outs/raw_feature_bc_matrix.h5 testing/${SAMPLE}_raw.h5
            gsutil cp gs://nihnialng-aging-brain/nisc/${SAMPLE}/outs/filtered_feature_bc_matrix.h5 testing/${SAMPLE}_filtered.h5
            docker run -v $(pwd)/testing:/data 56439f37d58e cellbender remove-background --input /data/${SAMPLE}_raw.h5 --output /data/${SAMPLE}_cellbender.h5
        fi
    done
done