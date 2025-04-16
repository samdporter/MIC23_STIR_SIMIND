#!/bin/bash

source ~/.bashrc

# This script expects the following environment variables:
#   ITERATION (optional; if not set, we treat it as the initial reconstruction)
#   DATA_DIR, OUTPUT_DIR, SCRIPTS_DIR, INITIAL_SUBSETS, INITIAL_EPOCHS, BASE_DIR, PYTHON

if [ -z "${ITERATION:-}" ]; then
    INDEX="0"
    ADDITIVE=""
else
    INDEX="${ITERATION}"
    # For iteration i, use the mean scatter file computed in the summing job.
    ADDITIVE="--additive_path=${OUTPUT_DIR}/mean_iter${ITERATION}_scatter.hs"
fi

$PYTHON ${SCRIPTS_DIR}/osem.py \
    --num_subsets=${INITIAL_SUBSETS} \
    --num_epochs=${INITIAL_EPOCHS} \
    --data_path=${DATA_DIR} \
    --output_path=${OUTPUT_DIR} \
    $ADDITIVE \
    --smoothing=True \
    --index=${INDEX}
