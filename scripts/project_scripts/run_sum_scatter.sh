#!/bin/bash

source ~/.bashrc

# This script expects:
#   ITERATION, OUTPUT_DIR, PYTHON, DATA_DIR
# It sums scatter outputs for the current iteration.

echo "Summing scatter outputs for iteration ${ITERATION}..."
echo "Input directory: ${OUTPUT_DIR}"

ITERATION_MINUS_ONE=$((ITERATION - 1))

$PYTHON sum_scatter.py \
    --input_dir=${OUTPUT_DIR} \
    --scatter_pattern="*_iter${ITERATION}_*_sca_w1.hs" \
    --total_pattern="*_iter${ITERATION}_*_tot_w1.hs" \
    --image_pattern="recon_osem_i*_s*_smoothed_${ITERATION_MINUS_ONE}.hv" \
    --output_file_prefix=${OUTPUT_DIR}/mean_iter${ITERATION} \
    --data_dir=${DATA_DIR} \
    --delete_files
