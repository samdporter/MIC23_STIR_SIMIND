#!/bin/bash

source ~/.bashrc

# This script expects the following environment variables:
#   ITERATION, DATA_DIR, OUTPUT_DIR, BASE_DIR, TOTAL_ACTIVITY, PYTHON, INITIAL_SUBSETS, INITIAL_EPOCHS
# SGE_TASK_ID is automatically provided.

# Construct a unique output prefix for each array task.
OUTPUT_PREFIX="output_iter${ITERATION}_${SGE_TASK_ID}"

# Determine the input image for simulation. For iteration i, use the reconstruction from iteration i-1.
if [ "${ITERATION}" -eq 1 ]; then
    INPUT_IMAGE="${OUTPUT_DIR}/recon_osem_i${INITIAL_EPOCHS}_s${INITIAL_SUBSETS}_smoothed_0.hv"
else
    PREV_ITER=$(( ITERATION - 1 ))
    INPUT_IMAGE="${OUTPUT_DIR}/recon_osem_i${INITIAL_EPOCHS}_s${INITIAL_SUBSETS}_smoothed_${PREV_ITER}.hv"
fi

$PYTHON ${BASE_DIR}/simulation_script.py \
    --total_activity=${TOTAL_ACTIVITY} \
    --time_per_projection=40 \
    --photon_multiplier=${PHOTON_MULTIPLIER} \
    --photopeak_energy=${PHOTON_ENERGY} \
    --window_lower=${WINDOW_LOWER} \
    --window_upper=${WINDOW_UPPER} \
    --source_type="y90_frey" \
    --collimator="ma-megp" \
    --kev_per_channel=20 \
    --max_energy=960 \
    --mu_map_path=${DATA_DIR}/umap_zoomed.hv \
    --image_path=${INPUT_IMAGE} \
    --measured_data_path=${DATA_DIR}/peak.hs \
    --output_prefix=${OUTPUT_PREFIX} \
    --output_dir=${OUTPUT_DIR} \
    --input_smc_file_path=${BASE_DIR}/input/input.smc \
    --simind_parent_dir=${BASE_DIR} \
    --scoring_routine=1 \
    --collimator_routine=1 \
    --photon_direction=3 \
    --crystal_thickness=15.9 \
    --crystal_half_length_radius=195 \
    --crystal_half_width=265 \
    --half_life=64.6 \
    --axial_slice=56

# remember to set collimator_routine=1 and photon_direction=3 for septal penetration
# time per projection = 20 for patients
