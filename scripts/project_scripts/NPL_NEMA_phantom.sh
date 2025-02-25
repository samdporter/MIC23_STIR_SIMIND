#!/bin/bash
set -euo pipefail

# Define common variables
PYTHON="python"
DATA_PATH="/home/storage/copied_data/data/phantom_data/for_cluster/SPECT"
SCRIPTS_DIR="/scripts/project_scripts"
BASE_DIR="/home/sam/working/STIR_users_MIC2023"
OUTPUT_DIR="simind_output"
INITIAL_SUBSETS=12
INITIAL_EPOCHS=10

# Initial reconstruction with index 0
$PYTHON scripts/project_scripts/osem.py \
    --num_subsets="${INITIAL_SUBSETS}" \
    --num_epochs="${INITIAL_EPOCHS}" \
    --data_path="${DATA_PATH}" \
    --smoothing=True \
    --index=0

for i in {0..4}; do
    IMAGE_PATH="${DATA_PATH}/recon_osem_i${INITIAL_EPOCHS}_s${INITIAL_SUBSETS}_smoothed_${i}.hv"
    OUTPUT_PREFIX="output_${i}"
    ADDITIVE_PATH="${BASE_DIR}/${OUTPUT_DIR}/${OUTPUT_PREFIX}_sca_w1.hs"


    # Note: currently set not to model septal penetration
    # So collimator_routine=0 (set to 1 for septal penetration)
    # And photon direction=2 (set to 3 for septal penetration)
    $PYTHON "${BASE_DIR}/STIR_demonstration.py" \
        --total_activity=187 \
        --time_per_projection=40 \
        --photon_multiplier=1000 \
        --photopeak_energy=150 \
        --window_lower=75 \
        --window_upper=225 \
        --source_type="y90_frey" \
        --collimator="ma-megp" \
        --kev_per_channel=25 \
        --max_energy=2180 \
        --mu_map_path="${BASE_DIR}/data/Y90/umap_zoomed.hv" \
        --image_path="${IMAGE_PATH}" \
        --measured_data_path="${BASE_DIR}/data/Y90/peak_1_projdata__f1g1d0b0.hs" \
        --measured_additive_path=None \
        --output_prefix="${OUTPUT_PREFIX}" \
        --output_dir="${OUTPUT_DIR}" \
        --input_smc_file_path="${BASE_DIR}/input/input.smc" \
        --scoring_routine=1 \
        --collimator_routine=1 \
        --photon_direction=3 \
        --crystal_thickness=15.9 \
        --crystal_half_length_radius=235 \
        --crystal_half_width=285 \
        --half_life=64.6 \
        --axial_slice=56

    $PYTHON "${BASE_DIR}/${SCRIPTS_DIR}/osem.py" \
        --num_subsets="${INITIAL_SUBSETS}" \
        --num_epochs="${INITIAL_EPOCHS}" \
        --additive_path="${ADDITIVE_PATH}" \
        --data_path="${DATA_PATH}" \
        --smoothing=True \
        --index=$((i+1))
done
