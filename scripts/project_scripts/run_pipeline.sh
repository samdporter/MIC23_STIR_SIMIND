#!/bin/bash

# Define common variables (adjust paths and parameters as needed)
PYTHON=python3
DATA_DIR="/home/sporter/synergistic_Y90/prepared_data/phantom_data/anthropomorphic_phantom_data/SPECT/cylindrical_208"
BASE_DIR="/home/sporter/synergistic_Y90/MIC23_STIR_SIMIND"
SCRIPTS_DIR="${BASE_DIR}/scripts/project_scripts"
SUFFIX="cylindrical_phantom_208"
OUTPUT_DIR="${BASE_DIR}/output/${SUFFIX}"
INITIAL_SUBSETS=12
INITIAL_EPOCHS=10
TOTAL_ACTIVITY=264.7 # 187 # 182.8
PHOTON_MULTIPLIER=1
NUM_ITERATIONS=1
NUM_ARRAY_JOBS=100  # Allow overriding with an environment variable or command-line argument

# Ensure the output directory exists
mkdir -p "${OUTPUT_DIR}"
# ensure the output directory is empty
rm -rf ${OUTPUT_DIR}/*

# Ensure the script is run from SCRIPTS_DIR
if [[ "$PWD" != "$SCRIPTS_DIR" ]]; then
    echo "Error: Please run this script from ${SCRIPTS_DIR}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# Function to extract the numeric job ID from qsub output
extract_job_id() {
    # Expects output like: "Your job 5334534.1-100:1 has been submitted."
    # Extract the third field and remove everything after the first dot.
    echo "$1" | awk '{print $3}' | cut -d. -f1
}

echo "Submitting initial OSEM reconstruction..."
INIT_OUT=$(qsub \
    -N init_osem_${SUFFIX} \
    -cwd -l h_rt=04:00:00,tmem=32G,h_vmem=32G,tscratch=10G \
    -j y -R y \
    -v DATA_DIR="${DATA_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",SCRIPTS_DIR="${SCRIPTS_DIR}",INITIAL_SUBSETS="${INITIAL_SUBSETS}",INITIAL_EPOCHS="${INITIAL_EPOCHS}",BASE_DIR="${BASE_DIR}",PYTHON="${PYTHON}" \
    "${SCRIPTS_DIR}/run_osem.sh")
INIT_JOB=$(extract_job_id "${INIT_OUT}")
echo "Initial OSEM job submitted with ID ${INIT_JOB}."

# Initialize dependency chain with the initial job
PREV_JOB=${INIT_JOB}

# Loop over iterations using scheduler dependencies
for i in $(seq 1 ${NUM_ITERATIONS}); do
    echo "Submitting jobs for iteration ${i}"
    
    SIM_OUT=$(qsub \
        -N sim_iter_${i}_${SUFFIX} \
        -cwd -l h_rt=96:00:00,tmem=16G,h_vmem=16G,tscratch=10G \
        -t 1-${NUM_ARRAY_JOBS} \
        -j y -R y \
        -hold_jid "${PREV_JOB}" \
        -v ITERATION="${i}",DATA_DIR="${DATA_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BASE_DIR="${BASE_DIR}",TOTAL_ACTIVITY="${TOTAL_ACTIVITY}",PYTHON="${PYTHON}",INITIAL_SUBSETS="${INITIAL_SUBSETS}",INITIAL_EPOCHS="${INITIAL_EPOCHS}" \
        "${SCRIPTS_DIR}/run_simulation_array.sh")
    SIM_JOB=$(extract_job_id "${SIM_OUT}")
    echo "Simulation job for iteration ${i} submitted with ID ${SIM_JOB}."
    
    SUM_OUT=$(qsub \
        -N sum_iter_${i}_${SUFFIX} \
        -cwd -l h_rt=01:00:00,tmem=32G,h_vmem=32G,tscratch=10G \
        -j y -R y \
        -hold_jid "${SIM_JOB}" \
        -v ITERATION="${i}",OUTPUT_DIR="${OUTPUT_DIR}",PYTHON="${PYTHON}",DATA_DIR="${DATA_DIR}" \
        "${SCRIPTS_DIR}/run_sum_scatter.sh")
    SUM_JOB=$(extract_job_id "${SUM_OUT}")
    echo "Scatter summing job for iteration ${i} submitted with ID ${SUM_JOB}."
    
    OSEM_OUT=$(qsub \
        -N osem_iter_${i}_${SUFFIX} \
        -cwd -l h_rt=04:00:00,tmem=32G,h_vmem=32G,tscratch=10G \
        -j y -R y \
        -hold_jid "${SUM_JOB}" \
        -v ITERATION="${i}",DATA_DIR="${DATA_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",SCRIPTS_DIR="${SCRIPTS_DIR}",INITIAL_SUBSETS="${INITIAL_SUBSETS}",INITIAL_EPOCHS="${INITIAL_EPOCHS}",BASE_DIR="${BASE_DIR}",PYTHON="${PYTHON}" \
        "${SCRIPTS_DIR}/run_osem.sh")
    OSEM_JOB=$(extract_job_id "${OSEM_OUT}")
    echo "OSEM update job for iteration ${i} submitted with ID ${OSEM_JOB}."
    
    # Set dependency for next iteration: simulation job waits on the current OSEM update
    PREV_JOB=${OSEM_JOB}
done
