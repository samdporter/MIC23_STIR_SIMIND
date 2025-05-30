#!/bin/bash
set -euo pipefail

# === Configuration ===
PYTHON=python3
DATA_DIR="/home/storage/prepared_data/phantom_data/anthropomorphic_phantom_data/SPECT/phantom_140"
BASE_DIR="/home/sam/working/STIR_users_MIC2023"
SCRIPTS_DIR="${BASE_DIR}/scripts/project_scripts"
SUFFIX="phantom_140"
OUTPUT_DIR="${BASE_DIR}/output/${SUFFIX}"
INITIAL_SUBSETS=12
INITIAL_EPOCHS=1
TOTAL_ACTIVITY=160.8 #145.5 cyl #187 nema #160.8 anthro # in MBq
PHOTON_MULTIPLIER=100
NUM_ITERATIONS=5
NUM_ARRAY_JOBS=1
WINDOW_LOWER=75
WINDOW_UPPER=225
PHOTON_ENERGY=150

# === Prep output dir ===
mkdir -p "${OUTPUT_DIR}"
rm -rf "${OUTPUT_DIR:?}/"*

# === Initial OSEM ===
echo "Running initial OSEM reconstruction..."
env \
  DATA_DIR="${DATA_DIR}" \
  OUTPUT_DIR="${OUTPUT_DIR}" \
  SCRIPTS_DIR="${SCRIPTS_DIR}" \
  INITIAL_SUBSETS="${INITIAL_SUBSETS}" \
  INITIAL_EPOCHS="${INITIAL_EPOCHS}" \
  BASE_DIR="${BASE_DIR}" \
  PYTHON="${PYTHON}" \
  "${SCRIPTS_DIR}/run_osem.sh"

# === Iterative loop ===
for (( i=1; i<=NUM_ITERATIONS; i++ )); do
  echo "=== Iteration ${i} ==="

  # Simulation “array” (1..NUM_ARRAY_JOBS)
  for (( task=1; task<=NUM_ARRAY_JOBS; task++ )); do
    echo "  Simulation task ${task}/${NUM_ARRAY_JOBS}..."
    env \
      ITERATION="${i}" \
      TASK_ID="${task}" \
      DATA_DIR="${DATA_DIR}" \
      OUTPUT_DIR="${OUTPUT_DIR}" \
      BASE_DIR="${BASE_DIR}" \
      TOTAL_ACTIVITY="${TOTAL_ACTIVITY}" \
      PYTHON="${PYTHON}" \
      INITIAL_SUBSETS="${INITIAL_SUBSETS}" \
      INITIAL_EPOCHS="${INITIAL_EPOCHS}" \
      PHOTON_MULTIPLIER="${PHOTON_MULTIPLIER}" \
      WINDOW_LOWER="${WINDOW_LOWER}" \
      WINDOW_UPPER="${WINDOW_UPPER}" \
      PHOTON_ENERGY="${PHOTON_ENERGY}" \
      "${SCRIPTS_DIR}/run_simulation_array.sh"
  done

  # Sum scatter
  echo "  Summing scatter for iteration ${i}..."
  env \
    ITERATION="${i}" \
    OUTPUT_DIR="${OUTPUT_DIR}" \
    PYTHON="${PYTHON}" \
    DATA_DIR="${DATA_DIR}" \
    "${SCRIPTS_DIR}/run_sum_scatter.sh"

  # OSEM update
  echo "  Running OSEM update for iteration ${i}..."
  env \
    ITERATION="${i}" \
    DATA_DIR="${DATA_DIR}" \
    OUTPUT_DIR="${OUTPUT_DIR}" \
    SCRIPTS_DIR="${SCRIPTS_DIR}" \
    INITIAL_SUBSETS="${INITIAL_SUBSETS}" \
    INITIAL_EPOCHS="${INITIAL_EPOCHS}" \
    BASE_DIR="${BASE_DIR}" \
    PYTHON="${PYTHON}" \
    "${SCRIPTS_DIR}/run_osem.sh"
done

echo "All iterations complete."
