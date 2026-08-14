#!/bin/bash
#SBATCH --job-name=surf-pressure
#SBATCH --account=def-plavoie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
# ================================
# run_compute_surface_pressure.sh
# ================================
# External, top-level entry point for the compute_surface_pressure_core
# package. This script lives OUTSIDE of the core package directory and is
# meant to be submitted directly with `sbatch run_compute_surface_pressure.sh`
# on a SLURM cluster, or run locally with `bash
# run_compute_surface_pressure.sh` for a quick interactive test.


set -euo pipefail

# ---- User-configurable paths -------------------------------------------------
SOL_DIR=${SOL_DIR:-"/project/rrg-plavoie/${USER}/FWH_solution"}
MESH_FILE=${MESH_FILE:-"/project/rrg-plavoie/${USER}/mesh/mesh.mesh.h5"}
WORKING_DIR=${WORKING_DIR:-"${SLURM_SUBMIT_DIR:-$(pwd)}/surface_pressure_output"}
SURFACE_PATCHES=${SURFACE_PATCHES:-"Airfoil_Surface"}

# Which post-processing steps to run (set to "true"/"false")
RUN_PSD=${RUN_PSD:-true}
RUN_CSD=${RUN_CSD:-false}
RUN_SURF_LINE=${RUN_SURF_LINE:-false}

# CSD-only settings (ignored if RUN_CSD=false)
MIC_DIR=${MIC_DIR:-"/project/rrg-plavoie/${USER}/microphones"}
MIC_FILE=${MIC_FILE:-"mic_data.hdf5"}
MIC_NUM=${MIC_NUM:-1}

FREQ_SELECT=${FREQ_SELECT:-"500 1000 1500 2000 2500 3000"}
BLOCK_SIZE=${BLOCK_SIZE:-1000}
NCHUNK=${NCHUNK:-4}

-------------------------------------------------------------------

# Directory containing this script -- also the directory containing
# run_compute_surface_pressure.py, one level above compute_surface_pressure_core/.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ARGS=(
    --sol-dir "${SOL_DIR}"
    --mesh-file "${MESH_FILE}"
    --surface-patches ${SURFACE_PATCHES}
    --working-dir "${WORKING_DIR}"
    --block-size "${BLOCK_SIZE}"
    --nchunk "${NCHUNK}"
    --freq-select ${FREQ_SELECT}
)

if [ "${RUN_PSD}" = true ]; then
    ARGS+=(--psd)
fi

if [ "${RUN_CSD}" = true ]; then
    ARGS+=(--csd --mic-dir "${MIC_DIR}" --mic-file "${MIC_FILE}" --mic-num "${MIC_NUM}")
fi

if [ "${RUN_SURF_LINE}" = true ]; then
    ARGS+=(--surf_line)
fi

echo "Running: python ${SCRIPT_DIR}/run_compute_surface_pressure.py ${ARGS[*]}"
python "${SCRIPT_DIR}/run_compute_surface_pressure.py" "${ARGS[@]}"