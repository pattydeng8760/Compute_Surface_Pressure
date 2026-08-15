#!/bin/bash
#SBATCH --job-name=surf-pressure
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --account=rrg-plavoie

# =============================================================================
# run_compute_surface_pressure.sh
# =============================================================================
#
# SLURM launcher for run_compute_surface_pressure.py.
#
# This script lives OUTSIDE compute_surface_pressure_core/ and passes the
# case-specific settings below to the Python runner as command-line arguments.
#
# Submit with:
#
#     sbatch run_compute_surface_pressure.sh
#
# The Python runner gives explicit CLI arguments precedence over the CONFIG
# dictionary in run_compute_surface_pressure.py. Therefore, when this shell
# script is used, the values defined here control the calculation.
# =============================================================================
source /project/rrg-plavoie/denggua1/pd_env.sh
# Required inputs
SOL_DIR=${SOL_DIR:-"/scratch/denggua1/Bombardier_LES/B_10AOA_LES/RUN_Fine_Jul25/FWH_Airfoil/FWH_Data_TTG"}
MESH_FILE=${MESH_FILE:-"/project/rrg-plavoie/denggua1/BBDB_10AOA/MESH_Fine_Jul25/Bombardier_10AOA_U30_Combine_Fine.mesh.h5"}

# General settings
WORKING_DIR=${WORKING_DIR:-"./surface_pressure_output"}

SURFACE_PATCHES=(
    "Airfoil_Surface"
    "Airfoil_Trailing_Edge"
    "Airfoil_Side_LE"
    "Airfoil_Side_Mid"
    "Airfoil_Side_TE"
)

FREQ_SELECT=(500 1000 1500 2000 2500 3000)

OPTION=${OPTION:-1}
NSKIP=${NSKIP:-1}
MAX_FILE=${MAX_FILE:-5000}
VAR=${VAR:-"pressure"}
RELOAD=${RELOAD:-false}

# Workflows
RUN_PSD=${RUN_PSD:-true}
RUN_CSD=${RUN_CSD:-true}
RUN_FFT=${RUN_FFT:-false}
RUN_SURF_LINE=${RUN_SURF_LINE:-false}

# PSD/CSD settings
BLOCK_SIZE=${BLOCK_SIZE:-1000}
NCHUNK=${NCHUNK:-4}
BAND_STAT=${BAND_STAT:-"mean"}
DF=${DF:-50.0}

# CSD settings
MIC_DIR=${MIC_DIR:-"/scratch/denggua1/Bombardier_LES/B_10AOA_LES/RUN_Fine_Jul25/FWH_Airfoil/FWHpostpro_TTG"}
MIC_FILE=${MIC_FILE:-"B_10AOA_U30_LES_Solid_Airfoil_Mic_Data.h5"}
MIC_NUM=${MIC_NUM:-68}

# Surface-line settings
AOA=${AOA:-10.0}
UINF=${UINF:-30.0}

ORIENTATION=( 0 0 1)

CUT_LOC_PERCENT=${CUT_LOC_PERCENT:-70.0}
Z_LOC_PERCENT=${Z_LOC_PERCENT:-70.0}
CUT_LENGTH_PERCENT=${CUT_LENGTH_PERCENT:-10.0}
DATA_SIZE=${DATA_SIZE:-1000}

AIRFOIL_FILE=${AIRFOIL_FILE:-""}
CAMBER_FILE=${CAMBER_FILE:-""}


# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}


# =============================================================================
# BUILD CLI
# =============================================================================

ARGS=(
    --sol-dir "${SOL_DIR}"
    --mesh-file "${MESH_FILE}"
    --surface-patches "${SURFACE_PATCHES[@]}"
    --working-dir "${WORKING_DIR}"
    --freq-select "${FREQ_SELECT[@]}"
    --option "${OPTION}"
    --nskip "${NSKIP}"
    --max-file "${MAX_FILE}"
    --var "${VAR}"
    --block-size "${BLOCK_SIZE}"
    --nchunk "${NCHUNK}"
    --band-stat "${BAND_STAT}"
    --df "${DF}"
)

if [[ "${RELOAD}" == "true" ]]; then
    ARGS+=(--reload)
fi

if [[ "${RUN_PSD}" == "true" ]]; then
    ARGS+=(--psd)
fi

if [[ "${RUN_CSD}" == "true" ]]; then
    ARGS+=(
        --csd
        --mic-dir "${MIC_DIR}"
        --mic-file "${MIC_FILE}"
        --mic-num "${MIC_NUM}"
    )
fi

if [[ "${RUN_FFT}" == "true" ]]; then
    ARGS+=(--fft)
fi

if [[ "${RUN_SURF_LINE}" == "true" ]]; then
    ARGS+=(
        --surf_line
        --AoA "${AOA}"
        --Uinf "${UINF}"
        --orientation "${ORIENTATION[@]}"
        --cut-loc-percent "${CUT_LOC_PERCENT}"
        --z-loc-percent "${Z_LOC_PERCENT}"
        --cut-length-percent "${CUT_LENGTH_PERCENT}"
        --data-size "${DATA_SIZE}"
    )

    if [[ -n "${AIRFOIL_FILE}" ]]; then
        ARGS+=(--airfoil-file "${AIRFOIL_FILE}")
    fi

    if [[ -n "${CAMBER_FILE}" ]]; then
        ARGS+=(--camber-file "${CAMBER_FILE}")
    fi
fi


# =============================================================================
# RUN
# =============================================================================

echo
echo "================================================================================"
echo "Surface Pressure Post-Processing"
echo "================================================================================"
echo "Script directory  : ${SCRIPT_DIR}"
echo "Solution directory: ${SOL_DIR}"
echo "Mesh file         : ${MESH_FILE}"
echo "Working directory : ${WORKING_DIR}"
echo "Surface patches   : ${SURFACE_PATCHES[*]}"
echo "Frequencies [Hz]  : ${FREQ_SELECT[*]}"
echo "PSD               : ${RUN_PSD}"
echo "CSD               : ${RUN_CSD}"
echo "FFT               : ${RUN_FFT}"
echo "Surface line      : ${RUN_SURF_LINE}"
echo "Reload            : ${RELOAD}"
echo "SLURM job ID      : ${SLURM_JOB_ID:-N/A}"
echo "Allocated CPUs    : ${SLURM_CPUS_PER_TASK:-N/A}"
echo "================================================================================"
echo

printf 'Running command:\npython %q' "${SCRIPT_DIR}/run_compute_surface_pressure.py"
printf ' %q' "${ARGS[@]}"
printf '\n\n'

python "${SCRIPT_DIR}/run_compute_surface_pressure.py" "${ARGS[@]}"
