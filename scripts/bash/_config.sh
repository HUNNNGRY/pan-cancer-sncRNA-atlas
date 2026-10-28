#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# _config.sh - shared configuration for the command-line steps of the atlas
#
# Every script in this directory sources this file first, so all paths and
# cluster settings are defined in one place.
#
#   export SNCRNA_HOME=/path/to/data/mirror
#   bash scripts/bash/03_pcc_mrna.sh
#
# Jobs are submitted with `sbatch`. On a machine without SLURM set
# SCHEDULER=none to run them in the foreground, or leave PARTITION empty.
# ---------------------------------------------------------------------------
set -euo pipefail

SNCRNA_HOME="${SNCRNA_HOME:-${HOME}/sncRNA-atlas-data}"

# repository layout
BASH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${BASH_DIR}/../.." && pwd)"

# data roots (keep in sync with scripts/config.R)
DATA_DIR="${DATA_DIR:-${SNCRNA_HOME}/projects/WCHSU-FTC}"
PROJECT_DIR="${PROJECT_DIR:-${SNCRNA_HOME}/projects/motif-RBP-EDA}"
TCGA_DIR="${TCGA_DIR:-${SNCRNA_HOME}/projects/tcga}"
CIBERSORT_DIR="${CIBERSORT_DIR:-${SNCRNA_HOME}/projects/lulab_tissue_pico}"
SHARED_UTILS_DIR="${SHARED_UTILS_DIR:-${SNCRNA_HOME}/shared_utils}"
SHARED_SCRIPTS_DIR="${SHARED_SCRIPTS_DIR:-${SNCRNA_HOME}/shared_scripts}"
EXSEEK_SCRIPTS_DIR="${EXSEEK_SCRIPTS_DIR:-${DATA_DIR}/exSeek-dev/scripts}"

# dataset and the 16 TCGA cancer types used in the manuscript
DST="${DST:-TCGA_small16}"
CANCERS=("BRCA" "COAD" "HNSC" "LIHC" "LUAD" "PRAD" "STAD" "THCA"
         "LUSC" "ESCA" "BLCA" "CHOL" "UCEC" "KIRC" "KIRP" "KICH")

# resource settings
CORES="${CORES:-12}"
PARTITION="${PARTITION:-}"       # e.g. "batch"; empty disables sbatch
SBATCH_EXTRA="${SBATCH_EXTRA:-}" # extra sbatch flags, e.g. "--mem=64G"
ENV_PCC="${ENV_PCC:-cfpeak}"     # conda env used by run-covarPCC.R
ENV_R="${ENV_R:-r41}"            # conda env used by the pathway helpers

mkdir -p "${PROJECT_DIR}/tmp" "${PROJECT_DIR}/logs"

# ---------------------------------------------------------------- activation
if command -v conda >/dev/null 2>&1; then
  CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
else
  CONDA_SH=""
fi

# activate_cmd <env> : shell snippet that activates a conda environment
activate_cmd() {
  if [ -n "${CONDA_SH}" ] && [ -f "${CONDA_SH}" ]; then
    echo "source '${CONDA_SH}' && conda activate '$1'"
  else
    echo "source activate '$1'"
  fi
}

# submit <job-name> <log-file> <command>
# Runs <command> on the cluster, or in the foreground when SLURM is disabled.
submit() {
  local name="$1" log="$2" cmd="$3"
  mkdir -p "$(dirname "${log}")"
  if [ "${SCHEDULER:-slurm}" = "none" ] || [ -z "${PARTITION}" ]; then
    echo "[run] ${name}"
    bash -c "${cmd}" 2>&1 | tee "${log}"
  else
    # shellcheck disable=SC2086
    sbatch -p "${PARTITION}" -N 1 -n "${CORES}" ${SBATCH_EXTRA} \
      --job-name="${name}" --output="${log}" --error="${log}" --wrap "${cmd}"
  fi
}
