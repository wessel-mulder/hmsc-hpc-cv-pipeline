#!/bin/bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=04:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bhr597@sund.ku.dk
#SBATCH --array=0-2
# Run from repository root: sbatch scripts/2-hmsc/S4b_HPC_Post_Processing_Launch.sh
set -euo pipefail
source scripts/2-hmsc/helpers/config.sh
: "${SLURM_ARRAY_TASK_ID:?Submit as an array job}"
module restore myhmscstack
source /opt/software/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate "$CH1_POSTPROCESS_CONDA_ENV"
models_file=$(mktemp)
trap 'rm -f "$models_file"' EXIT
Rscript scripts/2-hmsc/helpers/list-models.R > "$models_file"
models=()
while IFS= read -r model; do models+=("$model"); done < "$models_file"
(( SLURM_ARRAY_TASK_ID >= 0 && SLURM_ARRAY_TASK_ID < ${#models[@]} )) || exit 1
srun Rscript scripts/2-hmsc/S4b_HPC_Post_Processing.R "${models[$SLURM_ARRAY_TASK_ID]}"
