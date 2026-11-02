#!/bin/bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3GB
#SBATCH --time=08:00:00
#SBATCH --array=0-2
#SBATCH --mail-type=BEGIN,END,FAIL   # Alerts for start, finish, and crashes
#SBATCH --mail-user=bhr597@sund.ku.dk

set -eo pipefail

# Load Conda into the shell
module restore myhmscstack 
source /opt/software/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate my_hmsc_env

MODELS=(
  "2026-07-27_13-28-41_Atlas1_MinOccs5_CoverageGoodAverage"
  "2026-07-27_13-28-41_Atlas2_MinOccs5_CoverageGoodAverage"
  "2026-07-27_13-28-41_Atlas3_MinOccs5_CoverageGoodAverage"
)
THIN_VAL=10

MODEL="${MODELS[$SLURM_ARRAY_TASK_ID]}"
MODEL_FILE="HmscOutputs/${MODEL}/Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata"

if [[ ! -f "$MODEL_FILE" ]]; then
  echo "Required fitted model not found: $MODEL_FILE" >&2
  exit 1
fi

srun Rscript --verbose S6_Making_Predictions.R "$THIN_VAL" "$MODEL"