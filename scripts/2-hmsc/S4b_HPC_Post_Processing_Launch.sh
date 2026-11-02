#!/bin/bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=20
#SBATCH --mem=48GB
#SBATCH --time=24:00:00
#SBATCH --array=0-2
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=bhr597@sund.ku.dk

# Load Conda into the shell
module restore myhmscstack 
source /opt/software/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate my_hmsc_env

MODELS=(
  "2026-07-27_13-28-41_Atlas1_MinOccs5_CoverageGoodAverage"
  "2026-07-27_13-28-41_Atlas2_MinOccs5_CoverageGoodAverage"
  "2026-07-27_13-28-41_Atlas3_MinOccs5_CoverageGoodAverage"
)
THIN_VAL=1

srun Rscript --verbose S4b_HPC_Post_Processing.R \
  "$THIN_VAL" "${MODELS[$SLURM_ARRAY_TASK_ID]}"
