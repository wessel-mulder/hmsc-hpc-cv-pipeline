#!/bin/bash
#SBATCH --partition=gpuqueue
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1GB
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --array=0-3
# Submit via S2a_Fit_HMSC_Model_Launch.sh from the repository root.
MODE=full
source scripts/2-hmsc/helpers/run-gpu.sh
