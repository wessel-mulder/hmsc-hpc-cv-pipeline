#!/bin/bash
#SBATCH --partition=gpuqueue
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1GB
#SBATCH --gpus-per-node=1
#SBATCH --time=03:00:00
#SBATCH --array=1-40
# Submit via S4a_Fit_HMSC_Model_KFold_Launch.sh from the repository root.
MODE=cv
source scripts/2-hmsc/helpers/run-gpu.sh
