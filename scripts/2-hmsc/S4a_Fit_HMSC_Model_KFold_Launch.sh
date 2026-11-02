#!/bin/bash

# Configuration
SAMPLES=250
THIN=100
MODELS=(
  "2026-07-27_13-28-41_Atlas3_MinOccs5_CoverageGoodAverage"
)

for MODEL in "${MODELS[@]}"; do
  echo "Submitting $MODEL..."
  
  # Constructing the sbatch call
  sbatch \
    --job-name="$MODEL" \
    --export=NAME="$MODEL",SAMP=$SAMPLES,THIN=$THIN \
    S4a_Fit_HMSC_Model_KFold.sh
done