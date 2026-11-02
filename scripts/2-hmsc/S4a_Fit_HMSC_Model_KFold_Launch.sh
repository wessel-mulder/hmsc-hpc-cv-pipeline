#!/bin/bash

# Configuration
SAMPLES=250
THIN=100
MODELS=(
  "2026-03-13_06-58-56_Atlas1_MinOccs5_CoverageGoodAverage"
  "2026-03-13_06-58-56_Atlas2_MinOccs5_CoverageGoodAverage"
)

for MODEL in "${MODELS[@]}"; do
  echo "Submitting $MODEL..."
  
  # Constructing the sbatch call
  sbatch \
    --job-name="$MODEL" \
    --export=NAME="$MODEL",SAMP=$SAMPLES,THIN=$THIN \
    S4a_Fit_HMSC_Model_KFold.sh
done