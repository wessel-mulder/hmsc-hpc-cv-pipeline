#!/bin/bash

# Configuration
SAMPLES=250
THIN=100
MODELS=(
  "2026-01-27_14-52-46_All_All_Atlas3_MinOccs20"
  "2026-01-27_14-52-46_All_All_Atlas3_MinOccs5"
)

for MODEL in "${MODELS[@]}"; do
  echo "Submitting $MODEL..."
  
  # Constructing the sbatch call
  sbatch \
    --job-name="$MODEL" \
    --export=NAME="$MODEL",SAMP=$SAMPLES,THIN=$THIN \
    Fit_HMSC_Model_KFold.sh
done