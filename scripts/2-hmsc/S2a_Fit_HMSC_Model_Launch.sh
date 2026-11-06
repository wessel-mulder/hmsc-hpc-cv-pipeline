#!/bin/bash

# Configuration
SAMPLES=250
THIN=100
MODELS=(
  "2026-11-06_06-55-21_Atlas1_MinOccs5_CoverageGoodAverage"
  "2026-11-06_06-55-21_Atlas2_MinOccs5_CoverageGoodAverage"
  "2026-11-06_06-55-21_Atlas3_MinOccs5_CoverageGoodAverage"

)

for MODEL in "${MODELS[@]}"; do
  echo "Submitting $MODEL..."
  
  # Constructing the sbatch call
  sbatch \
    --job-name="$MODEL" \
    --export=NAME="$MODEL",SAMP=$SAMPLES,THIN=$THIN \
    scripts/2-hmsc/S2a_Fit_HMSC_Model.sh
done
