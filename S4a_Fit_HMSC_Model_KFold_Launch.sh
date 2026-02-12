#!/bin/bash

# Configuration
SAMPLES=250
THIN=100
MODELS=(
  "2026-02-10_07-03-53_All_All_Atlas1_MinOccs5"
  "2026-02-10_07-03-53_All_All_Atlas2_MinOccs5"
  "2026-02-10_07-03-53_All_All_Atlas3_MinOccs5"
)

for MODEL in "${MODELS[@]}"; do
  echo "Submitting $MODEL..."
  
  # Constructing the sbatch call
  sbatch \
    --job-name="$MODEL" \
    --export=NAME="$MODEL",SAMP=$SAMPLES,THIN=$THIN \
    S4a_Fit_HMSC_Model_KFold.sh
done