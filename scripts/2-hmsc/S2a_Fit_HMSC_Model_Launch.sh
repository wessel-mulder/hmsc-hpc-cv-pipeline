#!/bin/bash
set -euo pipefail

# Configuration
SAMPLES=250
THIN=100
MODELS=(
"2026-11-12_01-29-40_Atlas1_MinOccs5_CoverageGoodAverage_EnvironmentSpace_nf2"
"2026-11-12_01-29-40_Atlas2_MinOccs5_CoverageGoodAverage_EnvironmentSpace_nf2"
"2026-11-12_01-29-40_Atlas3_MinOccs5_CoverageGoodAverage_EnvironmentSpace_nf2"
)
DRY_RUN="${DRY_RUN:-0}"

# Wall-time limits per chain. The five-factor allocation includes headroom
# above the observed runtime of less than four hours. Benchmark these runs
# before tightening limits because spatial cost need not scale linearly.
WALL_TIMES=(
  "01:30:00"
  "01:30:00"
  "01:30:00"
)

# Check every required initialisation before submitting the first job.
if [[ "$DRY_RUN" == "0" ]]; then
  for MODEL in "${MODELS[@]}"; do
    INIT_FILE=$(printf "HmscOutputs/%s/Models/INITS/HPC_INIT_samples_%04d_thin_%02d_chains_4.rds" \
      "$MODEL" "$SAMPLES" "$THIN")
    if [[ ! -f "$INIT_FILE" ]]; then
      echo "Missing initialisation: $INIT_FILE" >&2
      echo "Run scripts/2-hmsc/S2_Model_Fitting_HPC_Version.R first." >&2
      exit 1
    fi
  done
fi

for INDEX in "${!MODELS[@]}"; do
  MODEL="${MODELS[$INDEX]}"
  WALL_TIME="${WALL_TIMES[$INDEX]}"
  COMMAND=(sbatch
    --job-name="$MODEL"
    --time="$WALL_TIME"
    --export=ALL,NAME="$MODEL",SAMP="$SAMPLES",THIN="$THIN"
    scripts/2-hmsc/S2a_Fit_HMSC_Model.sh
  )
  if [[ "$DRY_RUN" == "1" ]]; then
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
  else
    echo "Submitting $MODEL with $WALL_TIME per chain..."
    "${COMMAND[@]}"
  fi
done
