#!/bin/bash
set -euo pipefail

# Configuration
SAMPLES=250
THIN=100
MODEL_STAMP="2026-11-11_02-19-42"
MODEL_BASE="${MODEL_STAMP}_Atlas1_MinOccs5_CoverageGoodAverage"
# Preview all submissions locally with DRY_RUN=1 bash this-script.sh.
DRY_RUN="${DRY_RUN:-0}"
MODELS=(
  "${MODEL_BASE}_EnvironmentOnly"
  "${MODEL_BASE}_EnvironmentSpace_nf1"
  "${MODEL_BASE}_EnvironmentSpace_nf2"
  "${MODEL_BASE}_EnvironmentSpace_nf3"
  "${MODEL_BASE}_EnvironmentSpace_nf4"
  "${MODEL_BASE}_EnvironmentSpace_nf5"
)
# Wall-time limits per chain. The five-factor allocation includes headroom
# above the observed runtime of less than four hours. Benchmark these runs
# before tightening limits because spatial cost need not scale linearly.
WALL_TIMES=(
  "01:00:00"
  "02:00:00"
  "02:30:00"
  "03:00:00"
  "03:30:00"
  "05:00:00"
)

if [[ "$DRY_RUN" != "0" && "$DRY_RUN" != "1" ]]; then
  echo "DRY_RUN must be 0 or 1." >&2
  exit 1
fi

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
