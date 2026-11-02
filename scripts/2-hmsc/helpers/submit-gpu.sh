set -euo pipefail
source scripts/2-hmsc/helpers/config.sh
models_file=$(mktemp)
trap 'rm -f "$models_file"' EXIT
Rscript scripts/2-hmsc/helpers/list-models.R > "$models_file"
models=()
while IFS= read -r model; do models+=("$model"); done < "$models_file"
if [[ "$MODE" == full ]]; then
  script=scripts/2-hmsc/S2a_Fit_HMSC_Model.sh
  array="0-$((CH1_CHAINS-1))"
  logs=Logs
else
  script=scripts/2-hmsc/S4a_Fit_HMSC_Model_KFold.sh
  array="1-$((CH1_NFOLDS*CH1_CHAINS))"
  logs=Logs_CV
fi
# Preflight all models before submitting the first job.
for model in "${models[@]}"; do
  base="$CH1_MODEL_OUTPUT_DIR/$model/Models"
  if [[ "$MODE" == full ]]; then
    printf -v input '%s/INITS/HPC_INIT_samples_%04d_thin_%02d_chains_%d.rds' "$base" "$CH1_SAMPLES" "$CH1_THIN" "$CH1_CHAINS"
    [[ -s "$input" ]] || { echo "Missing $input" >&2; exit 1; }
  else
    for ((i=1;i<=CH1_NFOLDS*CH1_CHAINS;i++)); do
      printf -v input '%s/Temp/temp_samples_%04d_thin_%02d_thread_%d.rds' "$base" "$CH1_SAMPLES" "$CH1_THIN" "$i"
      [[ -s "$input" ]] || { echo "Missing $input" >&2; exit 1; }
    done
  fi
  mkdir -p "$base/$logs"
done
for model in "${models[@]}"; do
  sbatch --job-name="$model" --array="$array" --export="ALL,NAME=$model" \
    --output="$CH1_MODEL_OUTPUT_DIR/$model/Models/$logs/%x_%A_%a.log" \
    --error="$CH1_MODEL_OUTPUT_DIR/$model/Models/$logs/%x_%A_%a.log" "$script"
done
