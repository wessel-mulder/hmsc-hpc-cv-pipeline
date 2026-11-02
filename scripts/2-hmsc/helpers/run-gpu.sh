# Called by the full-fit and CV Slurm jobs, with MODE=full or MODE=cv.
set -euo pipefail
source scripts/2-hmsc/helpers/config.sh
: "${NAME:?Submit with the full model name in NAME}"
: "${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array task}"
: "${SLURM_JOB_ID:?Missing Slurm job ID}"
[[ "$NAME" =~ ^${CH1_RUN_DATE}(_[0-9]{2}-[0-9]{2}-[0-9]{2})?_Atlas[123]_MinOccs5_CoverageGoodAverage$ ]] || { echo 'Model name disagrees with run configuration' >&2; exit 1; }
[[ "$SLURM_ARRAY_TASK_ID" =~ ^[0-9]+$ ]] || exit 1
base="$CH1_MODEL_OUTPUT_DIR/$NAME/Models"
if [[ "$MODE" == full ]]; then
  (( SLURM_ARRAY_TASK_ID < CH1_CHAINS )) || exit 1
  printf -v input '%s/INITS/HPC_INIT_samples_%04d_thin_%02d_chains_%d.rds' "$base" "$CH1_SAMPLES" "$CH1_THIN" "$CH1_CHAINS"
  printf -v output '%s/Sampled/HPC_samples_%04d_thin_%02d_chain_%d.rds' "$base" "$CH1_SAMPLES" "$CH1_THIN" "$SLURM_ARRAY_TASK_ID"
  logs="$base/Logs"
else
  (( SLURM_ARRAY_TASK_ID >= 1 && SLURM_ARRAY_TASK_ID <= CH1_NFOLDS*CH1_CHAINS )) || exit 1
  printf -v input '%s/Temp/temp_samples_%04d_thin_%02d_thread_%d.rds' "$base" "$CH1_SAMPLES" "$CH1_THIN" "$SLURM_ARRAY_TASK_ID"
  printf -v output '%s/Temp/Sampled_HPC_samples_%04d_thin_%02d_thread_%d.rds' "$base" "$CH1_SAMPLES" "$CH1_THIN" "$SLURM_ARRAY_TASK_ID"
  logs="$base/Logs_CV"
fi
[[ -s "$input" ]] || { echo "Missing initialization: $input" >&2; exit 1; }
[[ ! -e "$output" || "$CH1_OVERWRITE" == true ]] || { echo "Output exists: $output" >&2; exit 1; }
mkdir -p "$logs"
module restore myhmscstack
source /opt/software/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate "$CH1_CONDA_ENV"
python3 - <<'PY'
import sys, tensorflow as tf
import importlib.metadata as metadata
print('Python:', sys.version, 'TensorFlow:', tf.__version__, 'hmsc:', metadata.version('hmsc'))
print('GPU:', tf.config.list_physical_devices('GPU'))
if not tf.config.list_physical_devices('GPU'):
    raise RuntimeError('No TensorFlow GPU available')
PY
metadata_file="$logs/run_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
{
  date; hostname
  env | sort | sed -n '/^CH1_/p'
  echo "Input: $input"; echo "Output: $output"
  md5sum "$input"
  python3 -m pip show hmsc tensorflow
} > "$metadata_file"
cp scripts/2-hmsc/run-config.tsv "$logs/config_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tsv"
cp "$0" "$logs/launcher_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.sh"
cp scripts/2-hmsc/helpers/run-gpu.sh "$logs/helper_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.sh"
nvidia-smi --query-gpu=timestamp,index,name,utilization.gpu,memory.used,memory.total --format=csv -l 10 > "$logs/gpu_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.csv" &
monitor_pid=$!
trap 'kill "$monitor_pid" 2>/dev/null || true' EXIT
# Python writes to a temporary name; failed jobs cannot leave apparently complete chains.
partial="${output}.partial.${SLURM_JOB_ID}"
sampler_args=(--input "$input" --output "$partial" --samples "$CH1_SAMPLES" --transient "$CH1_TRANSIENT" --thin "$CH1_THIN" --verbose 100)
if [[ "$MODE" == full ]]; then sampler_args+=(--chain "$SLURM_ARRAY_TASK_ID"); fi
srun "$(command -v python3)" -m hmsc.run_gibbs_sampler "${sampler_args[@]}"
[[ -s "$partial" ]] || { echo 'Sampler produced no output' >&2; exit 1; }
mv "$partial" "$output"
