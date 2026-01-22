#!/bin/bash
#SBATCH --partition=gpuqueue
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1GB
#SBATCH --gpus-per-node=1
#SBATCH --time=02:00:00
#SBATCH --array=1-20
#SBATCH --output=HmscOutputs/%x/Models/Logs_CV/%x_%A_%a.log
#SBATCH --error=HmscOutputs/%x/Models/Logs_CV/%x_%A_%a.log


# Load Conda into the shell
module restore myhmscstack 
source /opt/software/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate my_hmsc_env

echo "=== TensorFlow GPU sanity check ==="

python3 - << 'EOF'
import tensorflow as tf

print("TF version:", tf.__version__)
print("Built with CUDA:", tf.test.is_built_with_cuda())
print("GPUs:", tf.config.list_physical_devices("GPU"))
EOF

echo "=== End TF check ==="
echo "=== MODULE LIST ==="
module list || { echo "ERROR: module command not available"; exit 1; }
echo "==================="

### DATA PATH
Area='HmscOutputs'
model_name=$NAME
data_path=$(printf "%s/%s" $Area $model_name)

### WHICH MODEL 
samples=$SAMP
thin=$THIN
chains=1
verbose=100
transient=100000
screen_text=$(printf "#######\n Model Run starting\n Samples: %.4d\n Thining: %.2d Transient steps: %.d\n######\n" $samples $thin $transient)

input_path=$data_path/$(printf "Models/Temp/temp_samples_%.4d_thin_%.2d_thread_%.1d.rds" $samples $thin $SLURM_ARRAY_TASK_ID)
output_path=$data_path/$(printf "Models/Temp/Sampled_HPC_samples_%.4d_thin_%.2d_thread_%.1d.rds" $samples $thin $SLURM_ARRAY_TASK_ID)

# Create a logs directory in the model output
logs_dir=$data_path/Models/Logs_CV

# Start GPU logging in the background
gpu_log="$logs_dir/gpu_usage_JOBID_${SLURM_JOB_ID}_TASK_${SLURM_ARRAY_TASK_ID}.csv"

nvidia-smi --query-gpu=timestamp,index,name,utilization.gpu,memory.used,memory.total --format=csv -l 10 > $logs_dir/gpu_usage_${SLURM_ARRAY_TASK_ID}.log & NVIDIA_MONITOR_PID=$!

echo "GPU logging started: $gpu_log (PID=$NVIDIA_MONITOR_PID)"

# Only do this once per array job (save metadata)
if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then

  # Copy the bash script to the logs folder
  script_copy=$logs_dir/$(printf \
    "Fit_HMSC_Model_SAMP_%.4d_THIN_%.2d_JOBID_%s.sh" \
    $samples $thin $SLURM_JOB_ID)

  cp "$0" "$script_copy"

  # Save run parameters
  params_file=$logs_dir/$(printf \
    "run_params_SAMP_%.4d_THIN_%.2d_JOBID_%s.txt" \
    $samples $thin $SLURM_JOB_ID)

cat > $params_file << EOF
=== HMSC Model Fit Parameters ===
Samples: $samples
Thin: $thin
Chains: $chains
Transient: $transient
Job ID: $SLURM_JOB_ID
Submission Time: $(date)
Node: $(hostname)
Input: $input_path
Output: $output_path
=============================
EOF

echo "Run metadata saved to: $logs_dir"

fi

# Run sampling for all array tasks
echo "Run metadata saved to: $logs_dir"


srun python3 -m hmsc.run_gibbs_sampler \
  --input $input_path \
  --output $output_path \
  --samples $samples \
  --transient $transient \
  --thin $thin \
  --verbose $verbose 
