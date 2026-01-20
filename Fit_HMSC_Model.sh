#!/bin/bash
#SBATCH --job-name=hmsc-trial
#SBATCH --partition=gpuqueue
#SBATCH --cpus-per-gpu=16
#SBATCH --mem-per-cpu=1GB
#SBATCH --gpus-per-node=1
#SBATCH --time=00:00:05
#SBATCH --array=0-3
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=bhr597@sund.ku.dk

# Load Conda into the shell
module restore myhmscstack
source /opt/software/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate my_hmsc_env

### DATA PATH
Area='Hmsc Outputs'
model_name=$NAME
data_path=$(printf "%s/%s" $Area $model_name)

### WHICH MODEL 
samples=$SAMP
thin=$THIN
chains=4
verbose=100
transient=100000
screen_text=$(printf "#######\n Model Run starting\n Samples: %.4d\n Thining: %.2d Transient steps: %.d\n######\n" $samples $thin $transient)

input_path=$data_path/$(printf "Models/INIT/HPC_INIT_samples_%.4d_thin_%.2d_chains_%.1d.rds" $samples $thin $chains)
output_path=$data_path/$(printf "Models/Sampled/HPC_samples_%.4d_thin_%.2d_chain_%.1d.rds" $samples $thin $SLURM_ARRAY_TASK_ID)

echo $screen_text
echo $data_path
echo $input_path
echo $output_path

srun python Gibs_sampling_tensor_flow.py --samples $samples --transient $transient --thin $thin --verbose $verbose --input $input_path --output $output_path --chain $SLURM_ARRAY_TASK_ID