#!/bin/bash
#!/bin/bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=03:00:00
#SBATCH --array=0-1
#SBATCH --mail-type=BEGIN,END,FAIL   # Alerts for start, finish, and crashes
#SBATCH --mail-user=bhr597@sund.ku.dk

# Load Conda into the shell
module restore myhmscstack 
source /opt/software/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate my_hmsc_env

MODELS=(
  "2026-01-27_14-52-46_All_All_Atlas3_MinOccs20"
  "2026-01-27_14-52-46_All_All_Atlas3_MinOccs5"
)
THIN_VAL=10

srun Rscript --verbose S4b_HPC_Post_Processing.R $THIN_VAL $MODELS[$SLURM_ARRAY_TASK_ID]