#!/bin/bash
#SBATCH --job-name=Partition_Creation
#SBATCH  -M ukko
#SBATCH --partition=short 
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=R-%x.%j.out

module purge
module load R

srun Rscript --vanilla S4_Partition_Creation.R