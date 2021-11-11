#!/bin/bash

#SBATCH --account=ACCOUNT_HERE
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=EMAIL_HERE
#SBATCH --mail-type=ALL

#SBATCH --job-name=NAME_HERE
#SBATCH --output=DIRECTORY_HERE/%x-%j.out
#SBATCH --error=DIRECTORY_HERE/%x-%j.err

#SBATCH --array=0-99%50

module load python/3.7
module load scipy-stack

python3 example_cluster.py $SLURM_ARRAY_TASK_ID
