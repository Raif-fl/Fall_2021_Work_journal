#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --array=1-9
#SBATCH --mem=8G
#SBATCH --account=def-egreenbl
#SBATCH --job-name=Align_rna
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

# list should contain all of the
# Commands you wish to run.

echo "Starting task $SLURM_ARRAY_TASK_ID"
commands=$(sed -n "${SLURM_ARRAY_TASK_ID}p" step_15_rna_align_list)

# Then execute all of the commands in parrallel.
eval $commands