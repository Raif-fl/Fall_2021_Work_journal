#!/bin/bash
#SBATCH --time=01:10:00
#SBATCH --array=1-9
#SBATCH --account=def-egreenbl
#SBATCH --mem=1G
#SBATCH --job-name=clip_ends
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

# Command list should contain all of the
# Commands you wish to run.

module load fastx-toolkit

echo "Starting task $SLURM_ARRAY_TASK_ID"
commands=$(sed -n "${SLURM_ARRAY_TASK_ID}p" step_16_clip_list)

# Then execute all of the commands in parrallel.
eval $commands
