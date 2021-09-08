#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --array=1-11
#SBATCH --account=def-egreenbl
#SBATCH --job-name=bowtie_align
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

# File list should contain all of the
# Commands you wish to run.

echo "Starting task $SLURM_ARRAY_TASK_ID"
commands=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bowtie_list)

# Then start the download
eval $commands
