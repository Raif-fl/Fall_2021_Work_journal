#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --array=1-23
#SBATCH --account=def-egreenbl
#SBATCH --job-name=Download_SSR
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

# File list should contain all of the
# Commands you wish to run.
echo "Starting task $SLURM_ARRAY_TASK_ID"
Download=$(sed -n "${SLURM_ARRAY_TASK_ID}p" step_11_command_list)

# Then start the download
eval $Download
