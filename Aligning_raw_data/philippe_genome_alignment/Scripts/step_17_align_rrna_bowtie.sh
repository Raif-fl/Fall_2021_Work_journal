#!/bin/bash
#SBATCH --time=00:50:00
#SBATCH --array=1-8
#SBATCH --mem=1G
#SBATCH --account=def-egreenbl
#SBATCH --job-name=bowtie_align
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

# File list should contain all of the
# Commands you wish to run.

module load bowtie2

echo "Starting task $SLURM_ARRAY_TASK_ID"
commands=$(sed -n "${SLURM_ARRAY_TASK_ID}p" step_17_bowtie_list)

# Then execute all of the commands in parrallel.
eval $commands
