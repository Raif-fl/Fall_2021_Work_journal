#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --array=1-9
#SBATCH --account=def-egreenbl
#SBATCH --job-name=count_features
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

echo "Starting task $SLURM_ARRAY_TASK_ID"
commands=$(sed -n "${SLURM_ARRAY_TASK_ID}p" count_reads_list)

# Then start the download
eval $commands

# I am going to need to work a new file system into this, it just cannot do numbers.
