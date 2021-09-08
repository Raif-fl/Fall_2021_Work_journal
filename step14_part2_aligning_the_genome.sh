#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --array=1-4
#SBATCH --account=def-egreenbl
#SBATCH --job-name=align_genome_p2
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

STAR --genomeDir $HOME/projects/def-egreenbl/egreenbl/genomes/fly \
--runThreadN 8 \
--readFilesIn control_RPF_${SLURM_ARRAY_TASK_ID}.norrna.fastq \
--outFileNamePrefix control_RPF_${SLURM_ARRAY_TASK_I}_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard
