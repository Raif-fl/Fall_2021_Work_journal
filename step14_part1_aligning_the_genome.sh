#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --array=1-7
#SBATCH --mem=8G
#SBATCH --account=def-egreenbl
#SBATCH --job-name=align_genome_p1
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

STAR --genomeDir $HOME/projects/def-egreenbl/egreenbl/genomes/fly \
--runThreadN 8 \
--readFilesIn Fmr1_RPF_${SLURM_ARRAY_TASK_ID}.norrna.fastq \
--outFileNamePrefix Fmr1_RPF_${SLURM_ARRAY_TASK_ID}_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard & \
