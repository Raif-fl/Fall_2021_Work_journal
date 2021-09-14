#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=10
#SBATCH --account=def-egreenbl
#SBATCH --job-name=Index_genome
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

STAR --runThreadN 8 --limitGenomeGenerateRAM 20000000000 --runMode genomeGenerate --genomeDir . --genomeFastaFiles mm10.fa --sjdbGTFfile mm10.refGene.gtf --sjdbOverhang 100 --genomeSAindexNbases 13
