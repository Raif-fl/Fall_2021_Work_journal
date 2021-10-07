#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=10
#SBATCH --account=def-egreenbl
#SBATCH --job-name=Index_genome
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles GCF_000005845.2_ASM584v2_genomic.fa --sjdbGTFfile GCF_000005845.2_ASM584v2_genomic.gtf --sjdbOverhang 100 --genomeSAindexNbases 13
