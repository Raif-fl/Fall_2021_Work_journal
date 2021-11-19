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


STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.toplevel.fa --sjdbGTFfile Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.51.gtf --sjdbOverhang 100 --genomeSAindexNbases 13
