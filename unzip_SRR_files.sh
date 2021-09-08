#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=def-egreenbl
#SBATCH --job-name=Unzip_SSR
#SBATCH --output=output/%x-%j.out
#SBATCH --mail-user=keeganfl@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

gunzip *
