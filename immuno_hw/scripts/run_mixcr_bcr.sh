#!/bin/bash
#SBATCH --job-name=mixcr_bcr
#SBATCH --partition=AMD9554-common
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=outputs/mixcr_bcr_%j.log

source /home/STUDY/FBMF/studfbmf02_13/miniforge3/etc/profile.d/conda.sh
conda activate immune-mixcr

mixcr analyze generic-amplicon \
  --species hsa \
  --rna \
  --floating-left-alignment-boundary \
  --floating-right-alignment-boundary C \
  --assemble-clonotypes-by CDR3 \
  --threads 8 \
  data/bcr_raw/DChu_T7_pbmc_d7COVID2_2_R1.fastq.gz \
  data/bcr_raw/DChu_T7_pbmc_d7COVID2_2_R2.fastq.gz \
  outputs/bcr_mixcr/DChu_T7_pbmc_d7COVID2_2
