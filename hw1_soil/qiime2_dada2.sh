#!/bin/bash
#SBATCH --job-name=qiime2_dada
#SBATCH --cpus-per-task=20
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf02_13/metagenome_soil/dada2_%j.log

cd /home/STUDY/FBMF/studfbmf02_13/metagenome_soil

singularity exec \
  -H /home/STUDY/FBMF/studfbmf02_13/q2home:/home/qiime2 \
  -B /home/STUDY/FBMF \
  --env MPLCONFIGDIR=/home/STUDY/FBMF/studfbmf02_13/tmp,TMPDIR=/home/STUDY/FBMF/studfbmf02_13/tmp \
  /home/STUDY/FBMF/studfbmf02_13/qiime2-amplicon-2024.10.sif \
  qiime dada2 denoise-single \
    --i-demultiplexed-seqs qza/soil_reads.qza \
    --p-trim-left 25 \
    --p-trunc-len 200 \
    --p-max-ee 3 \
    --p-n-threads 20 \
    --p-pooling-method "pseudo" \
    --p-chimera-method "consensus" \
    --p-min-fold-parent-over-abundance 4 \
    --o-table qza/soil_ASV_table.qza \
    --o-representative-sequences qza/soil_rep_seq.qza \
    --o-denoising-stats qza/soil_reads.dada2.stats.qza
