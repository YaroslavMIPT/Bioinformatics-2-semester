#!/bin/bash
#SBATCH --job-name=q2_downstream
#SBATCH --partition=IXG6154-AI-common
#SBATCH --cpus-per-task=8
#SBATCH --mem=60gb
#SBATCH --time=02:00:00
#SBATCH --output=/home/STUDY/FBMF/studfbmf02_13/metagenome_soil/downstream_%j.log

set -e

cd /home/STUDY/FBMF/studfbmf02_13/metagenome_soil

SIF="/home/STUDY/FBMF/studfbmf02_13/qiime2-amplicon-2024.10.sif"
BIND_OPTS="-H /home/STUDY/FBMF/studfbmf02_13/q2home:/home/qiime2 -B /home/STUDY/FBMF --env MPLCONFIGDIR=/home/STUDY/FBMF/studfbmf02_13/tmp,TMPDIR=/home/STUDY/FBMF/studfbmf02_13/tmp"
METADATA="/home/STUDY/FBMF/bioinformatics/metagenomes/soil_hw/soil_metadata_full.tsv"
CLASSIFIER="/home/STUDY/FBMF/bioinformatics/metagenomes/soil_hw/2022.10.backbone.v4.nb.sklearn-1.4.2.qza"

if [ -f qza/rooted-tree.qza ]; then
  echo "=== Шаг 1: Филогенетическое дерево уже готово, пропускаем ==="
else
  echo "=== 1. Построение филогенетического дерева ==="
  singularity exec ${BIND_OPTS} ${SIF} qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences qza/soil_rep_seq.qza \
    --o-alignment qza/aligned-rep-seqs.qza \
    --o-masked-alignment qza/masked-aligned-rep-seqs.qza \
    --o-tree qza/unrooted-tree.qza \
    --o-rooted-tree qza/rooted-tree.qza \
    --p-n-threads 8
fi

echo "=== 2. Таксономическая классификация (sklearn, 4 потока) ==="
singularity exec ${BIND_OPTS} ${SIF} qiime feature-classifier classify-sklearn \
  --i-classifier ${CLASSIFIER} \
  --i-reads qza/soil_rep_seq.qza \
  --p-reads-per-batch 10000 \
  --p-n-jobs 4 \
  --o-classification qza/taxonomy.qza

singularity exec ${BIND_OPTS} ${SIF} qiime metadata tabulate \
  --m-input-file qza/taxonomy.qza \
  --o-visualization qzv/taxonomy.qzv

echo "=== 3. Построение графиков таксономии (Barplot) ==="
singularity exec ${BIND_OPTS} ${SIF} qiime taxa barplot \
  --i-table qza/soil_ASV_table.qza \
  --i-taxonomy qza/taxonomy.qza \
  --m-metadata-file ${METADATA} \
  --o-visualization qzv/taxa-bar-plots.qzv

echo "=== 4. Расчет разнообразия (sampling-depth = 4809) ==="
rm -rf core-metrics-results
singularity exec ${BIND_OPTS} ${SIF} qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qza/rooted-tree.qza \
  --i-table qza/soil_ASV_table.qza \
  --p-sampling-depth 4809 \
  --m-metadata-file ${METADATA} \
  --output-dir core-metrics-results

echo "=== Все вычисления успешно завершены! ==="
