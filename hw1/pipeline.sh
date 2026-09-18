#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Пайплайн анализа 16S рРНК метагеномных данных мышей (QIIME 2 amplicon-2024.10)
# Автор: Ярослав (prikhno.iam@phystech.edu)
# ==============================================================================

# 1. Переменные среды и пути
SIF_IMAGE="/home/STUDY/FBMF/studfbmf02_13/qiime2-amplicon-2024.10.sif"
BIND_PATH="/home/STUDY/FBMF"
Q2_HOME="/home/STUDY/FBMF/studfbmf02_13/q2home:/home/qiime2"
TMP_DIR="/home/STUDY/FBMF/studfbmf02_13/tmp"

WORKDIR="/home/STUDY/FBMF/studfbmf02_13/metagenome"
METADATA="/home/STUDY/FBMF/bioinformatics/metagenomes/mice/metadata.tsv"
ASV_TABLE="${WORKDIR}/qza/mice_ASV_table.qza"
REP_SEQS="${WORKDIR}/qza/mice_rep_seqs.qza"
TAXONOMY="${WORKDIR}/qza/mice_taxonomy.qza"

mkdir -p "${WORKDIR}/core_metrics" "${WORKDIR}/qzv" "${TMP_DIR}"

run_q2() {
  singularity exec \
    -H "${Q2_HOME}" \
    -B "${BIND_PATH}" \
    --env MPLCONFIGDIR="${TMP_DIR}",TMPDIR="${TMP_DIR}" \
    "${SIF_IMAGE}" "$@"
}

echo "=== 1. Построение филогенетического дерева ==="
if [ ! -f "${WORKDIR}/qza/rooted_tree.qza" ]; then
  run_q2 qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "${REP_SEQS}" \
    --o-alignment "${WORKDIR}/qza/aligned_rep_seqs.qza" \
    --o-masked-alignment "${WORKDIR}/qza/masked_aligned_rep_seqs.qza" \
    --o-tree "${WORKDIR}/qza/unrooted_tree.qza" \
    --o-rooted-tree "${WORKDIR}/qza/rooted_tree.qza"
fi

echo "=== 2. Расчет базовых метрик разнообразия (Rarefaction depth = 1251) ==="
run_q2 qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${WORKDIR}/qza/rooted_tree.qza" \
  --i-table "${ASV_TABLE}" \
  --p-sampling-depth 1251 \
  --m-metadata-file "${METADATA}" \
  --output-dir "${WORKDIR}/core_metrics"

echo "=== 3. Статистический анализ альфа-разнообразия (Kruskal-Wallis) ==="
run_q2 qiime diversity alpha-group-significance \
  --i-alpha-diversity "${WORKDIR}/core_metrics/faith_pd_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${WORKDIR}/core_metrics/faith_pd_significance.qzv"

run_q2 qiime diversity alpha-group-significance \
  --i-alpha-diversity "${WORKDIR}/core_metrics/shannon_vector.qza" \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${WORKDIR}/core_metrics/shannon_significance.qzv"

echo "=== 4. Статистический анализ бета-разнообразия (PERMANOVA) ==="
for metric in unweighted_unifrac weighted_unifrac bray_curtis; do
  run_q2 qiime diversity beta-group-significance \
    --i-distance-matrix "${WORKDIR}/core_metrics/${metric}_distance_matrix.qza" \
    --m-metadata-file "${METADATA}" \
    --m-metadata-column group \
    --p-pairwise \
    --o-visualization "${WORKDIR}/core_metrics/${metric}_group_significance.qzv"
done

echo "=== 5. Построение кривых альфа-разрежения (Alpha Rarefaction) ==="
run_q2 qiime diversity alpha-rarefaction \
  --i-table "${ASV_TABLE}" \
  --i-phylogeny "${WORKDIR}/qza/rooted_tree.qza" \
  --p-max-depth 4000 \
  --m-metadata-file "${METADATA}" \
  --o-visualization "${WORKDIR}/core_metrics/alpha_rarefaction.qzv"

echo "=== 6. Коллапсирование таблицы ASV до уровня родов (L6) ==="
run_q2 qiime taxa collapse \
  --i-table "${ASV_TABLE}" \
  --i-taxonomy "${TAXONOMY}" \
  --p-level 6 \
  --o-collapsed-table "${WORKDIR}/qza/mice_table_l6.qza"

echo "=== 7. Дифференциальный анализ таксонов (ANCOM-BC) ==="
run_q2 qiime composition ancombc \
  --i-table "${WORKDIR}/qza/mice_table_l6.qza" \
  --m-metadata-file "${METADATA}" \
  --p-formula 'group' \
  --o-differentials "${WORKDIR}/qza/ancombc_l6.qza"

run_q2 qiime composition tabulate \
  --i-data "${WORKDIR}/qza/ancombc_l6.qza" \
  --o-visualization "${WORKDIR}/qzv/ancombc_l6.qzv"

echo "=== Анализ завершен успешно. ==="
