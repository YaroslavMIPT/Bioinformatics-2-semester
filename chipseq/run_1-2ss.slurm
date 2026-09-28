#!/bin/bash
#SBATCH --job-name=chip_1-2ss
#SBATCH --output=logs/chip_1-2ss_%j.out
#SBATCH --error=logs/chip_1-2ss_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=01:30:00
#SBATCH --partition IXG6154-AI-common

mkdir -p results/signal results/peaks results/figures results/tables logs

echo "=== 1. Построение BigWig для ChIP ==="
bamCoverage \
  -b data/1-2ss/foxd3_1-2ss_chip.nuclear.bam \
  -o results/signal/foxd3_1-2ss_chip.CPM.bw \
  --normalizeUsing CPM \
  --binSize 10 \
  --numberOfProcessors 8

echo "=== 2. Построение BigWig для Input ==="
bamCoverage \
  -b data/1-2ss/foxd3_1-2ss_input.nuclear.bam \
  -o results/signal/foxd3_1-2ss_input.CPM.bw \
  --normalizeUsing CPM \
  --binSize 10 \
  --numberOfProcessors 8

echo "=== 3. Сравнение ChIP vs Input (log2) ==="
bigwigCompare \
  -b1 results/signal/foxd3_1-2ss_chip.CPM.bw \
  -b2 results/signal/foxd3_1-2ss_input.CPM.bw \
  --operation log2 \
  --pseudocount 1 \
  --binSize 10 \
  -o results/signal/foxd3_1-2ss.log2_IP_Input.bw

echo "=== 4. Поиск пиков (MACS3) ==="
macs3 callpeak \
  -t data/1-2ss/foxd3_1-2ss_chip.nuclear.bam \
  -c data/1-2ss/foxd3_1-2ss_input.nuclear.bam \
  -f BAMPE \
  -g 1360000000 \
  -n foxd3_1-2ss \
  -q 0.01 \
  --outdir results/peaks

echo "=== 5. Сбор матрицы для Heatmap/Profile ==="
computeMatrix reference-point \
  --referencePoint center \
  -R results/peaks/foxd3_1-2ss_summits.bed \
  -S results/signal/foxd3_1-2ss_chip.CPM.bw results/signal/foxd3_1-2ss_input.CPM.bw \
  -b 2000 -a 2000 \
  --binSize 20 \
  -o results/figures/foxd3_1-2ss.matrix.gz \
  -p 8

echo "=== 6. Отрисовка графиков ==="
plotHeatmap \
  -m results/figures/foxd3_1-2ss.matrix.gz \
  -out results/figures/foxd3_1-2ss.heatmap.pdf \
  --samplesLabel "Foxd3 ChIP" "Input" \
  --refPointLabel "summit"

plotProfile \
  -m results/figures/foxd3_1-2ss.matrix.gz \
  -out results/figures/foxd3_1-2ss.profile.pdf \
  --samplesLabel "Foxd3 ChIP" "Input" \
  --refPointLabel "summit"

echo "=== Стадия 1-2ss успешно завершена ==="
