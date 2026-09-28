#!/bin/bash
set -e

calc_frip_and_peaks() {
    local stage=$1
    local bam="data/${stage}/foxd3_${stage}_chip.nuclear.bam"
    local peaks="results/peaks/foxd3_${stage}_peaks.narrowPeak"

    if [ ! -f "$peaks" ]; then
        echo "[WAIT] Файл ${peaks} ещё не готов."
        return
    fi

    local num_peaks=$(wc -l < "$peaks")
    local total_reads=$(samtools view -c -F 0x904 "$bam")
    local reads_in_peaks=$(samtools view -c -F 0x904 -L "$peaks" "$bam")
    local frip=$(awk -v rip="$reads_in_peaks" -v tot="$total_reads" 'BEGIN {printf "%.8f", rip/tot}')

    # Топ-пик по -log10(q)
    local top_peak=$(sort -k9,9nr "$peaks" | head -n 1)
    local chrom=$(echo "$top_peak" | awk '{print $1}')
    local start=$(echo "$top_peak" | awk '{print $2}')
    local end=$(echo "$top_peak" | awk '{print $3}')
    local signalValue=$(echo "$top_peak" | awk '{print $7}')
    local pval=$(echo "$top_peak" | awk '{print $8}')
    local qval=$(echo "$top_peak" | awk '{print $9}')
    local real_q=$(awk -v q="$qval" 'BEGIN {printf "%.2e", 10^(-q)}')

    echo "=========================================="
    echo "Стадия: ${stage}"
    echo "Число пиков: ${num_peaks}"
    echo "Число выравниваний (ChIP): ${total_reads}"
    echo "Число ридов в пиках: ${reads_in_peaks}"
    echo "FRiP: ${frip}"
    echo "Топ-пик координаты: ${chrom}:${start}-${end}"
    echo "  signalValue: ${signalValue}"
    echo "  -log10(p-value): ${pval}"
    echo "  -log10(q-value): ${qval}"
    echo "  Настоящий q-value: ${real_q}"
    echo "=========================================="
}

calc_frip_and_peaks "75ep"
calc_frip_and_peaks "1-2ss"
calc_frip_and_peaks "14ss"
