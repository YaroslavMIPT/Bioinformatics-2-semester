#!/bin/bash
set -e

PEAKS_DIR="results/peaks"
P_75ep="${PEAKS_DIR}/foxd3_75ep_peaks.narrowPeak"
P_1_2ss="${PEAKS_DIR}/foxd3_1-2ss_peaks.narrowPeak"
P_5_6ss="teacher/foxd3_5-6ss_peaks.narrowPeak"
P_14ss="${PEAKS_DIR}/foxd3_14ss_peaks.narrowPeak"

printf "%-15s\t%-8s\t%-8s\t%-15s\t%-18s\t%-14s\n" "Pair (A vs B)" "Peaks_A" "Peaks_B" "Overlap(A in B)" "Stage-spec(A\\B)" "Jaccard_Index"

compare_pair() {
    local name=$1
    local file_a=$2
    local file_b=$3

    local count_a=$(wc -l < "$file_a")
    local count_b=$(wc -l < "$file_b")
    local overlap_a=$(bedtools intersect -u -a "$file_a" -b "$file_b" | wc -l)
    local specific_a=$((count_a - overlap_a))
    local jaccard=$(bedtools jaccard -a "$file_a" -b "$file_b" | awk 'NR==2 {printf "%.5f", $3}')

    printf "%-15s\t%-8d\t%-8d\t%-15d\t%-18d\t%-14s\n" "$name" "$count_a" "$count_b" "$overlap_a" "$specific_a" "$jaccard"
}

compare_pair "75ep vs 1-2ss" "$P_75ep" "$P_1_2ss"
compare_pair "75ep vs 5-6ss" "$P_75ep" "$P_5_6ss"
compare_pair "75ep vs 14ss"  "$P_75ep" "$P_14ss"
compare_pair "1-2ss vs 5-6ss" "$P_1_2ss" "$P_5_6ss"
compare_pair "1-2ss vs 14ss"  "$P_1_2ss" "$P_14ss"
compare_pair "5-6ss vs 14ss"  "$P_5_6ss" "$P_14ss"
