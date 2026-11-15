#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 FULL_LENGTH_FASTQ OUTPUT_DIR SAMPLE_ID" >&2
    exit 2
fi

full_length_fastq=$1
output_dir=$2
sample_id=$3
flair_bin=${FLAIR_BIN:-flair}
threads=${FLAIR_THREADS:-8}
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA.}"

for input in "$full_length_fastq" "$REFERENCE_FASTA"; do
    [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 1; }
done
command -v "$flair_bin" >/dev/null 2>&1 || { echo "FLAIR executable not found: $flair_bin" >&2; exit 1; }

mkdir -p "$output_dir"
output_prefix="$output_dir/$sample_id"

"$flair_bin" align \
    -g "$REFERENCE_FASTA" \
    -r "$full_length_fastq" \
    -o "$output_prefix" \
    -t "$threads"

for output in \
    "${output_prefix}.aligned.bam" \
    "${output_prefix}.aligned.bam.bai" \
    "${output_prefix}.aligned.bed"; do
    [[ -s "$output" ]] || { echo "Expected FLAIR align output is missing or empty: $output" >&2; exit 1; }
done
