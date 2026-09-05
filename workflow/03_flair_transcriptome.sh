#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 LONGREAD_BAM STAR_SJ_OUT_TAB OUTPUT_DIR SAMPLE_ID" >&2
    exit 2
fi

longread_bam=$1
star_junctions=$2
output_dir=$3
sample_id=$4
flair_bin=${FLAIR_BIN:-flair}
samtools_bin=${SAMTOOLS_BIN:-samtools}
threads=${FLAIR_THREADS:-8}
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA.}"
: "${REFERENCE_GTF:?Set REFERENCE_GTF.}"

for input in "$longread_bam" "$star_junctions" "$REFERENCE_FASTA" "$REFERENCE_GTF"; do
    [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 1; }
done
command -v "$flair_bin" >/dev/null 2>&1 || { echo "FLAIR executable not found: $flair_bin" >&2; exit 1; }
command -v "$samtools_bin" >/dev/null 2>&1 || { echo "samtools executable not found: $samtools_bin" >&2; exit 1; }

if [[ ! -s "${longread_bam}.bai" ]]; then
    "$samtools_bin" index -@ "$threads" "$longread_bam"
fi

mkdir -p "$output_dir"
output_prefix="$output_dir/$sample_id"

"$flair_bin" transcriptome \
    -g "$REFERENCE_FASTA" \
    -f "$REFERENCE_GTF" \
    --junction_tab "$star_junctions" \
    -b "$longread_bam" \
    -o "$output_prefix" \
    -t "$threads" \
    -s 5 \
    --parallelmode byregion \
    --filter nosubset \
    --stringent \
    --check_splice

for output in \
    "${output_prefix}.isoforms.bed" \
    "${output_prefix}.isoforms.gtf" \
    "${output_prefix}.isoforms.fa"; do
    [[ -s "$output" ]] || { echo "Expected FLAIR output is missing or empty: $output" >&2; exit 1; }
done
