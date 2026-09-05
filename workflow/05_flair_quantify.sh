#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 READS_MANIFEST COMBINED_FASTA COMBINED_BED OUTPUT_PREFIX TEMP_DIR" >&2
    exit 2
fi

reads_manifest=$1
combined_fasta=$2
combined_bed=$3
output_prefix=$4
temp_dir=$5
flair_bin=${FLAIR_BIN:-flair}
threads=${FLAIR_QUANT_THREADS:-16}

for input in "$reads_manifest" "$combined_fasta" "$combined_bed"; do
    [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 1; }
done
if awk -F '\t' 'NF != 4 {exit 1}' "$reads_manifest"; then
    :
else
    echo "The FLAIR reads manifest must contain four tab-separated columns on every row." >&2
    exit 1
fi
command -v "$flair_bin" >/dev/null 2>&1 || { echo "FLAIR executable not found: $flair_bin" >&2; exit 1; }

mkdir -p "$(dirname "$output_prefix")" "$temp_dir"

"$flair_bin" quantify \
    -r "$reads_manifest" \
    -i "$combined_fasta" \
    --isoform_bed "$combined_bed" \
    -o "$output_prefix" \
    -t "$threads" \
    --temp_dir "$temp_dir" \
    --sample_id_only \
    --tpm \
    --stringent \
    --check_splice

[[ -s "${output_prefix}.counts.tsv" ]] || {
    echo "FLAIR count matrix is missing or empty: ${output_prefix}.counts.tsv" >&2
    exit 1
}
[[ -s "${output_prefix}.tpm.tsv" ]] || {
    echo "FLAIR TPM matrix is missing or empty: ${output_prefix}.tpm.tsv" >&2
    exit 1
}
