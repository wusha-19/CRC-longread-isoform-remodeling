#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 INPUT_FASTQ OUTPUT_DIR SAMPLE_ID" >&2
    exit 2
fi

input_fastq=$1
output_dir=$2
sample_id=$3
pychopper_bin=${PYCHOPPER_BIN:-pychopper}
threads=${PYCHOPPER_THREADS:-20}

[[ -s "$input_fastq" ]] || { echo "Missing or empty input: $input_fastq" >&2; exit 1; }
command -v "$pychopper_bin" >/dev/null 2>&1 || {
    echo "Pychopper executable not found: $pychopper_bin" >&2
    exit 1
}

mkdir -p "$output_dir"
report_file="$output_dir/${sample_id}_report.pdf"
unclassified_file="$output_dir/${sample_id}_unclassified.fq"
rescued_file="$output_dir/${sample_id}_rescued.fq"
output_fastq="$output_dir/${sample_id}_full_length_output.fq"
log_file="$output_dir/pychopper_${sample_id}.log"

if [[ -e "$output_fastq" && "${ALLOW_OVERWRITE:-0}" != 1 ]]; then
    echo "Output already exists; set ALLOW_OVERWRITE=1 to replace it: $output_fastq" >&2
    exit 1
fi

"$pychopper_bin" \
    -t "$threads" \
    -r "$report_file" \
    -u "$unclassified_file" \
    -w "$rescued_file" \
    "$input_fastq" \
    "$output_fastq" \
    2>&1 | tee "$log_file"

[[ -s "$output_fastq" ]] || { echo "Pychopper output is missing or empty" >&2; exit 1; }
