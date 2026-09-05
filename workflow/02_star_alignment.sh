#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 R1_FASTQ_GZ R2_FASTQ_GZ OUTPUT_DIR SAMPLE_ID" >&2
    exit 2
fi

r1=$1
r2=$2
output_dir=$3
sample_id=$4
star_bin=${STAR_BIN:-STAR}
threads=${STAR_THREADS:-8}
: "${STAR_GENOME_DIR:?Set STAR_GENOME_DIR to the STAR index built for the declared reference.}"

for input in "$r1" "$r2"; do
    [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 1; }
done
[[ -d "$STAR_GENOME_DIR" ]] || { echo "STAR index not found: $STAR_GENOME_DIR" >&2; exit 1; }
command -v "$star_bin" >/dev/null 2>&1 || { echo "STAR executable not found: $star_bin" >&2; exit 1; }

mkdir -p "$output_dir"
prefix="$output_dir/${sample_id}_"
tmp_dir="$output_dir/tmp_${sample_id}_$$"
bam="${prefix}Aligned.sortedByCoord.out.bam"
sj="${prefix}SJ.out.tab"
log_final="${prefix}Log.final.out"

if [[ -e "$bam" && "${ALLOW_OVERWRITE:-0}" != 1 ]]; then
    echo "Output already exists; set ALLOW_OVERWRITE=1 to replace it: $bam" >&2
    exit 1
fi

"$star_bin" \
    --runThreadN "$threads" \
    --genomeDir "$STAR_GENOME_DIR" \
    --readFilesIn "$r1" "$r2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$prefix" \
    --outTmpDir "$tmp_dir" \
    --outSAMtype BAM SortedByCoordinate

[[ -s "$bam" ]] || { echo "STAR BAM is missing or empty: $bam" >&2; exit 1; }
[[ -s "$sj" ]] || { echo "STAR junction table is missing or empty: $sj" >&2; exit 1; }
[[ -s "$log_final" ]] || { echo "STAR final log is missing or empty: $log_final" >&2; exit 1; }
