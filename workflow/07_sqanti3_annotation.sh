#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 FILTERED_GTF FILTERED_FASTA FILTERED_BED OUTPUT_DIR" >&2
    exit 2
fi

isoform_gtf=$1
isoform_fasta=$2
isoform_bed=$3
output_root=$4
python_bin=${PYTHON_BIN:-python}
cpus=${SQANTI3_CPUS:-24}
chunks=${SQANTI3_CHUNKS:-1}
: "${REFERENCE_GTF:?Set REFERENCE_GTF.}"
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA.}"
: "${SQANTI3_DIR:?Set SQANTI3_DIR.}"
: "${SQANTI3_POLYA_MOTIFS:?Set SQANTI3_POLYA_MOTIFS.}"

sqanti3_qc=${SQANTI3_QC_BIN:-$SQANTI3_DIR/sqanti3_qc.py}
sqanti3_filter=${SQANTI3_FILTER_BIN:-$SQANTI3_DIR/sqanti3_filter.py}
sqanti3_rescue=${SQANTI3_RESCUE_BIN:-$SQANTI3_DIR/sqanti3_rescue.py}

for input in \
    "$isoform_gtf" "$isoform_fasta" "$isoform_bed" \
    "$REFERENCE_GTF" "$REFERENCE_FASTA" "$SQANTI3_POLYA_MOTIFS" \
    "$sqanti3_qc" "$sqanti3_filter" "$sqanti3_rescue"; do
    [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 1; }
done
command -v "$python_bin" >/dev/null 2>&1 || { echo "Python executable not found: $python_bin" >&2; exit 1; }

qc_dir="$output_root/01_QC"
filter_dir="$output_root/02_Filter_rules"
rescue_dir="$output_root/03_Rescue_automatic"
final_qc_dir="$output_root/04_Final_QC"
log_dir="$output_root/logs"
mkdir -p "$qc_dir" "$filter_dir" "$rescue_dir" "$final_qc_dir" "$log_dir"

"$python_bin" "$sqanti3_qc" \
    --isoforms "$isoform_gtf" \
    --refGTF "$REFERENCE_GTF" \
    --refFasta "$REFERENCE_FASTA" \
    --polyA_motif_list "$SQANTI3_POLYA_MOTIFS" \
    -t "$cpus" \
    -n "$chunks" \
    --force_id_ignore \
    --skipORF \
    --report html \
    -o combined \
    -d "$qc_dir" \
    >"$log_dir/01_qc.log" 2>&1

qc_class="$qc_dir/combined_classification.txt"
qc_gtf="$qc_dir/combined_corrected.gtf"
qc_fasta="$qc_dir/combined_corrected.fasta"
for output in "$qc_class" "$qc_gtf" "$qc_fasta"; do
    [[ -s "$output" ]] || { echo "SQANTI3 QC output is missing or empty: $output" >&2; exit 1; }
done

"$python_bin" "$sqanti3_filter" rules \
    --sqanti_class "$qc_class" \
    --filter_gtf "$qc_gtf" \
    --filter_isoforms "$qc_fasta" \
    -o combined_rules_filtered \
    -d "$filter_dir" \
    >"$log_dir/02_filter.log" 2>&1

filter_class=$(find "$filter_dir" -maxdepth 1 -type f -name '*RulesFilter*classification.txt' -print -quit)
filter_gtf=$(find "$filter_dir" -maxdepth 1 -type f -name '*.filtered.gtf' -print -quit)
[[ -n "$filter_class" && -s "$filter_class" ]] || { echo "SQANTI3 filtered classification was not produced" >&2; exit 1; }
[[ -n "$filter_gtf" && -s "$filter_gtf" ]] || { echo "SQANTI3 filtered GTF was not produced" >&2; exit 1; }

"$python_bin" "$sqanti3_rescue" rules \
    --filter_class "$filter_class" \
    --refGTF "$REFERENCE_GTF" \
    --refFasta "$REFERENCE_FASTA" \
    --rescue_isoforms "$qc_fasta" \
    --rescue_gtf "$filter_gtf" \
    --refClassif "$qc_class" \
    --mode automatic \
    -e none \
    -o combined_rescued \
    -d "$rescue_dir" \
    --cpus "$cpus" \
    >"$log_dir/03_rescue.log" 2>&1

rescue_gtf=$(find "$rescue_dir" -maxdepth 1 -type f -name '*rescued.gtf' -print -quit)
[[ -n "$rescue_gtf" && -s "$rescue_gtf" ]] || { echo "SQANTI3 rescued GTF was not produced" >&2; exit 1; }

"$python_bin" "$sqanti3_qc" \
    --isoforms "$rescue_gtf" \
    --refGTF "$REFERENCE_GTF" \
    --refFasta "$REFERENCE_FASTA" \
    --polyA_motif_list "$SQANTI3_POLYA_MOTIFS" \
    -t "$cpus" \
    -n "$chunks" \
    --force_id_ignore \
    --skipORF \
    --report html \
    -o final_rescued \
    -d "$final_qc_dir" \
    >"$log_dir/04_final_qc.log" 2>&1

final_class="$final_qc_dir/final_rescued_classification.txt"
[[ -s "$final_class" ]] || { echo "Final SQANTI3 classification is missing or empty: $final_class" >&2; exit 1; }
