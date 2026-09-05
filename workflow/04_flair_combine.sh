#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 COMBINE_MANIFEST OUTPUT_PREFIX" >&2
    exit 2
fi
manifest=$1
output_prefix=$2
flair_bin=${FLAIR_BIN:-flair}

[[ -s "$manifest" ]] || { echo "Missing or empty manifest: $manifest" >&2; exit 1; }
command -v "$flair_bin" >/dev/null 2>&1 || { echo "FLAIR executable not found: $flair_bin" >&2; exit 1; }
mkdir -p "$(dirname "$output_prefix")"

"$flair_bin" combine \
    -m "$manifest" \
    -o "$output_prefix" \
    -c

if [[ ! -s "${output_prefix}.isoforms.bed" && ! -s "${output_prefix}.bed" ]]; then
    echo "Combined FLAIR BED is missing or empty for prefix: $output_prefix" >&2
    exit 1
fi
