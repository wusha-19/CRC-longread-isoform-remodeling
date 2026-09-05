#!/usr/bin/env python3
"""Filter the combined FLAIR atlas using the retained expression-support rule."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Optional

import pandas as pd


def read_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    return frame.rename(columns={frame.columns[0]: "isoform_id"}).assign(
        isoform_id=lambda data: data["isoform_id"].astype(str)
    )


def get_attr_value(attributes: str, key: str) -> Optional[str]:
    quoted = re.search(rf'{key}\s+"([^"]+)"', attributes)
    if quoted:
        return quoted.group(1)
    unquoted = re.search(rf"{key}=([^;]+)", attributes)
    return unquoted.group(1).strip() if unquoted else None


def filter_bed(source: Path, destination: Path, keep_ids: set[str]) -> tuple[int, int]:
    total = kept = 0
    with source.open(encoding="utf-8") as input_handle, destination.open("w", encoding="utf-8") as output_handle:
        for line in input_handle:
            if not line.strip():
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3] in keep_ids:
                output_handle.write(line)
                kept += 1
    return total, kept


def filter_fasta(source: Path, destination: Path, keep_ids: set[str]) -> tuple[int, int]:
    total = kept = 0
    write_record = False
    with source.open(encoding="utf-8") as input_handle, destination.open("w", encoding="utf-8") as output_handle:
        for line in input_handle:
            if line.startswith(">"):
                total += 1
                transcript_id = line[1:].strip().split()[0]
                write_record = transcript_id in keep_ids
                if write_record:
                    output_handle.write(line)
                    kept += 1
            elif write_record:
                output_handle.write(line)
    return total, kept


def filter_gtf(source: Path, destination: Path, keep_ids: set[str]) -> tuple[int, int, int, int]:
    total_features = kept_features = total_transcripts = kept_transcripts = 0
    with source.open(encoding="utf-8") as input_handle, destination.open("w", encoding="utf-8") as output_handle:
        for line in input_handle:
            if line.startswith("#"):
                output_handle.write(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            total_features += 1
            feature_type = fields[2]
            transcript_id = get_attr_value(fields[8], "transcript_id")
            gene_id = get_attr_value(fields[8], "gene_id")
            if transcript_id is None:
                continue
            if feature_type == "transcript":
                total_transcripts += 1
            candidate_ids = {transcript_id}
            if gene_id is not None:
                candidate_ids.add(f"{transcript_id}_{gene_id}")
            if candidate_ids & keep_ids:
                output_handle.write(line)
                kept_features += 1
                if feature_type == "transcript":
                    kept_transcripts += 1
    return total_features, kept_features, total_transcripts, kept_transcripts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--tpm", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--fa", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--min_count_per_sample", type=float, default=2)
    parser.add_argument("--min_expressed_samples", type=int, default=8)
    parser.add_argument("--min_total_count", type=float, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for path in (args.counts, args.tpm, args.bed, args.gtf, args.fa):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    counts = read_matrix(args.counts)
    tpm = read_matrix(args.tpm)
    if counts["isoform_id"].duplicated().any():
        raise ValueError("Duplicate isoform identifiers in count matrix")
    if set(counts["isoform_id"]) != set(tpm["isoform_id"]):
        raise ValueError("Count and TPM matrices contain different isoform identifiers")

    sample_columns = [column for column in counts.columns if column != "isoform_id"]
    count_matrix = counts[sample_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    n_expressed = (count_matrix >= args.min_count_per_sample).sum(axis=1)
    total_count = count_matrix.sum(axis=1)
    max_count = count_matrix.max(axis=1)
    mean_count_nonzero = count_matrix.where(count_matrix > 0).mean(axis=1).fillna(0)
    keep = (n_expressed >= args.min_expressed_samples) & (total_count >= args.min_total_count)

    support = pd.DataFrame(
        {
            "isoform_id": counts["isoform_id"],
            "n_expressed_samples": n_expressed,
            "total_count": total_count,
            "max_count": max_count,
            "mean_count_nonzero": mean_count_nonzero,
            "keep": keep,
        }
    )
    keep_ids = set(support.loc[support["keep"], "isoform_id"].astype(str))

    keep_path = args.outdir / f"{args.prefix}.keep_isoform_ids.txt"
    keep_path.write_text("".join(f"{item}\n" for item in sorted(keep_ids)), encoding="utf-8")
    support.to_csv(args.outdir / f"{args.prefix}.expression_support_summary.tsv", sep="\t", index=False)
    counts[counts["isoform_id"].isin(keep_ids)].to_csv(
        args.outdir / f"{args.prefix}.counts.tsv", sep="\t", index=False
    )
    tpm[tpm["isoform_id"].isin(keep_ids)].to_csv(
        args.outdir / f"{args.prefix}.tpm.tsv", sep="\t", index=False
    )

    bed_stats = filter_bed(args.bed, args.outdir / f"{args.prefix}.bed", keep_ids)
    fasta_stats = filter_fasta(args.fa, args.outdir / f"{args.prefix}.fa", keep_ids)
    gtf_stats = filter_gtf(args.gtf, args.outdir / f"{args.prefix}.gtf", keep_ids)

    print(
        "Expression-support rule: "
        f"count >= {args.min_count_per_sample:g} in >= {args.min_expressed_samples} samples; "
        f"total count >= {args.min_total_count:g}"
    )
    print(f"Retained isoforms: {len(keep_ids)}/{len(counts)}")
    print(f"BED retained: {bed_stats[1]}/{bed_stats[0]}")
    print(f"FASTA retained: {fasta_stats[1]}/{fasta_stats[0]}")
    print(f"GTF transcripts retained: {gtf_stats[3]}/{gtf_stats[2]}")
    if bed_stats[1] != len(keep_ids) or fasta_stats[1] != len(keep_ids) or gtf_stats[3] != len(keep_ids):
        print("WARNING: annotation transcript counts do not exactly match the retained identifier set.")


if __name__ == "__main__":
    main()
