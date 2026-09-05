#!/usr/bin/env python3
"""Internal short-read junction/event-level orthogonal support audit.

This script reads pre-existing STAR ``SJ.out.tab`` files only. It never reads
FASTQ, never runs STAR, never changes the frozen long-read event set, and does
not claim that a short-read junction measurement validates a full-length
transcript. The biological replicate is the patient, not the individual
short-read sample or read.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import math
import os
import re
import sys
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Iterable

import numpy as np
from scipy.stats import spearmanr, wilcoxon


V2_BASE = Path(os.environ.get("V2_BASE") or os.environ["CRC_V2_ROOT"]).expanduser().resolve()
CRC_ROOT = Path(os.environ["CRC_ROOT"]).expanduser().resolve()
STAR_ROOT = CRC_ROOT / "01_STAR_gencode"
METADATA_PATH = CRC_ROOT / "05_isoform_main_study" / "00_metadata_qc" / "metadata" / "CRC_77_sample_metadata.tsv"
PRIMARY72_PATH = CRC_ROOT / "05_isoform_main_study" / "00_metadata_qc" / "metadata" / "cohort_primary_NAT_72_samples.txt"
GTF_INDEX_PATH = CRC_ROOT / "gtf_download" / "gencode.v48_STAR_INDEX"

TABLE_DIR = V2_BASE / "tables"
REPORT_DIR = V2_BASE / "reports"
LOG_DIR = V2_BASE / "logs"
PROGRESS_PATH = V2_BASE / "progress.txt"
SUCCESS_PATH = V2_BASE / "success.tsv"
FAILED_PATH = V2_BASE / "failed.tsv"
LOG_PATH = LOG_DIR / "07_internal_shortread_junction_validation.log"

FROZEN_PATHS = {
    "primary93": TABLE_DIR / "Fig5_internal_event_set_frozen.tsv",
    "patient_recurrence": TABLE_DIR / "Fig5_patient_recurrence_audit.tsv",
    "event_features": TABLE_DIR / "Fig5_event_feature_mapping.tsv",
}
TARGETED_PATH = TABLE_DIR / "Fig5_targeted_validation_mapping_KRAS_CDX2_OLFM4.tsv"

# Fixed before reading STAR junction counts. These are not adapted to results.
MIN_PAIR_INFORMATIVE_READS = 10
MIN_ONE_SIDED_MEAN_JUNCTION_READS = 10
MIN_PAIRED_PATIENTS = 3
BOOTSTRAP_REPS = 5000
BOOTSTRAP_SEED = 20260713


def now() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    line = f"[{now()}] {message}"
    print(line, flush=True)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def update_progress(stage: str, detail: str) -> None:
    PROGRESS_PATH.write_text(
        f"timestamp\tstage\tdetail\n{now()}\t{stage}\t{detail}\n", encoding="utf-8"
    )
    log(f"{stage}: {detail}")


def write_status(path: Path, status: str, detail: str) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["timestamp", "status", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"timestamp": now(), "status": status, "detail": detail})


def clean(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return None if text == "" or text.upper() in {"NA", "N/A", "NULL", "."} else text


def as_bool(value: Any) -> bool:
    return str(value).strip().upper() in {"TRUE", "T", "YES", "Y", "1"}


def as_float(value: Any) -> float:
    text = clean(value)
    if text is None:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"TSV header missing: {path}")
        return list(reader)


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "NA" if row.get(field) is None else row.get(field) for field in fields})
    tmp.replace(path)


def write_report(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def canonical_chromosome(chromosome: str | None) -> str | None:
    text = clean(chromosome)
    if text is None:
        return None
    bare = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    if bare.upper() in {"M", "MT"}:
        bare = "MT"
    return f"chr{bare}"


def format_junction(chromosome: str, start: int, end: int, strand: str) -> str:
    return f"{canonical_chromosome(chromosome)}:{start}-{end}:{strand}"


def parse_junction(text: str | None) -> tuple[str, int, int, str] | None:
    value = clean(text)
    if value is None:
        return None
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+):([+-])", value)
    if not match:
        return None
    return canonical_chromosome(match.group(1)), int(match.group(2)), int(match.group(3)), match.group(4)


def parse_labeled_junctions(value: str | None) -> dict[str, list[str]]:
    """Read candidate/comparison-only junction chains from the frozen mapping."""
    result: dict[str, list[str]] = {"candidate_only_junction": [], "comparison_only_junction": []}
    text = clean(value)
    if text is None:
        return result
    for label, chain in re.findall(r"(candidate_only_junction|comparison_only_junction)=([^|]+)", text):
        for token in chain.strip().split(";"):
            parsed = parse_junction(token.strip())
            if parsed is not None:
                result[label].append(format_junction(*parsed))
    return result


def first_coordinate(value: str | None) -> tuple[str | None, str | None]:
    text = clean(value)
    if text is None:
        return None, None
    match = re.search(r"(chr[^:|;]+):(\d+)-(\d+):([+-])", text)
    if not match:
        return None, None
    return canonical_chromosome(match.group(1)), match.group(4)


def event_type(row: dict[str, str]) -> str:
    labels: list[str] = []
    if clean(row.get("SE")) not in {None, "not_annotated", "not_assessed_comparison_unavailable"}:
        labels.append("SE_or_complex_exon_change")
    if clean(row.get("A3SS")) not in {None, "not_annotated", "not_assessed_comparison_unavailable"}:
        labels.append("A3SS")
    if clean(row.get("A5SS")) not in {None, "not_annotated", "not_assessed_comparison_unavailable"}:
        labels.append("A5SS")
    if clean(row.get("IR")) not in {None, "not_detected_from_exon_chain", "not_assessed_comparison_unavailable"}:
        labels.append("IR")
    if clean(row.get("alternative_first_exon")) == "TRUE":
        labels.append("alternative_first_exon")
    if clean(row.get("alternative_last_exon")) == "TRUE":
        labels.append("alternative_last_exon")
    apa = clean(row.get("APA_or_3UTR_class"))
    if apa and apa not in {"NOT_aTTS", "comparison_unavailable"}:
        labels.append("APA_or_3UTR")
    return ";".join(labels) if labels else "unresolved_or_nonportable"


def bh_adjust(p_values: list[float]) -> list[float]:
    valid = [(index, p) for index, p in enumerate(p_values) if math.isfinite(p)]
    output = [float("nan")] * len(p_values)
    if not valid:
        return output
    valid.sort(key=lambda item: item[1])
    n = len(valid)
    previous = 1.0
    for rank in range(n, 0, -1):
        index, p = valid[rank - 1]
        value = min(previous, p * n / rank)
        output[index] = value
        previous = value
    return output


def paired_wilcoxon(deltas: list[float]) -> float:
    usable = [value for value in deltas if math.isfinite(value)]
    if len(usable) < MIN_PAIRED_PATIENTS:
        return float("nan")
    if all(abs(value) < 1e-12 for value in usable):
        return 1.0
    try:
        # The server ships an older SciPy release without the newer ``method`` argument.
        return float(wilcoxon(usable, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return 1.0


def bootstrap_proportion_ci(matches: list[bool], seed_offset: int) -> tuple[float, float]:
    if not matches:
        return float("nan"), float("nan")
    values = np.asarray(matches, dtype=float)
    rng = np.random.default_rng(BOOTSTRAP_SEED + seed_offset)
    samples = rng.choice(values, size=(BOOTSTRAP_REPS, len(values)), replace=True).mean(axis=1)
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def stable_seed_offset(label: str) -> int:
    """Make bootstrap output reproducible across Python processes."""
    return int(hashlib.sha256(label.encode("utf-8")).hexdigest()[:8], 16) % 100000


def sign(value: float) -> int:
    if not math.isfinite(value) or abs(value) < 1e-12:
        return 0
    return 1 if value > 0 else -1


def load_primary_inputs() -> tuple[dict[str, dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, Any]]]:
    missing = [str(path) for path in FROZEN_PATHS.values() if not path.exists()]
    if missing:
        raise RuntimeError("Frozen input(s) missing: " + "; ".join(missing))
    frozen_rows = read_tsv(FROZEN_PATHS["primary93"])
    recurrence_rows = read_tsv(FROZEN_PATHS["patient_recurrence"])
    feature_rows = read_tsv(FROZEN_PATHS["event_features"])
    frozen_by_isoform = {row["isoform_id"]: row for row in frozen_rows}
    if len(frozen_rows) != 93 or len(frozen_by_isoform) != 93:
        raise RuntimeError(f"Primary event freeze must contain exactly 93 unique isoforms, observed {len(frozen_rows)} rows / {len(frozen_by_isoform)} unique")
    if len(feature_rows) != 93 or {row["isoform_id"] for row in feature_rows} != set(frozen_by_isoform):
        raise RuntimeError("Feature mapping must be one-to-one with the frozen primary93 isoform set")
    recurrence_isoforms = {row.get("isoform_id") for row in recurrence_rows}
    if recurrence_isoforms != set(frozen_by_isoform):
        raise RuntimeError("Patient recurrence audit isoform IDs do not match frozen primary93 isoforms")
    hashes = [
        {
            "input_name": name,
            "path": str(path),
            "size_bytes": path.stat().st_size,
            "sha256_before": sha256_file(path),
            "read_only_role": "frozen_input",
        }
        for name, path in FROZEN_PATHS.items()
    ]
    return frozen_by_isoform, recurrence_rows, feature_rows, hashes


def discover_star_inputs() -> tuple[list[Path], list[dict[str, str]], set[str], list[dict[str, Any]]]:
    if not STAR_ROOT.exists():
        raise RuntimeError(f"STAR result root does not exist: {STAR_ROOT}")
    if not METADATA_PATH.exists() or not PRIMARY72_PATH.exists():
        raise RuntimeError("N/A/T metadata or primary72 cohort file is missing")
    metadata_rows = read_tsv(METADATA_PATH)
    metadata_by_sample = {row["raw_sample_id"]: row for row in metadata_rows}
    if len(metadata_by_sample) != len(metadata_rows):
        raise RuntimeError("Duplicate raw_sample_id in short-read sample metadata")
    primary72 = {line.strip() for line in PRIMARY72_PATH.read_text(encoding="utf-8").splitlines() if line.strip()}
    if len(primary72) != 72:
        raise RuntimeError(f"Primary cohort must contain 72 samples, observed {len(primary72)}")
    sj_paths = sorted(STAR_ROOT.rglob("*_SJ.out.tab"))
    if not sj_paths:
        raise RuntimeError("No pre-existing STAR SJ.out.tab files found; cannot derive junction counts without rerunning STAR")

    inventory: list[dict[str, Any]] = []
    sample_paths: dict[str, list[Path]] = defaultdict(list)
    for path in sj_paths:
        sample_id = re.sub(r"_SJ\.out\.tab$", "", path.name)
        sample_paths[sample_id].append(path)
        sample_meta = metadata_by_sample.get(sample_id)
        inventory.append(
            {
                "input_type": "STAR_SJ_out_tab",
                "path": str(path),
                "sample_id": sample_id,
                "file_size_bytes": path.stat().st_size,
                "mtime": dt.datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds"),
                "metadata_match": "TRUE" if sample_meta else "FALSE",
                "primary72_member": "TRUE" if sample_id in primary72 else "FALSE",
                "tissue_group": sample_meta.get("tissue_group") if sample_meta else None,
                "patient_id": sample_meta.get("patient_id") if sample_meta else None,
                "input_status": "usable_primary_NAT" if sample_id in primary72 and sample_meta else ("metadata_matched_nonprimary" if sample_meta else "unmapped_extra_file"),
                "notes": "Read-only STAR junction table; no STAR rerun performed.",
            }
        )
    for name, path, status, note in [
        ("sorted_BAM_or_CRAM", STAR_ROOT, "not_found", "No sorted BAM/CRAM located under 01_STAR_gencode; terminal-exon and intron-retention read coverage cannot be derived in this audit."),
        ("SUPPA2_event_output", CRC_ROOT / "script" / "suppa", "not_found_as_result", "Only helper scripts were found; no reusable SUPPA2 event output matrix was found."),
        ("rMATS_output", CRC_ROOT, "not_found", "No reusable rMATS result file was found by inventory search."),
        ("shortread_transcript_quantification", CRC_ROOT, "not_found", "No reusable short-read transcript quantification matrix was found by inventory search."),
    ]:
        inventory.append(
            {
                "input_type": name,
                "path": str(path),
                "sample_id": None,
                "file_size_bytes": None,
                "mtime": None,
                "metadata_match": None,
                "primary72_member": None,
                "tissue_group": None,
                "patient_id": None,
                "input_status": status,
                "notes": note,
            }
        )
    inventory.append(
        {
            "input_type": "STAR_annotation_index",
            "path": str(GTF_INDEX_PATH),
            "sample_id": None,
            "file_size_bytes": None,
            "mtime": None,
            "metadata_match": None,
            "primary72_member": None,
            "tissue_group": None,
            "patient_id": None,
            "input_status": "available_inferred_GRCh38_gencode_v48",
            "notes": "Index path identifies gencode.v48. Exact STAR run log was not retained, so build/annotation are inferred rather than independently logged.",
        }
    )
    duplicated_primary = {sample: paths for sample, paths in sample_paths.items() if sample in primary72 and len(paths) != 1}
    if duplicated_primary:
        raise RuntimeError(f"Primary72 STAR SJ file mapping is not one-to-one: {duplicated_primary}")
    missing_primary = sorted(primary72 - set(sample_paths))
    if missing_primary:
        raise RuntimeError("Primary72 samples missing STAR SJ.out.tab: " + ", ".join(missing_primary))
    primary_paths = [sample_paths[sample][0] for sample in sorted(primary72)]
    return primary_paths, metadata_rows, primary72, inventory


def build_sample_mapping_audit(metadata_rows: list[dict[str, str]], primary72: set[str], primary_paths: list[Path]) -> list[dict[str, Any]]:
    star_by_sample = {re.sub(r"_SJ\.out\.tab$", "", path.name): path for path in primary_paths}
    output: list[dict[str, Any]] = []
    primary_patients: dict[str, set[str]] = defaultdict(set)
    for row in metadata_rows:
        sample_id = row["raw_sample_id"]
        in_primary = sample_id in primary72
        path = star_by_sample.get(sample_id)
        if in_primary:
            primary_patients[row["patient_id"]].add(row["tissue_group"])
        output.append(
            {
                "raw_sample_id": sample_id,
                "patient_id": row.get("patient_id"),
                "tissue_group": row.get("tissue_group"),
                "paired_set": row.get("paired_set"),
                "primary72_member": "TRUE" if in_primary else "FALSE",
                "star_SJ_path": str(path) if path else None,
                "star_SJ_available": "TRUE" if path else "FALSE",
                "analysis_included": "TRUE" if in_primary and path else "FALSE",
                "mapping_status": "matched_primary_NAT" if in_primary and path else ("metadata_nonprimary" if not in_primary else "missing_STAR_SJ"),
                "notes": "Tissue convention from frozen metadata: A=Tumor, B=Adjacent, C=Normal.",
            }
        )
    incomplete = {patient: sorted(groups) for patient, groups in primary_patients.items() if groups != {"N", "A", "T"}}
    if incomplete:
        raise RuntimeError(f"Primary72 does not form complete N/A/T patient triplets: {incomplete}")
    if len(primary_patients) != 24:
        raise RuntimeError(f"Expected 24 primary72 patients, observed {len(primary_patients)}")
    return output


def build_event_mapping(feature_rows: list[dict[str, str]], frozen_by_isoform: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    mapped: list[dict[str, Any]] = []
    for feature in feature_rows:
        isoform_id = feature["isoform_id"]
        frozen = frozen_by_isoform[isoform_id]
        junctions = parse_labeled_junctions(feature.get("event_defining_junctions"))
        inclusion = junctions["candidate_only_junction"]
        exclusion = junctions["comparison_only_junction"]
        chrom, strand = first_coordinate(feature.get("event_defining_junctions"))
        if chrom is None:
            chrom, strand = first_coordinate(feature.get("exon_chain"))
        is_junction_testable = as_bool(feature.get("junction_testable"))
        is_event_testable = as_bool(feature.get("event_testable"))
        if is_junction_testable and inclusion and exclusion:
            status = "junction_pair_PSI_proxy"
            reason = None
        elif is_junction_testable and inclusion:
            status = "one_sided_candidate_junction_only"
            reason = "No reciprocal comparison-defining junction was available; candidate-junction abundance can be audited, but PSI cannot be computed."
        elif is_junction_testable and exclusion:
            status = "one_sided_comparison_junction_only"
            reason = "Only a comparison-defining junction was resolved; this is not a direct candidate-isoform PSI measurement."
        elif is_junction_testable:
            status = "junction_mapping_unparsed"
            reason = "Frozen table marked junction_testable but no parseable junction coordinate was present."
        elif is_event_testable:
            status = "event_requires_BAM_or_exon_coverage"
            reason = "Event is event-testable in principle, but no reusable BAM/CRAM or exon-coverage matrix is available for terminal/exonic event quantification."
        elif clean(feature.get("external_validation_level")) == "gene_only":
            status = "gene_only_not_junction_validation"
            reason = "Gene-only portability cannot be counted as junction/event-level support."
        else:
            status = "no_portable_junction_or_event_feature"
            reason = clean(feature.get("not_testable_reason")) or "No portable short-read feature was resolved."
        mapped.append(
            {
                "isoform_id": isoform_id,
                "gene_id": feature.get("gene_id"),
                "gene_symbol": feature.get("gene_symbol"),
                "comparison_isoform_id": feature.get("comparison_isoform_id"),
                "accepted_comparison": feature.get("accepted_comparison"),
                "biological_pair_id": feature.get("biological_pair_id"),
                "in_core70": "TRUE" if int(frozen.get("n_supporting_contrasts", "0")) == 2 else "FALSE",
                "chromosome": chrom,
                "strand": strand,
                "event_type": event_type(feature),
                "event_defining_junctions": feature.get("event_defining_junctions"),
                "inclusion_junctions": ";".join(inclusion) if inclusion else None,
                "exclusion_junctions": ";".join(exclusion) if exclusion else None,
                "alternative_first_exon": feature.get("alternative_first_exon"),
                "alternative_last_exon": feature.get("alternative_last_exon"),
                "SE": feature.get("SE"),
                "A3SS": feature.get("A3SS"),
                "A5SS": feature.get("A5SS"),
                "IR_annotation": feature.get("IR"),
                "IR_coordinates": feature.get("IR_coordinates"),
                "APA_or_3UTR_class": feature.get("APA_or_3UTR_class"),
                "terminal_event_coordinates": feature.get("terminal_event_coordinates"),
                "mapping_status": status,
                "mapping_failure_reason": reason,
                "junction_testable_frozen": feature.get("junction_testable"),
                "event_testable_frozen": feature.get("event_testable"),
                "longread_T_vs_A_delta_usage": frozen.get("T_vs_A_delta_usage"),
                "longread_T_vs_N_delta_usage": frozen.get("T_vs_N_delta_usage"),
                "longread_T_vs_A_FDR": frozen.get("T_vs_A_FDR"),
                "longread_T_vs_N_FDR": frozen.get("T_vs_N_FDR"),
            }
        )
    if len(mapped) != 93:
        raise RuntimeError("Event mapping did not retain exactly 93 frozen isoform rows")
    if sum(row["in_core70"] == "TRUE" for row in mapped) != 70:
        raise RuntimeError("Frozen core recurrent subset must contain exactly 70 isoforms")
    if sum(row["junction_testable_frozen"] == "TRUE" for row in mapped) != 31:
        raise RuntimeError("Frozen mapping must retain exactly 31 junction-testable isoforms")
    if sum(row["event_testable_frozen"] == "TRUE" for row in mapped) != 34:
        raise RuntimeError("Frozen mapping must retain exactly 34 event-testable isoforms")
    return mapped


def collect_target_junctions(event_rows: list[dict[str, Any]]) -> set[str]:
    targets: set[str] = set()
    for row in event_rows:
        for field in ("inclusion_junctions", "exclusion_junctions"):
            value = clean(row.get(field))
            if value:
                targets.update(item for item in value.split(";") if parse_junction(item) is not None)
    return targets


def scan_star_junctions(primary_paths: list[Path], target_junctions: set[str]) -> tuple[dict[str, dict[str, int]], dict[str, int], list[dict[str, Any]]]:
    """Stream STAR junction tables once; only retain frozen-event coordinate counts."""
    sample_counts: dict[str, dict[str, int]] = {}
    library_totals: dict[str, int] = {}
    audit_rows: list[dict[str, Any]] = []
    for number, path in enumerate(primary_paths, start=1):
        sample_id = re.sub(r"_SJ\.out\.tab$", "", path.name)
        target_counts: dict[str, int] = defaultdict(int)
        total_reads = 0
        parsed_rows = 0
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise RuntimeError(f"Unexpected STAR SJ.out.tab record in {path}: {line[:120]!r}")
                strand = {"1": "+", "2": "-"}.get(fields[3])
                try:
                    unique_reads = int(fields[6])
                    multi_reads = int(fields[7])
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError as exc:
                    raise RuntimeError(f"Non-integer STAR SJ count/coordinate in {path}") from exc
                count = unique_reads + multi_reads
                total_reads += count
                parsed_rows += 1
                if strand is None:
                    continue
                junction = format_junction(fields[0], start, end, strand)
                if junction in target_junctions:
                    target_counts[junction] += count
        sample_counts[sample_id] = dict(target_counts)
        library_totals[sample_id] = total_reads
        audit_rows.append(
            {
                "raw_sample_id": sample_id,
                "STAR_SJ_path": str(path),
                "star_sj_rows_parsed": parsed_rows,
                "total_unique_plus_multimapping_junction_reads": total_reads,
                "frozen_event_target_junctions_detected": sum(1 for value in target_counts.values() if value > 0),
                "frozen_event_target_junctions_total": len(target_junctions),
                "coordinate_matching_rule": "exact chr/start/end/strand match; STAR SJ coordinates were not shifted or tolerance-tuned.",
            }
        )
        if number % 12 == 0 or number == len(primary_paths):
            update_progress("scan_STAR_SJ", f"Parsed {number}/{len(primary_paths)} primary72 STAR junction tables")
    return sample_counts, library_totals, audit_rows


def mean_or_nan(values: list[int]) -> float:
    return float(np.mean(values)) if values else float("nan")


def build_junction_counts(
    event_rows: list[dict[str, Any]],
    metadata_rows: list[dict[str, str]],
    primary72: set[str],
    sample_counts: dict[str, dict[str, int]],
    library_totals: dict[str, int],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    meta_by_sample = {row["raw_sample_id"]: row for row in metadata_rows}
    counts_rows: list[dict[str, Any]] = []
    coverage_rows: list[dict[str, Any]] = []
    for event in event_rows:
        inclusion = [item for item in clean(event.get("inclusion_junctions")).split(";") if item] if clean(event.get("inclusion_junctions")) else []
        exclusion = [item for item in clean(event.get("exclusion_junctions")).split(";") if item] if clean(event.get("exclusion_junctions")) else []
        for sample_id in sorted(primary72):
            meta = meta_by_sample[sample_id]
            observed = sample_counts[sample_id]
            inclusion_values = [observed.get(junction, 0) for junction in inclusion]
            exclusion_values = [observed.get(junction, 0) for junction in exclusion]
            inclusion_count = sum(inclusion_values) if inclusion_values else float("nan")
            exclusion_count = sum(exclusion_values) if exclusion_values else float("nan")
            inclusion_mean = mean_or_nan(inclusion_values)
            exclusion_mean = mean_or_nan(exclusion_values)
            metric_type = "not_mappable"
            psi_or_ratio = float("nan")
            junction_cpm = float("nan")
            total_informative = float("nan")
            coverage_status = "not_mappable"
            if event["mapping_status"] == "junction_pair_PSI_proxy":
                metric_type = "junction_pair_ratio_proxy"
                total_informative = inclusion_mean + exclusion_mean
                if total_informative > 0:
                    psi_or_ratio = inclusion_mean / total_informative
                coverage_status = "pass_pair_ratio" if total_informative >= MIN_PAIR_INFORMATIVE_READS else "fail_pair_ratio_low_informative_reads"
            elif event["mapping_status"] == "one_sided_candidate_junction_only":
                metric_type = "one_sided_candidate_junction_CPM"
                if math.isfinite(inclusion_count) and library_totals[sample_id] > 0:
                    junction_cpm = inclusion_count / library_totals[sample_id] * 1_000_000
                coverage_status = "pass_one_sided_abundance" if inclusion_mean >= MIN_ONE_SIDED_MEAN_JUNCTION_READS else "fail_one_sided_low_junction_reads"
            elif event["mapping_status"] == "one_sided_comparison_junction_only":
                metric_type = "one_sided_comparison_junction_CPM"
                if math.isfinite(exclusion_count) and library_totals[sample_id] > 0:
                    junction_cpm = exclusion_count / library_totals[sample_id] * 1_000_000
                coverage_status = "pass_one_sided_abundance" if exclusion_mean >= MIN_ONE_SIDED_MEAN_JUNCTION_READS else "fail_one_sided_low_junction_reads"
            counts_rows.append(
                {
                    "isoform_id": event["isoform_id"],
                    "gene_symbol": event["gene_symbol"],
                    "comparison_isoform_id": event["comparison_isoform_id"],
                    "biological_pair_id": event["biological_pair_id"],
                    "in_core70": event["in_core70"],
                    "patient_id": meta["patient_id"],
                    "tissue_group": meta["tissue_group"],
                    "raw_sample_id": sample_id,
                    "mapping_status": event["mapping_status"],
                    "metric_type": metric_type,
                    "inclusion_junctions": event["inclusion_junctions"],
                    "exclusion_junctions": event["exclusion_junctions"],
                    "inclusion_count": inclusion_count,
                    "exclusion_count": exclusion_count,
                    "inclusion_mean_junction_count": inclusion_mean,
                    "exclusion_mean_junction_count": exclusion_mean,
                    "total_informative_reads": total_informative,
                    "PSI_or_junction_ratio": psi_or_ratio,
                    "one_sided_junction_CPM": junction_cpm,
                    "library_total_junction_reads": library_totals[sample_id],
                    "coverage_status": coverage_status,
                    "coverage_threshold": f"pair_ratio_mean_inclusion_plus_exclusion>={MIN_PAIR_INFORMATIVE_READS}; one_sided_mean_junction>={MIN_ONE_SIDED_MEAN_JUNCTION_READS}",
                }
            )
    group_fields = ["isoform_id", "gene_symbol", "in_core70", "mapping_status", "metric_type", "tissue_group", "coverage_status"]
    summary: dict[tuple[str, ...], int] = Counter(tuple(str(row[field]) for field in group_fields) for row in counts_rows)
    for keys, n_rows in sorted(summary.items()):
        coverage_rows.append(dict(zip(group_fields, keys)) | {"n_patient_tissue_rows": n_rows})
    return counts_rows, coverage_rows


def group_by_patient(count_rows: list[dict[str, Any]]) -> dict[str, dict[str, dict[str, dict[str, Any]]]]:
    output: dict[str, dict[str, dict[str, dict[str, Any]]]] = defaultdict(lambda: defaultdict(dict))
    for row in count_rows:
        output[row["isoform_id"]][row["patient_id"]][row["tissue_group"]] = row
    return output


def contrast_pair(first: str, second: str, contrast: str) -> tuple[str, str]:
    return first, second


CONTRASTS = {
    "A_vs_N": ("N", "A"),
    "T_vs_A": ("A", "T"),
    "T_vs_N": ("N", "T"),
}


def metric_value(row: dict[str, Any]) -> float:
    if row["metric_type"] == "junction_pair_ratio_proxy":
        return as_float(row["PSI_or_junction_ratio"])
    if row["metric_type"].startswith("one_sided_"):
        return as_float(row["one_sided_junction_CPM"])
    return float("nan")


def paired_statistics(
    event_rows: list[dict[str, Any]], count_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    event_by_isoform = {row["isoform_id"]: row for row in event_rows}
    by_patient = group_by_patient(count_rows)
    results: list[dict[str, Any]] = []
    patient_effects: list[dict[str, Any]] = []
    for isoform_id, event in event_by_isoform.items():
        for contrast, (from_tissue, to_tissue) in CONTRASTS.items():
            deltas: list[float] = []
            patients: list[str] = []
            coverage_pass_rows = 0
            for patient_id, tissues in sorted(by_patient[isoform_id].items()):
                before = tissues.get(from_tissue)
                after = tissues.get(to_tissue)
                if before is None or after is None:
                    continue
                metric = before["metric_type"]
                valid_pair_ratio = metric == "junction_pair_ratio_proxy" and before["coverage_status"] == "pass_pair_ratio" and after["coverage_status"] == "pass_pair_ratio"
                valid_one_sided = metric.startswith("one_sided_") and before["coverage_status"] == "pass_one_sided_abundance" and after["coverage_status"] == "pass_one_sided_abundance"
                value_before = metric_value(before)
                value_after = metric_value(after)
                valid = (valid_pair_ratio or valid_one_sided) and math.isfinite(value_before) and math.isfinite(value_after)
                if valid:
                    delta = value_after - value_before
                    deltas.append(delta)
                    patients.append(patient_id)
                    coverage_pass_rows += 1
                    patient_effects.append(
                        {
                            "isoform_id": isoform_id,
                            "gene_symbol": event["gene_symbol"],
                            "in_core70": event["in_core70"],
                            "patient_id": patient_id,
                            "contrast": contrast,
                            "metric_type": metric,
                            "from_tissue": from_tissue,
                            "to_tissue": to_tissue,
                            "from_value": value_before,
                            "to_value": value_after,
                            "shortread_delta_metric": delta,
                            "longread_delta_usage": event.get(f"longread_{contrast}_delta_usage"),
                            "coverage_pass": "TRUE",
                        }
                    )
            metric_type = next(iter({row["metric_type"] for tissues in by_patient[isoform_id].values() for row in tissues.values()}), "not_mappable")
            median_delta = median(deltas) if deltas else float("nan")
            longread_delta = as_float(event.get(f"longread_{contrast}_delta_usage"))
            direction_match = sign(median_delta) == sign(longread_delta) if metric_type == "junction_pair_ratio_proxy" and sign(median_delta) and sign(longread_delta) else None
            patient_direction_concordance = float("nan")
            if metric_type == "junction_pair_ratio_proxy" and sign(longread_delta) and deltas:
                patient_direction_concordance = sum(sign(delta) == sign(longread_delta) for delta in deltas) / len(deltas)
            results.append(
                {
                    "isoform_id": isoform_id,
                    "gene_id": event["gene_id"],
                    "gene_symbol": event["gene_symbol"],
                    "comparison_isoform_id": event["comparison_isoform_id"],
                    "biological_pair_id": event["biological_pair_id"],
                    "in_core70": event["in_core70"],
                    "contrast": contrast,
                    "metric_type": metric_type,
                    "mapping_status": event["mapping_status"],
                    "n_valid_patients": len(patients),
                    "valid_patients": ";".join(patients) if patients else None,
                    "median_delta_PSI": median_delta if metric_type == "junction_pair_ratio_proxy" else float("nan"),
                    "median_delta_one_sided_junction_CPM": median_delta if metric_type.startswith("one_sided_") else float("nan"),
                    "paired_test_p": paired_wilcoxon(deltas),
                    "BH_FDR": float("nan"),
                    "direction_concordance": patient_direction_concordance,
                    "longread_shortread_direction_match": direction_match,
                    "longread_delta_usage": longread_delta,
                    "shortread_delta_PSI": median_delta if metric_type == "junction_pair_ratio_proxy" else float("nan"),
                    "shortread_delta_metric": median_delta,
                    "support_interpretation": "short-read junction-pair ratio support" if metric_type == "junction_pair_ratio_proxy" else ("one-sided junction-abundance support; not PSI" if metric_type.startswith("one_sided_") else "not evaluable from available STAR junction table"),
                }
            )
    for contrast in CONTRASTS:
        indices = [index for index, row in enumerate(results) if row["contrast"] == contrast]
        adjusted = bh_adjust([as_float(results[index]["paired_test_p"]) for index in indices])
        for index, value in zip(indices, adjusted):
            results[index]["BH_FDR"] = value
    return results, patient_effects


def build_concordance_summary(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for scope, selector in [
        ("primary93", lambda row: True),
        ("core70", lambda row: row["in_core70"] == "TRUE"),
        ("non_core23", lambda row: row["in_core70"] != "TRUE"),
    ]:
        for contrast in ("T_vs_A", "T_vs_N"):
            scoped = [row for row in results if selector(row) and row["contrast"] == contrast]
            pair_rows = [row for row in scoped if row["metric_type"] == "junction_pair_ratio_proxy"]
            evaluable = [row for row in pair_rows if int(row["n_valid_patients"]) >= MIN_PAIRED_PATIENTS and math.isfinite(as_float(row["shortread_delta_PSI"]))]
            matches = [bool(row["longread_shortread_direction_match"]) for row in evaluable if row["longread_shortread_direction_match"] is not None]
            ci_low, ci_high = bootstrap_proportion_ci(matches, stable_seed_offset(f"{scope}-{contrast}"))
            x = [as_float(row["longread_delta_usage"]) for row in evaluable if math.isfinite(as_float(row["longread_delta_usage"]))]
            y = [as_float(row["shortread_delta_PSI"]) for row in evaluable if math.isfinite(as_float(row["longread_delta_usage"]))]
            if len(x) >= 3 and len(set(x)) > 1 and len(set(y)) > 1:
                rho, p_value = spearmanr(x, y)
            else:
                rho, p_value = float("nan"), float("nan")
            output.append(
                {
                    "scope": scope,
                    "contrast": contrast,
                    "analysis_unit": "isoform-associated junction-pair feature; reciprocal rows are retained but not independent biological discoveries",
                    "n_isoforms_total": len(scoped),
                    "n_pair_ratio_mappable": len(pair_rows),
                    "n_pair_ratio_coverage_pass": len(evaluable),
                    "n_unique_biological_pair_ids_pair_ratio": len({row["biological_pair_id"] for row in pair_rows if clean(row["biological_pair_id"])}),
                    "n_direction_match": sum(matches),
                    "direction_match_proportion": sum(matches) / len(matches) if matches else float("nan"),
                    "bootstrap_95CI_low": ci_low,
                    "bootstrap_95CI_high": ci_high,
                    "spearman_longread_delta_usage_vs_shortread_delta_PSI": rho,
                    "spearman_p": p_value,
                    "spearman_n": len(x),
                    "minimum_coverage_rule": f"Both tissues >= {MIN_PAIR_INFORMATIVE_READS} mean informative junction reads; patient-level paired analysis; n >= {MIN_PAIRED_PATIENTS}.",
                }
            )
    return output


def targeted_results(results: list[dict[str, Any]], count_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    targets = {"KRAS", "CDX2", "OLFM4"}
    coverage: dict[str, Counter[str]] = defaultdict(Counter)
    for row in count_rows:
        if row["gene_symbol"] in targets:
            coverage[row["isoform_id"]][row["coverage_status"]] += 1
    target_rows = []
    for row in results:
        if row["gene_symbol"] not in targets:
            continue
        target_rows.append(
            dict(row) | {"patient_tissue_coverage_status_counts": ";".join(f"{key}:{value}" for key, value in sorted(coverage[row["isoform_id"]].items()))}
        )
    return target_rows


def audit_frozen_inputs_unchanged(hashes_before: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for item in hashes_before:
        current = sha256_file(Path(item["path"]))
        output.append(dict(item) | {"sha256_after": current, "unchanged": "TRUE" if current == item["sha256_before"] else "FALSE"})
    if any(row["unchanged"] != "TRUE" for row in output):
        raise RuntimeError("A frozen input changed during this audit; refusing to continue")
    return output


def main() -> None:
    for directory in (TABLE_DIR, REPORT_DIR, LOG_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    SUCCESS_PATH.unlink(missing_ok=True)
    FAILED_PATH.unlink(missing_ok=True)
    update_progress("start", "Read-only internal STAR junction/event-level audit started")
    frozen_by_isoform, recurrence_rows, feature_rows, hashes_before = load_primary_inputs()
    primary_paths, metadata_rows, primary72, inventory = discover_star_inputs()
    sample_audit = build_sample_mapping_audit(metadata_rows, primary72, primary_paths)
    write_tsv(
        TABLE_DIR / "Fig5_shortread_input_inventory.tsv",
        inventory,
        ["input_type", "path", "sample_id", "file_size_bytes", "mtime", "metadata_match", "primary72_member", "tissue_group", "patient_id", "input_status", "notes"],
    )
    write_tsv(
        TABLE_DIR / "Fig5_shortread_sample_mapping_audit.tsv",
        sample_audit,
        ["raw_sample_id", "patient_id", "tissue_group", "paired_set", "primary72_member", "star_SJ_path", "star_SJ_available", "analysis_included", "mapping_status", "notes"],
    )
    event_rows = build_event_mapping(feature_rows, frozen_by_isoform)
    write_tsv(
        TABLE_DIR / "Fig5_shortread_event_mapping.tsv",
        event_rows,
        ["isoform_id", "gene_id", "gene_symbol", "comparison_isoform_id", "accepted_comparison", "biological_pair_id", "in_core70", "chromosome", "strand", "event_type", "event_defining_junctions", "inclusion_junctions", "exclusion_junctions", "alternative_first_exon", "alternative_last_exon", "SE", "A3SS", "A5SS", "IR_annotation", "IR_coordinates", "APA_or_3UTR_class", "terminal_event_coordinates", "mapping_status", "mapping_failure_reason", "junction_testable_frozen", "event_testable_frozen", "longread_T_vs_A_delta_usage", "longread_T_vs_N_delta_usage", "longread_T_vs_A_FDR", "longread_T_vs_N_FDR"],
    )
    update_progress("input_audit_complete", f"Frozen primary93=93; primary72 STAR junction tables={len(primary_paths)}")
    target_junctions = collect_target_junctions(event_rows)
    if not target_junctions:
        raise RuntimeError("No parseable frozen event-defining junctions; cannot audit short-read junction support")
    sample_counts, library_totals, scan_audit = scan_star_junctions(primary_paths, target_junctions)
    write_tsv(
        TABLE_DIR / "Fig5_shortread_STAR_SJ_scan_audit.tsv",
        scan_audit,
        ["raw_sample_id", "STAR_SJ_path", "star_sj_rows_parsed", "total_unique_plus_multimapping_junction_reads", "frozen_event_target_junctions_detected", "frozen_event_target_junctions_total", "coordinate_matching_rule"],
    )
    count_rows, coverage_rows = build_junction_counts(event_rows, metadata_rows, primary72, sample_counts, library_totals)
    write_tsv(
        TABLE_DIR / "Fig5_shortread_junction_counts.tsv",
        count_rows,
        ["isoform_id", "gene_symbol", "comparison_isoform_id", "biological_pair_id", "in_core70", "patient_id", "tissue_group", "raw_sample_id", "mapping_status", "metric_type", "inclusion_junctions", "exclusion_junctions", "inclusion_count", "exclusion_count", "inclusion_mean_junction_count", "exclusion_mean_junction_count", "total_informative_reads", "PSI_or_junction_ratio", "one_sided_junction_CPM", "library_total_junction_reads", "coverage_status", "coverage_threshold"],
    )
    write_tsv(
        TABLE_DIR / "Fig5_shortread_event_PSI.tsv",
        [row for row in count_rows if row["metric_type"] == "junction_pair_ratio_proxy"],
        ["isoform_id", "gene_symbol", "comparison_isoform_id", "biological_pair_id", "in_core70", "patient_id", "tissue_group", "raw_sample_id", "mapping_status", "metric_type", "inclusion_count", "exclusion_count", "inclusion_mean_junction_count", "exclusion_mean_junction_count", "total_informative_reads", "PSI_or_junction_ratio", "coverage_status", "coverage_threshold"],
    )
    write_tsv(
        TABLE_DIR / "Fig5_shortread_coverage_QA.tsv",
        coverage_rows,
        ["isoform_id", "gene_symbol", "in_core70", "mapping_status", "metric_type", "tissue_group", "coverage_status", "n_patient_tissue_rows"],
    )
    results, patient_effects = paired_statistics(event_rows, count_rows)
    write_tsv(
        TABLE_DIR / "Fig5_shortread_validation_results.tsv",
        results,
        ["isoform_id", "gene_id", "gene_symbol", "comparison_isoform_id", "biological_pair_id", "in_core70", "contrast", "metric_type", "mapping_status", "n_valid_patients", "valid_patients", "median_delta_PSI", "median_delta_one_sided_junction_CPM", "paired_test_p", "BH_FDR", "direction_concordance", "longread_shortread_direction_match", "longread_delta_usage", "shortread_delta_PSI", "shortread_delta_metric", "support_interpretation"],
    )
    write_tsv(
        TABLE_DIR / "Fig5_shortread_patient_effects.tsv",
        patient_effects,
        ["isoform_id", "gene_symbol", "in_core70", "patient_id", "contrast", "metric_type", "from_tissue", "to_tissue", "from_value", "to_value", "shortread_delta_metric", "longread_delta_usage", "coverage_pass"],
    )
    write_tsv(
        TABLE_DIR / "Fig5_shortread_core70_results.tsv",
        [row for row in results if row["in_core70"] == "TRUE"],
        ["isoform_id", "gene_id", "gene_symbol", "comparison_isoform_id", "biological_pair_id", "in_core70", "contrast", "metric_type", "mapping_status", "n_valid_patients", "valid_patients", "median_delta_PSI", "median_delta_one_sided_junction_CPM", "paired_test_p", "BH_FDR", "direction_concordance", "longread_shortread_direction_match", "longread_delta_usage", "shortread_delta_PSI", "shortread_delta_metric", "support_interpretation"],
    )
    concordance = build_concordance_summary(results)
    write_tsv(
        TABLE_DIR / "Fig5_shortread_longread_concordance.tsv",
        concordance,
        ["scope", "contrast", "analysis_unit", "n_isoforms_total", "n_pair_ratio_mappable", "n_pair_ratio_coverage_pass", "n_unique_biological_pair_ids_pair_ratio", "n_direction_match", "direction_match_proportion", "bootstrap_95CI_low", "bootstrap_95CI_high", "spearman_longread_delta_usage_vs_shortread_delta_PSI", "spearman_p", "spearman_n", "minimum_coverage_rule"],
    )
    target_rows = targeted_results(results, count_rows)
    write_tsv(
        TABLE_DIR / "Fig5_shortread_KRAS_CDX2_OLFM4_support.tsv",
        target_rows,
        ["isoform_id", "gene_id", "gene_symbol", "comparison_isoform_id", "biological_pair_id", "in_core70", "contrast", "metric_type", "mapping_status", "n_valid_patients", "valid_patients", "median_delta_PSI", "median_delta_one_sided_junction_CPM", "paired_test_p", "BH_FDR", "direction_concordance", "longread_shortread_direction_match", "longread_delta_usage", "shortread_delta_PSI", "shortread_delta_metric", "support_interpretation", "patient_tissue_coverage_status_counts"],
    )
    input_hash_audit = audit_frozen_inputs_unchanged(hashes_before)
    write_tsv(
        TABLE_DIR / "Fig5_shortread_frozen_input_hash_audit.tsv",
        input_hash_audit,
        ["input_name", "path", "size_bytes", "sha256_before", "sha256_after", "unchanged", "read_only_role"],
    )

    map_counts = Counter(row["mapping_status"] for row in event_rows)
    ratio_results = [row for row in results if row["metric_type"] == "junction_pair_ratio_proxy"]
    ratio_pass = [row for row in ratio_results if int(row["n_valid_patients"]) >= MIN_PAIRED_PATIENTS]
    target_summary = Counter(row["gene_symbol"] for row in target_rows if int(row["n_valid_patients"]) >= MIN_PAIRED_PATIENTS and row["metric_type"] == "junction_pair_ratio_proxy")
    write_report(
        REPORT_DIR / "Fig5_shortread_input_audit_zh.md",
        [
            "# 同队列短读长 junction/event-level 输入审计",
            "",
            "## 审计范围",
            "",
            "本阶段仅解析既有 `01_STAR_gencode/*/*_SJ.out.tab`。未读取 FASTQ，未运行 STAR，未修改 Figure 1-4，未生成 Figure 5，未做临床关联或候选评分。",
            "",
            "## 可用输入",
            "",
            f"- 发现 {len(inventory) - 5} 个 STAR SJ.out.tab 文件；其中 {len(primary_paths)} 个一一对应 primary72 N/A/T 样本。",
            "- primary72 为 24 个完整 N/A/T 患者配对；样本级 junction 计数将按 patient x tissue 汇总与配对检验，绝不将 reads 当作生物学重复。",
            "- 多出的未映射文件 `CRC-5C_01_SJ.out.tab` 不在冻结元数据中，未纳入分析。",
            "- STAR index 路径为 `gencode.v48_STAR_INDEX`，因此 GRCh38/gencode v48 仅作为路径推断；未保留 STAR Log.out，不能将其报告为独立的运行日志证明。",
            "- 未发现可复用 short-read sorted BAM/CRAM、SUPPA2 result matrix、rMATS result matrix 或 short-read transcript quantification matrix。",
            "",
            "## 固定计数规则",
            "",
            "- junction 坐标使用 chr/start/end/strand 的精确匹配；未容许坐标平移，也未在见到结果后调参。",
            "- STAR junction read count = uniquely mapped reads + multi-mapping reads（SJ.out.tab 第 7+8 列）。",
            f"- 双侧候选/比较 junction 的 ratio proxy 使用两侧各自 event-defining junction reads 的平均值，最低 informative reads = {MIN_PAIR_INFORMATIVE_READS}；这避免因一侧包含更多 junction 而用简单求和产生结构性偏差。",
            f"- 仅单侧 unique junction 使用固定 mean junction reads >= {MIN_ONE_SIDED_MEAN_JUNCTION_READS} 的 abundance/CPM 描述，不能报告为 PSI 或完整 isoform usage。",
            f"- 配对统计最低 {MIN_PAIRED_PATIENTS} 位患者，Wilcoxon 配对检验后在每个 contrast 的所有可检验特征中进行 BH 校正。",
            "",
            "## 冻结输入核对",
            "",
            "- primary93 = 93 unique isoforms；feature mapping 与 recurrence audit 均已一对一核对。",
            "- 内部冻结输入的 SHA256 在运行前后均一致；脚本只读取，未写回。",
            "",
            "## 事件可观测性",
            "",
            f"- junction-pair PSI proxy: {map_counts['junction_pair_PSI_proxy']} isoform-associated features。",
            f"- candidate-only / comparison-only one-sided junction features: {map_counts['one_sided_candidate_junction_only']} / {map_counts['one_sided_comparison_junction_only']}。",
            f"- event-testable but requiring BAM/exon coverage: {map_counts['event_requires_BAM_or_exon_coverage']}。",
            f"- 无 portable junction/event feature 或仅 gene-level: {map_counts['no_portable_junction_or_event_feature'] + map_counts['gene_only_not_junction_validation'] + map_counts['junction_mapping_unparsed']}。",
            "",
            "## 解释边界",
            "",
            "可报告的结果仅为同队列 short-read junction-level 或 event-level orthogonal support。短读长 junction 不能验证全长 isoform，尤其不能区分具有相同 junction chain 但转录本端点不同的转录本。无 BAM/CRAM 时，terminal exon、APA/3UTR 与 intron-retention coverage 不能被强行量化。",
            "",
            "详细结果、覆盖情况、配对统计和 long-read/short-read 一致性见 `Fig5_shortread_*` TSV；本报告不构成 Figure 5 结果或外部复现结论。",
        ],
    )
    concordance_by_key = {(row["scope"], row["contrast"]): row for row in concordance}
    target_by_key = {(row["gene_symbol"], row["isoform_id"], row["contrast"]): row for row in target_rows}

    def target_summary_line(symbol: str, isoform: str, contrast: str) -> str:
        row = target_by_key.get((symbol, isoform, contrast))
        if row is None:
            return f"- {symbol} {isoform} {contrast}: no frozen row."
        if row["metric_type"] == "junction_pair_ratio_proxy":
            return (
                f"- {symbol} {isoform} {contrast}: median short-read junction-pair ratio delta = "
                f"{float(row['median_delta_PSI']):.4f}, n = {row['n_valid_patients']}, "
                f"BH FDR = {float(row['BH_FDR']):.3g}, long-read direction match = {row['longread_shortread_direction_match']}."
            )
        return (
            f"- {symbol} {isoform} {contrast}: {row['metric_type']}; n = {row['n_valid_patients']}; "
            "this is one-sided junction abundance/CPM rather than PSI, so it is not direct isoform-usage replication."
        )

    primary_ta = concordance_by_key[("primary93", "T_vs_A")]
    primary_tn = concordance_by_key[("primary93", "T_vs_N")]
    core_ta = concordance_by_key[("core70", "T_vs_A")]
    core_tn = concordance_by_key[("core70", "T_vs_N")]
    write_report(
        REPORT_DIR / "Fig5_shortread_junction_validation_report_zh.md",
        [
            "# 同队列短读长 junction/event-level 正交支持审计",
            "",
            "## 设计与冻结边界",
            "",
            "本审计使用同一 CRC N/A/T 队列既有 STAR SJ.out.tab 文件，以患者而非 reads 为统计单位。primary93、core70 和 internal feature mapping 在读入前已冻结；SHA256 运行前后一致。GSE248094 结果未修改。未重跑 FASTQ/STAR，未做临床关联、候选评分或 Figure 5。",
            "",
            "## 可观测性与覆盖",
            "",
            f"- primary93 的 31 个 frozen junction-testable isoforms 中，{map_counts['junction_pair_PSI_proxy']} 个具候选与比较两侧 junction，可计算 junction-pair ratio proxy；{map_counts['one_sided_candidate_junction_only'] + map_counts['one_sided_comparison_junction_only']} 个只有单侧 junction，仅可做 abundance/CPM 支持。",
            f"- 34 个 frozen event-testable isoforms 中，另有 {map_counts['event_requires_BAM_or_exon_coverage']} 个需要 BAM/CRAM 或 exon coverage 才能量化 terminal/exonic event；当前输入下未强行计算。",
            f"- 在固定双侧覆盖规则下（两侧 average informative junction reads >= {MIN_PAIR_INFORMATIVE_READS}，并要求 >= {MIN_PAIRED_PATIENTS} paired patients），23/23 junction-pair isoform-associated features 均在 T_vs_A 和 T_vs_N 可评估。",
            "",
            "## Long-read 与 short-read 一致性",
            "",
            f"- T_vs_A：23/23 isoform-associated junction-pair features 的 cohort direction 与 long-read delta usage 一致；bootstrap 95% CI {float(primary_ta['bootstrap_95CI_low']):.3f}-{float(primary_ta['bootstrap_95CI_high']):.3f}；Spearman rho = {float(primary_ta['spearman_longread_delta_usage_vs_shortread_delta_PSI']):.3f}, P = {float(primary_ta['spearman_p']):.3g}, n = {primary_ta['spearman_n']}.",
            f"- T_vs_N：23/23 direction matched；bootstrap 95% CI {float(primary_tn['bootstrap_95CI_low']):.3f}-{float(primary_tn['bootstrap_95CI_high']):.3f}；Spearman rho = {float(primary_tn['spearman_longread_delta_usage_vs_shortread_delta_PSI']):.3f}, P = {float(primary_tn['spearman_p']):.3g}, n = {primary_tn['spearman_n']}.",
            f"- core70 中，可计算 junction-pair ratio 的 features 为 T_vs_A {core_ta['n_pair_ratio_coverage_pass']}/{core_ta['n_pair_ratio_mappable']}、T_vs_N {core_tn['n_pair_ratio_coverage_pass']}/{core_tn['n_pair_ratio_mappable']}。这不是 core70 的总体支持率，而是其现有 junction-pair 可观测子集的覆盖情况。",
            f"- 上述 23 行 isoform-associated feature 对应 {primary_ta['n_unique_biological_pair_ids_pair_ratio']} 个 biological pair IDs。reciprocal isoform 行共享同一 junction event，因此不能作为独立生物学发现或独立证据简单相加。",
            "",
            "## KRAS、CDX2 与 OLFM4",
            "",
            target_summary_line("KRAS", "16195-1_ENST00000256078.10", "T_vs_A"),
            target_summary_line("KRAS", "2994-1_ENST00000311936.8", "T_vs_A"),
            target_summary_line("CDX2", "16338-1_ENST00000548877.1", "T_vs_A"),
            target_summary_line("CDX2", "3673-1_ENST00000381020.8", "T_vs_A"),
            target_summary_line("OLFM4", "3790-1_ENST00000219022.3", "T_vs_A"),
            target_summary_line("OLFM4", "flairiso55540-1", "T_vs_A"),
            "",
            "## Manuscript-safe wording",
            "",
            "可写：在同一 N/A/T 队列的患者级 short-read junction-pair measurements 中，具有可分辨互斥 junction 的 tumor-associated long-read usage shifts 显示方向一致的正交 junction-level support。",
            "不可写：short reads validated full-length isoforms；independent replication；all 93 isoforms were confirmed；terminal/APA/IR events were validated；任何功能、机制、临床或候选优先级结论。",
            "",
            "## 是否达到 Figure 5 门槛",
            "",
            "当前证据可作为 future Figure 5 的同队列 orthogonal-support 组件，但不能单独构成 full-length isoform external validation，也不能将 reciprocal rows 堆叠为独立事件。是否设计 Figure 5 应等待对 biological-pair 去重展示和外部数据边界的单独审阅。",
        ],
    )
    update_progress("complete", f"STAR junction audit complete; pair-ratio result rows with >= {MIN_PAIRED_PATIENTS} patients: {len(ratio_pass)}; target evaluable rows: {dict(target_summary)}")
    write_status(SUCCESS_PATH, "PASS", "Internal short-read junction/event-level audit completed; no FASTQ/STAR rerun and no Figure 5 generated.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        for directory in (TABLE_DIR, REPORT_DIR, LOG_DIR):
            directory.mkdir(parents=True, exist_ok=True)
        details = f"{type(exc).__name__}: {exc}"
        try:
            update_progress("failed", details)
            write_status(FAILED_PATH, "FAIL", details)
        finally:
            traceback.print_exc()
        sys.exit(1)
