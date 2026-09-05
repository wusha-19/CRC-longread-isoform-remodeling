#!/usr/bin/env python3
import datetime as dt
import json
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ["CRC_ROOT"]).expanduser().resolve()
PHASE00 = ROOT / "05_isoform_main_study/00_metadata_qc"
OUT = ROOT / "05_isoform_main_study/01_master_catalog"

STRICT_DIR = ROOT / "FLAIR/03_3_combine_S5_expr_filtered_strict/strict_ge8samples_count2_total20"
SQANTI_DIR = ROOT / "FLAIR/03_6_sqanti3_S5_strict_ge8_count2_total20_fast/04_Final_QC"

COUNT_PATH = STRICT_DIR / "CRC_combined.S5.strict_ge8samples_count2_total20.counts.tsv"
TPM_PATH = STRICT_DIR / "CRC_combined.S5.strict_ge8samples_count2_total20.tpm.tsv"
CLASS_PATH = SQANTI_DIR / "final_rescued_classification.txt"
GTF_PATH = SQANTI_DIR / "final_rescued_corrected.gtf"
FASTA_PATH = SQANTI_DIR / "final_rescued_corrected.fasta"

META_PATH = PHASE00 / "metadata/CRC_77_sample_metadata.tsv"
PRIMARY_LIST = PHASE00 / "metadata/cohort_primary_NAT_72_samples.txt"
SENS_LIST = PHASE00 / "metadata/cohort_sensitivity_NAT_75_samples.txt"
FILTER_IMPACT = PHASE00 / "tables/isoform_filter_cohort_impact.tsv"
FILTER_CAT = PHASE00 / "tables/isoform_filter_cohort_impact_by_sqanti3_category.tsv"
PHASE00_CONCORD = PHASE00 / "tables/isoform_id_concordance.tsv"

EXPECTED_FINAL_N = 66805
EXPECTED_PRIMARY_N = 63623


def mkdirs():
    for sub in ["code", "config", "tables", "matrices", "reports", "logs", "figures", "provenance", "annotation"]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def q(path):
    return "'" + str(path).replace("'", "'\"'\"'") + "'"


def shell(cmd):
    return subprocess.run(cmd, shell=True, text=True, capture_output=True)


def count_lines(path):
    cp = shell(f"wc -l < {q(path)}")
    return int(cp.stdout.strip()) if cp.returncode == 0 and cp.stdout.strip() else None


def file_audit(role, path, id_col=None):
    path = Path(path)
    status = "OK" if path.exists() and path.stat().st_size > 0 else "MISSING"
    out = {
        "role": role,
        "path": str(path),
        "size": path.stat().st_size if path.exists() else 0,
        "modified_time": dt.datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds") if path.exists() else "NA",
        "rows": "NA",
        "columns": "NA",
        "id_count": "NA",
        "sample_count": "NA",
        "format": path.suffix.lstrip(".") or "text",
        "status": status,
    }
    if not path.exists():
        return out
    try:
        with path.open(errors="replace") as f:
            header = f.readline().rstrip("\n").split("\t")
        out["columns"] = len(header)
        out["rows"] = max(0, (count_lines(path) or 1) - 1)
        sample_cols = [c for c in header if re.match(r"^CRC[-_]\d+[A-Z][-_]\d+$", c)]
        out["sample_count"] = len(sample_cols)
        if id_col and id_col in header:
            vals = pd.read_csv(path, sep="\t", usecols=[id_col])[id_col].astype(str)
            out["id_count"] = vals.nunique()
    except Exception:
        if path.suffix in [".gtf", ".fa", ".fasta"]:
            out["rows"] = count_lines(path)
    return out


def read_list(path):
    return [x.strip() for x in Path(path).read_text().splitlines() if x.strip()]


def to_sqanti_id(series):
    return series.astype(str).str.replace(r"_ENSG[0-9A-Za-z_.-]+$", "", regex=True)


def parse_gtf_ids(path):
    ids = set()
    pat = re.compile(r'transcript_id "([^"]+)"')
    with open(path, errors="replace") as f:
        for line in f:
            if "transcript_id" not in line:
                continue
            m = pat.search(line)
            if m:
                ids.add(m.group(1))
    return ids


def parse_fasta_ids(path):
    ids = set()
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def subset_gtf(path, keep_ids, out_path):
    pat = re.compile(r'transcript_id "([^"]+)"')
    n = 0
    with open(path, errors="replace") as inp, open(out_path, "w") as out:
        for line in inp:
            m = pat.search(line)
            if m and m.group(1) in keep_ids:
                out.write(line)
                n += 1
    return n


def subset_fasta(path, keep_ids, out_path):
    nseq = 0
    write = False
    with open(path, errors="replace") as inp, open(out_path, "w") as out:
        for line in inp:
            if line.startswith(">"):
                tid = line[1:].strip().split()[0]
                write = tid in keep_ids
                if write:
                    nseq += 1
                    out.write(line)
            elif write:
                out.write(line)
    return nseq


def concordance_table(sources):
    rows = []
    names = list(sources)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            A, B = sources[a], sources[b]
            rows.append({
                "source_A": a,
                "source_B": b,
                "n_A": len(A),
                "n_B": len(B),
                "n_intersection": len(A & B),
                "n_only_A": len(A - B),
                "n_only_B": len(B - A),
                "status": "OK" if len(A - B) == 0 and len(B - A) == 0 else "CHECK",
                "notes": "2,454 final SQANTI3 isoforms absent from strict count" if {a, b} == {"final_sqanti3_classification", "strict_count_matrix"} else "",
            })
    return pd.DataFrame(rows)


def main():
    mkdirs()
    started = dt.datetime.now().isoformat(timespec="seconds")

    meta = pd.read_csv(META_PATH, sep="\t")
    primary = read_list(PRIMARY_LIST)
    sensitivity = read_list(SENS_LIST)
    primary_set = set(primary)
    sensitivity_set = set(sensitivity)

    cls = pd.read_csv(CLASS_PATH, sep="\t")
    if len(cls) != EXPECTED_FINAL_N:
        raise SystemExit(f"ERROR: final SQANTI3 classification rows={len(cls)} expected={EXPECTED_FINAL_N}")
    cls_ids = set(cls["isoform"].astype(str))

    counts_raw = pd.read_csv(COUNT_PATH, sep="\t")
    counts_raw["sqanti_isoform_id"] = to_sqanti_id(counts_raw["isoform_id"])
    if counts_raw["sqanti_isoform_id"].duplicated().any():
        raise SystemExit("ERROR: strict count sqanti_isoform_id is not unique after ENSG suffix removal")
    sample_cols = [c for c in counts_raw.columns if re.match(r"^CRC[-_]\d+[A-Z][-_]\d+$", c)]
    if len(sample_cols) != 77:
        raise SystemExit(f"ERROR: strict count sample columns={len(sample_cols)} expected=77")

    missing_meta = set(sample_cols) - set(meta["raw_sample_id"])
    if missing_meta:
        raise SystemExit(f"ERROR: sample columns missing from metadata: {sorted(missing_meta)[:10]}")
    if len(primary) != 72:
        raise SystemExit(f"ERROR: primary sample list length={len(primary)} expected=72")
    if len(sensitivity) != 75:
        raise SystemExit(f"ERROR: sensitivity sample list length={len(sensitivity)} expected=75")
    if any(s.endswith("D_01") or re.match(r"^CRC[-_]\d+D[-_]\d+$", s) for s in primary + sensitivity):
        raise SystemExit("ERROR: D sample found in primary/sensitivity list")
    if any(s.startswith("CRC_20") or s.startswith("CRC-20") for s in primary):
        raise SystemExit("ERROR: CRC_20 found in primary list")
    if not any(s.startswith("CRC_20") or s.startswith("CRC-20") for s in sensitivity):
        raise SystemExit("ERROR: CRC_20 missing from sensitivity list")

    filter_df = pd.read_csv(FILTER_IMPACT, sep="\t")
    filter_df["isoform_id"] = filter_df["isoform_id"].astype(str)

    strict_present_ids = set(counts_raw["sqanti_isoform_id"]) & cls_ids
    count_expr = counts_raw[counts_raw["sqanti_isoform_id"].isin(cls_ids)].copy()
    count_expr = count_expr.set_index("sqanti_isoform_id", drop=False)

    # Matrix outputs.
    primary_ids = filter_df.loc[filter_df["pass_primary72"].astype(bool), "isoform_id"].astype(str).tolist()
    sens_ids = filter_df.loc[filter_df["pass_NAT75"].astype(bool), "isoform_id"].astype(str).tolist()
    if len(primary_ids) != EXPECTED_PRIMARY_N:
        raise SystemExit(f"ERROR: primary72 isoforms={len(primary_ids)} expected={EXPECTED_PRIMARY_N}")

    count_matrix = count_expr.set_index("sqanti_isoform_id")
    primary_counts = count_matrix.loc[primary_ids, primary]
    sens_counts = count_matrix.loc[sens_ids, sensitivity]
    if primary_counts.shape != (EXPECTED_PRIMARY_N, 72):
        raise SystemExit(f"ERROR: primary count matrix shape={primary_counts.shape}")
    if sens_counts.shape[1] != 75:
        raise SystemExit(f"ERROR: sensitivity count matrix columns={sens_counts.shape[1]}")
    primary_counts.index.name = "isoform_id"
    sens_counts.index.name = "isoform_id"
    primary_counts.to_csv(OUT / "matrices/counts_primary72_isoforms.tsv", sep="\t")
    sens_counts.to_csv(OUT / "matrices/counts_sensitivity75_isoforms.tsv", sep="\t")

    # TPM outputs if reliable matching matrix exists.
    tpm_generated = False
    if TPM_PATH.exists():
        tpm_raw = pd.read_csv(TPM_PATH, sep="\t")
        if "isoform_id" in tpm_raw.columns:
            tpm_raw["sqanti_isoform_id"] = to_sqanti_id(tpm_raw["isoform_id"])
            if not tpm_raw["sqanti_isoform_id"].duplicated().any():
                tpm_expr = tpm_raw[tpm_raw["sqanti_isoform_id"].isin(cls_ids)].set_index("sqanti_isoform_id")
                tpm_expr.loc[primary_ids, primary].rename_axis("isoform_id").to_csv(OUT / "matrices/tpm_primary72_isoforms.tsv", sep="\t")
                tpm_expr.loc[sens_ids, sensitivity].rename_axis("isoform_id").to_csv(OUT / "matrices/tpm_sensitivity75_isoforms.tsv", sep="\t")
                tpm_generated = True

    meta_cols = ["raw_sample_id", "normalized_sample_id", "patient_id", "tissue_group", "paired_set", "CRC20_status", "adenoma_exploratory"]
    meta_primary = meta.set_index("raw_sample_id").loc[primary, meta_cols[1:]].reset_index().rename(columns={"raw_sample_id": "sample_id"})
    meta_sens = meta.set_index("raw_sample_id").loc[sensitivity, meta_cols[1:]].reset_index().rename(columns={"raw_sample_id": "sample_id"})
    meta_primary.to_csv(OUT / "matrices/metadata_primary72.tsv", sep="\t", index=False)
    meta_sens.to_csv(OUT / "matrices/metadata_sensitivity75.tsv", sep="\t", index=False)

    # Master catalog.
    cat = pd.DataFrame({"isoform_id": cls["isoform"].astype(str)})
    colmap = {
        "associated_gene": "associated_gene",
        "associated_transcript": "associated_transcript",
        "structural_category": "sqanti3_category",
        "subcategory": "sqanti3_subcategory",
        "chrom": "chrom",
        "strand": "strand",
        "length": "length",
        "exons": "exon_count",
        "RTS_stage": "RTS_stage_if_available",
        "all_canonical": "all_canonical",
    }
    cls_tmp = cls.rename(columns={k: v for k, v in colmap.items() if k in cls.columns})
    keep_cols = ["isoform"] + [v for v in colmap.values() if v in cls_tmp.columns]
    cat = cat.merge(cls_tmp[keep_cols], left_on="isoform_id", right_on="isoform", how="left").drop(columns=["isoform"])
    cat["gene_id"] = cat.get("associated_gene", pd.Series(["NA"] * len(cat))).fillna("NA")
    cat["gene_name"] = "NA"
    for missing in ["ORF_length_if_available", "NMD_if_available", "RTS_stage_if_available", "all_canonical"]:
        if missing not in cat.columns:
            cat[missing] = "NA"

    fi = filter_df.rename(columns={
        "support_ge2_all77": "support_ge2_all77",
        "total_count_all77": "strict_count_total_all77",
        "pass_all77": "pass_strict_all77",
        "support_ge2_NAT75": "support_ge2_NAT75",
        "total_count_NAT75": "strict_count_total_NAT75",
        "pass_NAT75": "pass_strict_NAT75",
        "support_ge2_primary72": "support_ge2_primary72",
        "total_count_primary72": "strict_count_total_primary72",
        "pass_primary72": "pass_strict_primary72",
    })
    fi_cols = [
        "isoform_id", "support_ge2_all77", "strict_count_total_all77", "pass_strict_all77",
        "support_ge2_NAT75", "strict_count_total_NAT75", "pass_strict_NAT75",
        "support_ge2_primary72", "strict_count_total_primary72", "pass_strict_primary72",
        "pass_change_due_to_D", "pass_change_due_to_CRC20",
    ]
    fi_cols = [c for c in fi_cols if c in fi.columns]
    cat = cat.merge(fi[fi_cols], on="isoform_id", how="left")
    cat["strict_count_present"] = cat["isoform_id"].isin(strict_present_ids)
    for c in ["support_ge2_all77", "strict_count_total_all77", "support_ge2_NAT75", "strict_count_total_NAT75", "support_ge2_primary72", "strict_count_total_primary72"]:
        cat[c] = cat[c].fillna(0).astype(int)
    for c in ["pass_strict_all77", "pass_strict_NAT75", "pass_strict_primary72", "pass_change_due_to_D", "pass_change_due_to_CRC20"]:
        cat[c] = cat[c].fillna(False).astype(bool)
    cat["min_sample_count_available"] = "NA"
    cat["analysis_set_primary72"] = cat["strict_count_present"] & cat["pass_strict_primary72"]
    cat["analysis_set_sensitivity75"] = cat["strict_count_present"] & cat["pass_strict_NAT75"]
    cat["recommended_for_DTE_primary"] = cat["analysis_set_primary72"]

    primary_gene_counts = cat.loc[cat["analysis_set_primary72"], "gene_id"].value_counts()
    dtu_genes = set(primary_gene_counts[primary_gene_counts >= 2].index)
    cat["recommended_for_DTU_primary"] = cat["analysis_set_primary72"] & cat["gene_id"].isin(dtu_genes)
    cat["notes"] = np.where(~cat["strict_count_present"], "final SQANTI3 isoform absent from strict count after ID normalization; retained in catalog but excluded from expression analysis", "")

    ordered = [
        "isoform_id", "gene_id", "gene_name", "sqanti3_category", "sqanti3_subcategory",
        "associated_transcript", "associated_gene", "chrom", "strand", "length", "exon_count",
        "ORF_length_if_available", "NMD_if_available", "RTS_stage_if_available", "all_canonical",
        "min_sample_count_available", "strict_count_present", "strict_count_total_all77",
        "support_ge2_all77", "pass_strict_all77", "strict_count_total_NAT75",
        "support_ge2_NAT75", "pass_strict_NAT75", "strict_count_total_primary72",
        "support_ge2_primary72", "pass_strict_primary72", "lost_after_removing_D",
        "lost_after_removing_CRC20", "analysis_set_primary72", "analysis_set_sensitivity75",
        "recommended_for_DTE_primary", "recommended_for_DTU_primary", "notes",
    ]
    cat = cat.rename(columns={"pass_change_due_to_D": "lost_after_removing_D", "pass_change_due_to_CRC20": "lost_after_removing_CRC20"})
    for c in ordered:
        if c not in cat.columns:
            cat[c] = "NA"
    cat[ordered].to_csv(OUT / "tables/CRC_isoform_master_catalog.tsv", sep="\t", index=False)
    if len(cat) != EXPECTED_FINAL_N:
        raise SystemExit(f"ERROR: master catalog rows={len(cat)} expected={EXPECTED_FINAL_N}")
    if int(cat["analysis_set_primary72"].sum()) != EXPECTED_PRIMARY_N:
        raise SystemExit(f"ERROR: analysis_set_primary72={int(cat['analysis_set_primary72'].sum())} expected={EXPECTED_PRIMARY_N}")

    # Summaries.
    summary_cat = cat.groupby("sqanti3_category", dropna=False).agg(
        n_final_sqanti3=("isoform_id", "count"),
        n_strict_count_present=("strict_count_present", "sum"),
        n_pass_all77=("pass_strict_all77", "sum"),
        n_pass_NAT75=("pass_strict_NAT75", "sum"),
        n_pass_primary72=("pass_strict_primary72", "sum"),
        n_analysis_set_primary72=("analysis_set_primary72", "sum"),
        n_analysis_set_sensitivity75=("analysis_set_sensitivity75", "sum"),
    ).reset_index()
    summary_cat.to_csv(OUT / "tables/master_catalog_summary_by_sqanti3_category.tsv", sep="\t", index=False)

    summary_gene = cat.groupby(["gene_id", "gene_name"], dropna=False).agg(
        n_final_isoforms=("isoform_id", "count"),
        n_primary72_isoforms=("analysis_set_primary72", "sum"),
        n_sensitivity75_isoforms=("analysis_set_sensitivity75", "sum"),
        eligible_for_DTU_primary=("recommended_for_DTU_primary", "any"),
    ).reset_index()
    summary_gene.to_csv(OUT / "tables/master_catalog_summary_by_gene.tsv", sep="\t", index=False)

    # ID lists and annotation subsets.
    (OUT / "tables/primary72_isoform_ids.txt").write_text("\n".join(primary_ids) + "\n")
    (OUT / "tables/sensitivity75_isoform_ids.txt").write_text("\n".join(sens_ids) + "\n")
    gtf_primary_n = subset_gtf(GTF_PATH, set(primary_ids), OUT / "annotation/primary72_isoforms.gtf")
    gtf_sens_n = subset_gtf(GTF_PATH, set(sens_ids), OUT / "annotation/sensitivity75_isoforms.gtf")
    fa_primary_n = subset_fasta(FASTA_PATH, set(primary_ids), OUT / "annotation/primary72_isoforms.fa")
    fa_sens_n = subset_fasta(FASTA_PATH, set(sens_ids), OUT / "annotation/sensitivity75_isoforms.fa")

    # Concordance.
    gtf_ids = parse_gtf_ids(GTF_PATH)
    fasta_ids = parse_fasta_ids(FASTA_PATH)
    sources = {
        "final_sqanti3_classification": cls_ids,
        "final_gtf": gtf_ids,
        "final_fasta": fasta_ids,
        "strict_count_matrix": strict_present_ids,
        "primary72_matrix": set(primary_ids),
        "sensitivity75_matrix": set(sens_ids),
        "master_catalog": set(cat["isoform_id"]),
    }
    concord = concordance_table(sources)
    concord.to_csv(OUT / "tables/isoform_id_concordance_phase01.tsv", sep="\t", index=False)

    # Input audit.
    audits = [
        file_audit("Phase00 metadata", META_PATH),
        file_audit("Phase00 primary cohort", PRIMARY_LIST),
        file_audit("Phase00 sensitivity cohort", SENS_LIST),
        file_audit("Phase00 filter impact", FILTER_IMPACT, "isoform_id"),
        file_audit("Phase00 filter impact by SQANTI3 category", FILTER_CAT),
        file_audit("Phase00 ID concordance", PHASE00_CONCORD),
        file_audit("final SQANTI3 classification", CLASS_PATH, "isoform"),
        file_audit("final SQANTI3 corrected GTF", GTF_PATH),
        file_audit("final SQANTI3 corrected FASTA", FASTA_PATH),
        file_audit("strict count matrix", COUNT_PATH, "isoform_id"),
        file_audit("strict TPM matrix", TPM_PATH, "isoform_id"),
    ]
    pd.DataFrame(audits).to_csv(OUT / "tables/phase01_input_audit.tsv", sep="\t", index=False)

    # Plot data tables.
    summary_cat.to_csv(OUT / "figures/isoform_retention_by_sqanti3_category.data.tsv", sep="\t", index=False)
    summary_gene.to_csv(OUT / "figures/isoforms_per_gene_primary72.data.tsv", sep="\t", index=False)
    lib = primary_counts.sum(axis=0).reset_index()
    lib.columns = ["sample_id", "library_size_after_filter"]
    lib = lib.merge(meta[["raw_sample_id", "patient_id", "tissue_group", "CRC20_status"]], left_on="sample_id", right_on="raw_sample_id", how="left")
    lib.to_csv(OUT / "figures/primary72_library_size_after_filter.data.tsv", sep="\t", index=False)

    # Reports.
    missing_count = EXPECTED_FINAL_N - len(strict_present_ids)
    report = f"""# Phase 01 master isoform catalog report

## Inputs

See `tables/phase01_input_audit.tsv` for file size, timestamp, row/column counts and status.

Main expression source: `{COUNT_PATH}`.
Main annotation source: `{CLASS_PATH}`.

## Cohort definitions

- Primary cohort: 72 N/A/T samples, 24 patients, excluding CRC_20 and D samples.
- Sensitivity cohort: 75 N/A/T samples, 25 patients, including CRC_20, excluding D samples.
- CRC_20 and D handling were inherited from Phase 00 and not reinterpreted.

## Master catalog

- Master catalog rows: {len(cat)}
- Final SQANTI3 isoforms: {len(cls_ids)}
- Strict count-present final SQANTI3 isoforms: {len(strict_present_ids)}
- Final SQANTI3 isoforms absent from strict count: {missing_count}

## Expression analysis sets

- primary72 expression-analysis isoforms: {int(cat['analysis_set_primary72'].sum())}
- sensitivity75 expression-analysis isoforms: {int(cat['analysis_set_sensitivity75'].sum())}
- recommended DTE primary isoforms: {int(cat['recommended_for_DTE_primary'].sum())}
- recommended DTU primary isoforms: {int(cat['recommended_for_DTU_primary'].sum())}

## ID mismatch handling

Strict count IDs contain an ENSG suffix, whereas final SQANTI3 IDs do not. Phase 01 normalized strict count IDs by removing the trailing `_ENSG...` suffix. After normalization, {len(strict_present_ids)} final SQANTI3 isoforms had expression counts; {missing_count} final SQANTI3 isoforms did not and were retained in the master catalog with `strict_count_present = FALSE`.

## SQANTI3 category distribution

See `tables/master_catalog_summary_by_sqanti3_category.tsv`.

## Re-filtering decision

The structure-level SQANTI3 result is retained. For formal expression analyses, use the Phase 01 primary/sensitivity matrices, which reapply the support filter at the analysis layer according to the Phase 00 locked cohorts.

## GTF/FASTA subsets

Generated:

- `annotation/primary72_isoforms.gtf` ({gtf_primary_n} GTF lines)
- `annotation/sensitivity75_isoforms.gtf` ({gtf_sens_n} GTF lines)
- `annotation/primary72_isoforms.fa` ({fa_primary_n} sequences)
- `annotation/sensitivity75_isoforms.fa` ({fa_sens_n} sequences)

## Risks and limitations

- `gene_name` was not reliably available from the SQANTI3 classification and is set to `NA`.
- ORF and NMD fields were not available in the current SQANTI3 classification and are set to `NA`.
- DTU eligibility here is a preparatory flag based on at least two primary72 isoforms per gene; no DTU model was run.

## Next Phase

Proceed to Phase 02/03 only after accepting these frozen input matrices. Do not use D samples in N/A/T paired models.
"""
    (OUT / "reports/phase01_master_catalog_report.md").write_text(report)

    worklog = f"""# Phase 01 worklog

## Date

- Start: {started}
- End: {dt.datetime.now().isoformat(timespec='seconds')}

## Task goal

Build the master isoform catalog and frozen primary/sensitivity expression inputs without running DTE, DTU, SUPPA2 differential analysis, trend models, or candidate screening.

## Inputs

- `{COUNT_PATH}`
- `{TPM_PATH}`
- `{CLASS_PATH}`
- `{GTF_PATH}`
- `{FASTA_PATH}`
- Phase 00 metadata and cohort/filter-impact tables.

## Outputs

- `tables/CRC_isoform_master_catalog.tsv`
- `matrices/counts_primary72_isoforms.tsv`
- `matrices/counts_sensitivity75_isoforms.tsv`
- `matrices/tpm_primary72_isoforms.tsv` and `matrices/tpm_sensitivity75_isoforms.tsv` if generated: {tpm_generated}
- `matrices/metadata_primary72.tsv`
- `matrices/metadata_sensitivity75.tsv`
- `tables/isoform_id_concordance_phase01.tsv`
- GTF/FASTA subsets in `annotation/`

## Scripts

- `code/phase01_build_master_catalog.py`
- `code/phase01_make_figures.R`

## Commands

`python3 code/phase01_build_master_catalog.py`
`Rscript code/phase01_make_figures.R`

## Software versions

- Python: {sys.version.split()[0]}
- pandas: {pd.__version__}
- numpy: {np.__version__}

## Key parameters

- Expression support threshold inherited from Phase 00: count >= 2 in at least 8 samples and total count >= 20.
- Primary cohort: 72 samples.
- Sensitivity cohort: 75 samples.

## Run status

SUCCESS

## Errors and fixes

No upstream files were modified. ID suffix mismatch was handled by removing trailing `_ENSG...` from strict count isoform IDs.

## Result summary

- Master rows: {len(cat)}
- Primary isoforms: {int(cat['analysis_set_primary72'].sum())}
- Sensitivity isoforms: {int(cat['analysis_set_sensitivity75'].sum())}
- Final SQANTI3 absent from strict count: {missing_count}

## Next step

Review reports and frozen matrices before Phase 02/03.
"""
    (OUT / "reports/phase01_worklog.md").write_text(worklog)

    provenance = {
        "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
        "count_path": str(COUNT_PATH),
        "tpm_path": str(TPM_PATH),
        "classification_path": str(CLASS_PATH),
        "primary_samples": primary,
        "sensitivity_samples": sensitivity,
        "master_rows": len(cat),
        "primary_isoforms": int(cat["analysis_set_primary72"].sum()),
        "sensitivity_isoforms": int(cat["analysis_set_sensitivity75"].sum()),
        "strict_count_present_final_sqanti3": len(strict_present_ids),
        "final_sqanti3_absent_from_strict_count": missing_count,
    }
    (OUT / "provenance/phase01_provenance.json").write_text(json.dumps(provenance, indent=2))
    print("PHASE01_MASTER_COMPLETE", len(cat), int(cat["analysis_set_primary72"].sum()), int(cat["analysis_set_sensitivity75"].sum()), missing_count)


if __name__ == "__main__":
    main()
