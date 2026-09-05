#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import socket
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(os.environ["CRC_ROOT"]).expanduser().resolve()
PHASE01 = ROOT / "05_isoform_main_study/01_master_catalog"
PHASE02 = ROOT / "05_isoform_main_study/02_DTE_paired"
PHASE03 = ROOT / "05_isoform_main_study/03_DTU_paired"
PHASE03B = ROOT / "05_isoform_main_study/03B_DTU_review"
OUT = ROOT / "05_isoform_main_study/04_NAT_dynamic_patterns"

CONTRASTS = ["A_vs_N", "T_vs_A", "T_vs_N"]
TISSUES = ["N", "A", "T"]
MAIN_DTU = "DTU_FDR05_delta010"


def ensure_dirs() -> None:
    for sub in ["code", "tables", "figures", "figures/HIGH_DTU_NAT_plots", "reports", "logs", "matrices", "provenance"]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k, "NA")) for k in fieldnames})


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, (float, np.floating)):
        if math.isnan(float(value)):
            return "NA"
        return f"{float(value):.12g}"
    if value == "":
        return "NA"
    return str(value)


def fnum(value: object) -> float:
    try:
        if value in (None, "", "NA"):
            return math.nan
        return float(str(value))
    except ValueError:
        return math.nan


def bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def table_key(rows: list[dict[str, str]], keys: tuple[str, ...]) -> dict[tuple[str, ...], dict[str, str]]:
    return {tuple(row.get(k, "") for k in keys): row for row in rows}


def median_ratio_size_factors(count_matrix: np.ndarray) -> np.ndarray:
    positive = count_matrix > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        log_counts = np.where(positive, np.log(count_matrix), np.nan)
    geo = np.exp(np.nanmean(log_counts, axis=1))
    valid = np.isfinite(geo) & (geo > 0)
    if valid.sum() == 0:
        lib = count_matrix.sum(axis=0)
        med = np.median(lib[lib > 0])
        return lib / med
    ratios = count_matrix[valid, :] / geo[valid, None]
    factors = np.nanmedian(np.where(ratios > 0, ratios, np.nan), axis=0)
    if np.any(~np.isfinite(factors)) or np.any(factors <= 0):
        lib = count_matrix.sum(axis=0)
        med = np.median(lib[lib > 0])
        factors = lib / med
    return factors


def direction_consistency(values: np.ndarray, meta: list[dict[str, str]], group1: str, group0: str) -> float:
    by_patient: dict[str, dict[str, int]] = defaultdict(dict)
    for i, row in enumerate(meta):
        by_patient[row["patient_id"]][row["tissue_group"]] = i
    same = 0
    total = 0
    diffs = []
    for pid, idx in by_patient.items():
        if group1 in idx and group0 in idx:
            diff = float(values[idx[group1]] - values[idx[group0]])
            if not math.isnan(diff) and diff != 0:
                diffs.append(diff)
    if not diffs:
        return math.nan
    mean_diff = float(np.nanmean(diffs))
    if mean_diff == 0 or math.isnan(mean_diff):
        return math.nan
    for diff in diffs:
        total += 1
        same += int((diff > 0 and mean_diff > 0) or (diff < 0 and mean_diff < 0))
    return same / total if total else math.nan


def classify_expression(mean_n: float, mean_a: float, mean_t: float, dte: dict[str, dict[str, str]], cons: dict[str, float]) -> str:
    sig = {}
    lfc = {}
    for c in CONTRASTS:
        row = dte.get(c, {})
        padj = fnum(row.get("padj"))
        lfc[c] = fnum(row.get("log2FoldChange"))
        sig[c] = (not math.isnan(padj)) and padj < 0.05 and (not math.isnan(lfc[c])) and abs(lfc[c]) >= 0.58
    if (not sig["A_vs_N"]) and sig["T_vs_A"] and sig["T_vs_N"] and lfc["T_vs_A"] * lfc["T_vs_N"] > 0:
        return "Tumor_specific_up" if lfc["T_vs_N"] > 0 else "Tumor_specific_down"
    if mean_n < mean_a < mean_t and cons.get("T_vs_N", 0) >= 0.6:
        return "Progressive_up"
    if mean_n > mean_a > mean_t and cons.get("T_vs_N", 0) >= 0.6:
        return "Progressive_down"
    if sig["A_vs_N"] and lfc["A_vs_N"] > 0 and (math.isnan(lfc["T_vs_A"]) or lfc["T_vs_A"] >= -0.25):
        return "Adjacent_primed_up"
    if sig["A_vs_N"] and lfc["A_vs_N"] < 0 and (math.isnan(lfc["T_vs_A"]) or lfc["T_vs_A"] <= 0.25):
        return "Adjacent_primed_down"
    if mean_a > mean_n and mean_a > mean_t and abs(mean_a - max(mean_n, mean_t)) >= 0.25:
        return "Adjacent_peak"
    if mean_a < mean_n and mean_a < mean_t and abs(mean_a - min(mean_n, mean_t)) >= 0.25:
        return "Adjacent_trough"
    return "Stable_or_unclear"


def classify_usage(mean_n: float, mean_a: float, mean_t: float, dtu: dict[str, dict[str, str]], cons: dict[str, float]) -> str:
    sig = {}
    delta = {}
    for c in CONTRASTS:
        row = dtu.get(c, {})
        fdr = fnum(row.get("padj"))
        delta[c] = fnum(row.get("delta_usage"))
        sig[c] = row.get("dtu_status") == MAIN_DTU and (not math.isnan(fdr)) and fdr < 0.05 and abs(delta[c]) >= 0.10
    if (not sig["A_vs_N"]) and sig["T_vs_A"] and sig["T_vs_N"] and delta["T_vs_A"] * delta["T_vs_N"] > 0:
        return "Tumor_specific_up" if delta["T_vs_N"] > 0 else "Tumor_specific_down"
    if mean_n < mean_a < mean_t and cons.get("T_vs_N", 0) >= 0.6:
        return "Progressive_up"
    if mean_n > mean_a > mean_t and cons.get("T_vs_N", 0) >= 0.6:
        return "Progressive_down"
    if sig["A_vs_N"] and delta["A_vs_N"] > 0 and (math.isnan(delta["T_vs_A"]) or delta["T_vs_A"] >= -0.05):
        return "Adjacent_primed_up"
    if sig["A_vs_N"] and delta["A_vs_N"] < 0 and (math.isnan(delta["T_vs_A"]) or delta["T_vs_A"] <= 0.05):
        return "Adjacent_primed_down"
    if mean_a > mean_n and mean_a > mean_t and abs(mean_a - max(mean_n, mean_t)) >= 0.05:
        return "Adjacent_peak"
    if mean_a < mean_n and mean_a < mean_t and abs(mean_a - min(mean_n, mean_t)) >= 0.05:
        return "Adjacent_trough"
    return "Stable_or_unclear"


def load_counts(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[1:]
        ids = []
        data = []
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            ids.append(parts[0])
            data.append([float(x) for x in parts[1:]])
    return ids, samples, np.asarray(data, dtype=float)


def plot_bar(rows: list[dict[str, object]], group_key: str, value_key: str, out_pdf: Path, title: str) -> None:
    labels = [str(r[group_key]) for r in rows]
    values = [float(r[value_key]) for r in rows]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels, values, color="#4C78A8")
    ax.set_title(title)
    ax.set_ylabel("Isoforms")
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_overlap(n_dte: int, n_dtu: int, n_both: int, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 4))
    labels = ["DTE only", "DTE+DTU", "DTU only"]
    values = [n_dte - n_both, n_both, n_dtu - n_both]
    ax.bar(labels, values, color=["#72B7B2", "#F58518", "#E45756"])
    ax.set_ylabel("Isoforms")
    ax.set_title("DTE and DTU overlap")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_heatmap(matrix: np.ndarray, row_labels: list[str], col_labels: list[str], out_pdf: Path, title: str, cmap: str) -> None:
    fig_h = max(4, 0.28 * len(row_labels) + 1.5)
    fig, ax = plt.subplots(figsize=(5.5, fig_h))
    im = ax.imshow(matrix, aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(col_labels)), col_labels)
    ax.set_yticks(range(len(row_labels)), row_labels, fontsize=6)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def paired_line_plot(event: dict[str, str], expr_values: np.ndarray, usage_values: np.ndarray, meta: list[dict[str, str]], out_pdf: Path) -> None:
    by_patient = defaultdict(dict)
    for i, row in enumerate(meta):
        by_patient[row["patient_id"]][row["tissue_group"]] = i
    x = np.arange(3)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4), sharex=True)
    for ax, values, ylabel, title in [
        (axes[0], expr_values, "log2 normalized expression", "Expression"),
        (axes[1], usage_values, "isoform usage", "Usage"),
    ]:
        for pid, idx in by_patient.items():
            if all(t in idx for t in TISSUES):
                y = [values[idx[t]] for t in TISSUES]
                ax.plot(x, y, color="#9A9A9A", alpha=0.45, linewidth=0.8)
                ax.scatter(x, y, color="#4C78A8", s=12, alpha=0.75)
        means = [np.nanmean([values[idx[t]] for idx in by_patient.values() if t in idx]) for t in TISSUES]
        ax.plot(x, means, color="#E45756", marker="o", linewidth=2.2)
        ax.set_xticks(x, TISSUES)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    fig.suptitle(
        f"{event['gene_id']} | {event['isoform_id']} | {event['contrast']} | "
        f"gene-silent={event.get('gene_silent_DTU')} | dominant={event.get('dominant_isoform_switch')} | replicated={event.get('replicated_in_sensitivity75')}",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    inputs = [
        PHASE01 / "matrices/counts_primary72_isoforms.tsv",
        PHASE01 / "matrices/metadata_primary72.tsv",
        PHASE01 / "tables/CRC_isoform_master_catalog.tsv",
        PHASE03B / "tables/DTU_pre_candidate_review.tsv",
        PHASE03B / "tables/DTU_high_priority_events.tsv",
    ]
    for contrast in CONTRASTS:
        inputs.append(PHASE02 / "tables" / f"primary72_DTE_{contrast}_all.tsv")
        inputs.append(PHASE03 / "tables" / f"primary72_DTU_isoform_{contrast}.tsv")
    inputs.extend(
        [
            PHASE03 / "tables/primary72_gene_silent_DTU.tsv",
            PHASE03 / "tables/primary72_dominant_isoform_switch.tsv",
            PHASE03 / "tables/primary72_DTU_patient_direction_consistency.tsv",
        ]
    )
    for path in inputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Required input missing or empty: {path}")

    meta = read_tsv(PHASE01 / "matrices/metadata_primary72.tsv")
    counts_ids, samples, counts = load_counts(PHASE01 / "matrices/counts_primary72_isoforms.tsv")
    if samples != [r["sample_id"] for r in meta]:
        raise ValueError("Sample order mismatch between primary72 count matrix and metadata")
    if any(r["tissue_group"] == "D" for r in meta):
        raise ValueError("primary72 metadata contains D samples")
    if any(r.get("CRC20_status") == "CRC20" for r in meta):
        raise ValueError("primary72 metadata contains CRC_20")

    catalog_rows = read_tsv(PHASE01 / "tables/CRC_isoform_master_catalog.tsv")
    catalog = {r["isoform_id"]: r for r in catalog_rows}
    factors = median_ratio_size_factors(counts)
    expr = np.log2((counts / factors[None, :]) + 1.0)
    np.savetxt(OUT / "matrices/primary72_isoform_log2_sizefactor_normalized_expression.tsv.gz", expr, delimiter="\t")

    sample_by_tissue = {t: [i for i, r in enumerate(meta) if r["tissue_group"] == t] for t in TISSUES}
    dte_by_iso: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    dte_related = set()
    for c in CONTRASTS:
        for row in read_tsv(PHASE02 / "tables" / f"primary72_DTE_{c}_all.tsv"):
            dte_by_iso[row["isoform_id"]][c] = row
            if fnum(row.get("padj")) < 0.05 and abs(fnum(row.get("log2FoldChange"))) >= 1.0:
                dte_related.add(row["isoform_id"])

    dtu_by_iso: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    dtu_related = set()
    for c in CONTRASTS:
        for row in read_tsv(PHASE03 / "tables" / f"primary72_DTU_isoform_{c}.tsv"):
            dtu_by_iso[row["isoform_id"]][c] = row
            if row.get("dtu_status") == MAIN_DTU:
                dtu_related.add(row["isoform_id"])

    pre_rows = read_tsv(PHASE03B / "tables/DTU_pre_candidate_review.tsv")
    high_rows = read_tsv(PHASE03B / "tables/DTU_high_priority_events.tsv")
    pre_by_iso: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in pre_rows:
        pre_by_iso[row["isoform_id"]].append(row)
    high_ids = {r["isoform_id"] for r in high_rows}
    review_priority_by_iso = {}
    priority_rank = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
    for isoform_id, rows in pre_by_iso.items():
        review_priority_by_iso[isoform_id] = sorted(rows, key=lambda r: priority_rank.get(r.get("review_priority"), 9))[0]
    missing_review_ids = [x for x in dtu_related if x not in review_priority_by_iso]
    if missing_review_ids:
        raise ValueError(f"DTU main isoforms missing Phase03B review rows: {len(missing_review_ids)}")

    id_to_idx = {isoform_id: i for i, isoform_id in enumerate(counts_ids)}
    missing_count_ids = [x for x in high_ids | dte_related | dtu_related if x not in id_to_idx]
    if missing_count_ids:
        raise ValueError(f"Related isoforms missing count matrix: {len(missing_count_ids)}")

    expression_rows = []
    expression_by_iso = {}
    for isoform_id in counts_ids:
        idx = id_to_idx[isoform_id]
        vals = expr[idx, :]
        means = {t: float(np.nanmean(vals[sample_by_tissue[t]])) for t in TISSUES}
        medians = {t: float(np.nanmedian(vals[sample_by_tissue[t]])) for t in TISSUES}
        cons = {
            "A_vs_N": direction_consistency(vals, meta, "A", "N"),
            "T_vs_A": direction_consistency(vals, meta, "T", "A"),
            "T_vs_N": direction_consistency(vals, meta, "T", "N"),
        }
        pattern = classify_expression(means["N"], means["A"], means["T"], dte_by_iso.get(isoform_id, {}), cons)
        cat = catalog.get(isoform_id, {})
        row = {
            "gene_id": cat.get("gene_id", "NA"),
            "gene_name": cat.get("gene_name", "NA"),
            "isoform_id": isoform_id,
            "sqanti3_category": cat.get("sqanti3_category", "NA"),
            "mean_vst_N": means["N"],
            "mean_vst_A": means["A"],
            "mean_vst_T": means["T"],
            "median_vst_N": medians["N"],
            "median_vst_A": medians["A"],
            "median_vst_T": medians["T"],
            "delta_A_minus_N": means["A"] - means["N"],
            "delta_T_minus_A": means["T"] - means["A"],
            "delta_T_minus_N": means["T"] - means["N"],
            "paired_direction_consistency_A_vs_N": cons["A_vs_N"],
            "paired_direction_consistency_T_vs_A": cons["T_vs_A"],
            "paired_direction_consistency_T_vs_N": cons["T_vs_N"],
            "expression_dynamic_pattern": pattern,
        }
        expression_rows.append(row)
        expression_by_iso[isoform_id] = row

    expr_fields = list(expression_rows[0].keys())
    write_tsv(OUT / "tables/primary72_isoform_expression_NAT_dynamic_patterns.tsv", expression_rows, expr_fields)

    dtu_gene_ids = {catalog.get(i, {}).get("gene_id", "") for i in dtu_by_iso.keys()}
    dtu_gene_ids.discard("")
    gene_counts = defaultdict(lambda: np.zeros(len(samples), dtype=float))
    for i, isoform_id in enumerate(counts_ids):
        gid = catalog.get(isoform_id, {}).get("gene_id", "")
        if gid in dtu_gene_ids:
            gene_counts[gid] += counts[i, :]

    usage_rows = []
    usage_by_iso = {}
    for isoform_id in sorted(dtu_by_iso.keys()):
        idx = id_to_idx[isoform_id]
        gid = catalog.get(isoform_id, {}).get("gene_id", "")
        denom = gene_counts[gid]
        with np.errstate(divide="ignore", invalid="ignore"):
            usage = np.where(denom > 0, counts[idx, :] / denom, np.nan)
        if np.all(~np.isfinite(usage)):
            raise ValueError(f"All usage values are NA for {isoform_id}")
        means = {t: float(np.nanmean(usage[sample_by_tissue[t]])) for t in TISSUES}
        medians = {t: float(np.nanmedian(usage[sample_by_tissue[t]])) for t in TISSUES}
        cons = {
            "A_vs_N": direction_consistency(usage, meta, "A", "N"),
            "T_vs_A": direction_consistency(usage, meta, "T", "A"),
            "T_vs_N": direction_consistency(usage, meta, "T", "N"),
        }
        pattern = classify_usage(means["N"], means["A"], means["T"], dtu_by_iso.get(isoform_id, {}), cons)
        cat = catalog.get(isoform_id, {})
        row = {
            "gene_id": gid,
            "gene_name": cat.get("gene_name", "NA"),
            "isoform_id": isoform_id,
            "sqanti3_category": cat.get("sqanti3_category", "NA"),
            "mean_usage_N": means["N"],
            "mean_usage_A": means["A"],
            "mean_usage_T": means["T"],
            "median_usage_N": medians["N"],
            "median_usage_A": medians["A"],
            "median_usage_T": medians["T"],
            "delta_usage_A_minus_N": means["A"] - means["N"],
            "delta_usage_T_minus_A": means["T"] - means["A"],
            "delta_usage_T_minus_N": means["T"] - means["N"],
            "usage_direction_consistency_A_vs_N": cons["A_vs_N"],
            "usage_direction_consistency_T_vs_A": cons["T_vs_A"],
            "usage_direction_consistency_T_vs_N": cons["T_vs_N"],
            "usage_dynamic_pattern": pattern,
        }
        usage_rows.append(row)
        usage_by_iso[isoform_id] = row
    usage_fields = list(usage_rows[0].keys())
    write_tsv(OUT / "tables/primary72_isoform_usage_NAT_dynamic_patterns.tsv", usage_rows, usage_fields)

    integrated_ids = sorted(dte_related | dtu_related)
    integrated_rows = []
    for isoform_id in integrated_ids:
        cat = catalog.get(isoform_id, {})
        expr_row = expression_by_iso.get(isoform_id, {})
        usage_row = usage_by_iso.get(isoform_id, {})
        dte = dte_by_iso.get(isoform_id, {})
        dtu = dtu_by_iso.get(isoform_id, {})
        pre_best = review_priority_by_iso.get(isoform_id, {})
        dte_contrasts = [c for c in CONTRASTS if fnum(dte.get(c, {}).get("padj")) < 0.05 and abs(fnum(dte.get(c, {}).get("log2FoldChange"))) >= 1.0]
        dtu_contrasts = [c for c in CONTRASTS if dtu.get(c, {}).get("dtu_status") == MAIN_DTU]
        replicated = any(r.get("replicated_in_sensitivity75") == "TRUE" for r in pre_by_iso.get(isoform_id, []))
        gene_silent = any(r.get("gene_silent_DTU") == "TRUE" for r in pre_by_iso.get(isoform_id, []))
        dominant_switch = any(r.get("dominant_isoform_switch") == "TRUE" for r in pre_by_iso.get(isoform_id, []))
        review_pri = pre_best.get("review_priority", "NA")
        if review_pri == "HIGH" or ((gene_silent or dominant_switch) and replicated):
            next_pri = "HIGH"
        elif dtu_contrasts or (dte_contrasts and str(expr_row.get("expression_dynamic_pattern", "")).startswith("Tumor_specific")):
            next_pri = "MEDIUM"
        else:
            next_pri = "LOW"
        notes = []
        if review_pri == "NA" and dtu_contrasts:
            notes.append("DTU main event lacks Phase03B priority row")
        if not dtu_contrasts and dte_contrasts:
            notes.append("DTE-only related isoform")
        if dtu_contrasts and not dte_contrasts:
            notes.append("DTU-only related isoform")
        row = {
            "gene_id": cat.get("gene_id", "NA"),
            "gene_name": cat.get("gene_name", "NA"),
            "isoform_id": isoform_id,
            "sqanti3_category": cat.get("sqanti3_category", "NA"),
            "in_DTE": bool_text(bool(dte_contrasts)),
            "in_DTU": bool_text(bool(dtu_contrasts)),
            "DTE_contrasts": ",".join(dte_contrasts) if dte_contrasts else "NA",
            "DTU_contrasts": ",".join(dtu_contrasts) if dtu_contrasts else "NA",
            "expression_dynamic_pattern": expr_row.get("expression_dynamic_pattern", "NA"),
            "usage_dynamic_pattern": usage_row.get("usage_dynamic_pattern", "NA"),
            "primary72_review_priority": review_pri,
            "replicated_in_sensitivity75": bool_text(replicated),
            "gene_silent_DTU": bool_text(gene_silent),
            "dominant_isoform_switch": bool_text(dominant_switch),
            "patient_direction_consistency_expression": expr_row.get("paired_direction_consistency_T_vs_N", "NA"),
            "patient_direction_consistency_usage": usage_row.get("usage_direction_consistency_T_vs_N", "NA"),
            "mean_vst_N": expr_row.get("mean_vst_N", "NA"),
            "mean_vst_A": expr_row.get("mean_vst_A", "NA"),
            "mean_vst_T": expr_row.get("mean_vst_T", "NA"),
            "mean_usage_N": usage_row.get("mean_usage_N", "NA"),
            "mean_usage_A": usage_row.get("mean_usage_A", "NA"),
            "mean_usage_T": usage_row.get("mean_usage_T", "NA"),
            "A_vs_N_DTE_padj": dte.get("A_vs_N", {}).get("padj", "NA"),
            "T_vs_A_DTE_padj": dte.get("T_vs_A", {}).get("padj", "NA"),
            "T_vs_N_DTE_padj": dte.get("T_vs_N", {}).get("padj", "NA"),
            "A_vs_N_DTE_log2FC": dte.get("A_vs_N", {}).get("log2FoldChange", "NA"),
            "T_vs_A_DTE_log2FC": dte.get("T_vs_A", {}).get("log2FoldChange", "NA"),
            "T_vs_N_DTE_log2FC": dte.get("T_vs_N", {}).get("log2FoldChange", "NA"),
            "A_vs_N_DTU_FDR": dtu.get("A_vs_N", {}).get("padj", "NA"),
            "T_vs_A_DTU_FDR": dtu.get("T_vs_A", {}).get("padj", "NA"),
            "T_vs_N_DTU_FDR": dtu.get("T_vs_N", {}).get("padj", "NA"),
            "A_vs_N_delta_usage": dtu.get("A_vs_N", {}).get("delta_usage", "NA"),
            "T_vs_A_delta_usage": dtu.get("T_vs_A", {}).get("delta_usage", "NA"),
            "T_vs_N_delta_usage": dtu.get("T_vs_N", {}).get("delta_usage", "NA"),
            "next_stage_annotation_priority": next_pri,
            "notes": "; ".join(notes) if notes else "NA",
        }
        integrated_rows.append(row)

    integrated_fields = list(integrated_rows[0].keys())
    write_tsv(OUT / "tables/CRC_isoform_integrated_NAT_DTE_DTU_patterns.tsv", integrated_rows, integrated_fields)

    summary_rows = []
    for feature_type, rows, key in [
        ("expression", expression_rows, "expression_dynamic_pattern"),
        ("usage", usage_rows, "usage_dynamic_pattern"),
        ("integrated", integrated_rows, "expression_dynamic_pattern"),
    ]:
        patterns = sorted({str(r.get(key, "NA")) for r in rows})
        for pattern in patterns:
            sub = [r for r in rows if str(r.get(key, "NA")) == pattern]
            iso_ids = {str(r["isoform_id"]) for r in sub}
            genes = {str(r.get("gene_id", "NA")) for r in sub}
            summary_rows.append(
                {
                    "feature_type": feature_type,
                    "pattern": pattern,
                    "n_isoforms": len(iso_ids),
                    "n_genes": len(genes - {"NA", ""}),
                    "n_DTE": len(iso_ids & dte_related),
                    "n_DTU": len(iso_ids & dtu_related),
                    "n_gene_silent_DTU": sum(any(x.get("gene_silent_DTU") == "TRUE" for x in pre_by_iso.get(i, [])) for i in iso_ids),
                    "n_dominant_switch": sum(any(x.get("dominant_isoform_switch") == "TRUE" for x in pre_by_iso.get(i, [])) for i in iso_ids),
                    "n_HIGH_review": len(iso_ids & high_ids),
                }
            )
    summary_fields = ["feature_type", "pattern", "n_isoforms", "n_genes", "n_DTE", "n_DTU", "n_gene_silent_DTU", "n_dominant_switch", "n_HIGH_review"]
    write_tsv(OUT / "tables/NAT_dynamic_pattern_summary.tsv", summary_rows, summary_fields)

    next_rows = [r for r in integrated_rows if r["next_stage_annotation_priority"] in {"HIGH", "MEDIUM"}]
    write_tsv(OUT / "tables/NAT_next_stage_annotation_input.tsv", next_rows, integrated_fields)

    # Figure data and figures.
    expr_fig = [r for r in summary_rows if r["feature_type"] == "expression"]
    usage_fig = [r for r in summary_rows if r["feature_type"] == "usage"]
    write_tsv(OUT / "figures/NAT_expression_pattern_summary.data.tsv", expr_fig, summary_fields)
    write_tsv(OUT / "figures/NAT_usage_pattern_summary.data.tsv", usage_fig, summary_fields)
    plot_bar(expr_fig, "pattern", "n_isoforms", OUT / "figures/NAT_expression_pattern_summary.pdf", "Expression dynamic patterns")
    plot_bar(usage_fig, "pattern", "n_isoforms", OUT / "figures/NAT_usage_pattern_summary.pdf", "Usage dynamic patterns")

    n_both = len(dte_related & dtu_related)
    overlap_rows = [
        {"category": "DTE_only", "n_isoforms": len(dte_related - dtu_related)},
        {"category": "DTE_and_DTU", "n_isoforms": n_both},
        {"category": "DTU_only", "n_isoforms": len(dtu_related - dte_related)},
    ]
    write_tsv(OUT / "figures/DTE_DTU_overlap_summary.data.tsv", overlap_rows, ["category", "n_isoforms"])
    plot_overlap(len(dte_related), len(dtu_related), n_both, OUT / "figures/DTE_DTU_overlap_summary.pdf")

    high_order = [r["isoform_id"] for r in high_rows]
    high_labels = [f"{r['gene_id']}|{r['isoform_id']}|{r['contrast']}" for r in high_rows]
    expr_heat = np.asarray([[expression_by_iso[i]["mean_vst_N"], expression_by_iso[i]["mean_vst_A"], expression_by_iso[i]["mean_vst_T"]] for i in high_order], dtype=float)
    usage_heat = np.asarray([[usage_by_iso[i]["mean_usage_N"], usage_by_iso[i]["mean_usage_A"], usage_by_iso[i]["mean_usage_T"]] for i in high_order], dtype=float)
    write_tsv(
        OUT / "figures/HIGH_DTU_NAT_expression_heatmap.data.tsv",
        [{"event": lab, "mean_vst_N": row[0], "mean_vst_A": row[1], "mean_vst_T": row[2]} for lab, row in zip(high_labels, expr_heat)],
        ["event", "mean_vst_N", "mean_vst_A", "mean_vst_T"],
    )
    write_tsv(
        OUT / "figures/HIGH_DTU_NAT_usage_heatmap.data.tsv",
        [{"event": lab, "mean_usage_N": row[0], "mean_usage_A": row[1], "mean_usage_T": row[2]} for lab, row in zip(high_labels, usage_heat)],
        ["event", "mean_usage_N", "mean_usage_A", "mean_usage_T"],
    )
    plot_heatmap(expr_heat, high_labels, TISSUES, OUT / "figures/HIGH_DTU_NAT_expression_heatmap.pdf", "HIGH DTU expression across N/A/T", "viridis")
    plot_heatmap(usage_heat, high_labels, TISSUES, OUT / "figures/HIGH_DTU_NAT_usage_heatmap.pdf", "HIGH DTU usage across N/A/T", "magma")

    for event in high_rows:
        isoform_id = event["isoform_id"]
        idx = id_to_idx[isoform_id]
        gid = event["gene_id"]
        denom = gene_counts[gid]
        usage_vals = np.where(denom > 0, counts[idx, :] / denom, np.nan)
        out_pdf = OUT / "figures/HIGH_DTU_NAT_plots" / f"{event['contrast']}_{gid}_{isoform_id.replace('/', '_')}.pdf"
        paired_line_plot(event, expr[idx, :], usage_vals, meta, out_pdf)
    high_plot_manifest = [
        {"plot_file": str(p), "file_size_bytes": p.stat().st_size, "non_empty": bool_text(p.stat().st_size > 0)}
        for p in sorted((OUT / "figures/HIGH_DTU_NAT_plots").glob("*.pdf"))
    ]
    write_tsv(OUT / "tables/HIGH_DTU_NAT_plot_manifest.tsv", high_plot_manifest, ["plot_file", "file_size_bytes", "non_empty"])

    # Validation.
    if len(next_rows) == 0 or len(integrated_rows) == 0 or len(expression_rows) == 0 or len(usage_rows) == 0:
        raise ValueError("One or more required output tables would be empty")
    missing_high_next = sorted(high_ids - {r["isoform_id"] for r in next_rows})
    if missing_high_next:
        raise ValueError(f"HIGH events missing next-stage annotation input: {len(missing_high_next)}")

    report = f"""# Phase04 N-A-T dynamic pattern and integrated isoform remodeling map

## Inputs

- Phase01 primary72 counts and metadata: `{PHASE01}`
- Phase02 DTE tables: `{PHASE02 / 'tables'}`
- Phase03 DTU tables: `{PHASE03 / 'tables'}`
- Phase03B review tables: `{PHASE03B / 'tables'}`

## Methods

This phase did not rerun DTE or DTU. A complete Phase02 VST matrix was not available; only PCA VST coordinates were present. Therefore expression dynamics were computed as descriptive `log2(size-factor normalized count + 1)` values using DESeq-style median-ratio size factors. This is used only for N/A/T state-pattern description, not for new differential testing.

Isoform usage was calculated as `isoform_count / gene_total_count` within each sample for DTU-tested genes. Gene totals equal the sum of retained isoform counts for each gene in the primary72 master count matrix.

N/A/T patterns are cross-sectional ordered tissue-state patterns, not true temporal evolution.

## Pattern Definitions

- Tumor-specific: A_vs_N weak or not significant, with T_vs_A and T_vs_N significant and direction-concordant.
- Progressive: ordered means N < A < T or N > A > T with patient direction support.
- Adjacent-primed: A_vs_N already changed and T remains in a similar direction or close to A.
- Adjacent peak/trough: A higher or lower than both N and T.
- Stable_or_unclear: insufficient or conflicting evidence.

## Expression Dynamic Results

```text
{(OUT / 'figures/NAT_expression_pattern_summary.data.tsv').read_text()}
```

## Usage Dynamic Results

```text
{(OUT / 'figures/NAT_usage_pattern_summary.data.tsv').read_text()}
```

## DTE and DTU Integration

```text
{(OUT / 'figures/DTE_DTU_overlap_summary.data.tsv').read_text()}
```

Integrated table: `tables/CRC_isoform_integrated_NAT_DTE_DTU_patterns.tsv`

## HIGH DTU Events

HIGH DTU event heatmaps are saved as:

- `figures/HIGH_DTU_NAT_usage_heatmap.pdf`
- `figures/HIGH_DTU_NAT_expression_heatmap.pdf`

Per-event paired N/A/T expression and usage plots are in `figures/HIGH_DTU_NAT_plots/`.

## Next-stage Annotation Input

`tables/NAT_next_stage_annotation_input.tsv` contains all HIGH and MEDIUM events for later ORF/NMD/domain annotation. It is a prioritization input, not final wet-lab candidate selection.

## Boundaries

This phase does not perform ORF/NMD/domain annotation, SUPPA2, RT-PCR primer design, external data download, or final experimental candidate ranking.

## Next Step

Review HIGH and MEDIUM events, then start a separate ORF/NMD/domain annotation phase if approved.
"""
    (OUT / "reports/phase04_NAT_dynamic_report.md").write_text(report)

    worklog = f"""# Phase04 worklog

## 日期

{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## 任务目标

Build N/A/T ordered tissue-state dynamic patterns for isoform expression and usage, then integrate DTE, DTU and Phase03B review priority.

## 输入

Phase01 counts/metadata/catalog, Phase02 DTE, Phase03 DTU, Phase03B review tables.

## 输出

`{OUT}`

## 脚本

`code/phase04_nat_dynamic_patterns.py`

## 命令

```bash
cd {OUT}
python3 code/phase04_nat_dynamic_patterns.py 2>&1 | tee logs/phase04_nat_dynamic_patterns.log
```

## 软件环境

Python `{subprocess.check_output(['python3', '--version'], text=True).strip()}`; numpy `{np.__version__}`; matplotlib `{matplotlib.__version__}`.

## 参数

DTE significance: padj < 0.05 and |log2FC| >= 1 for DTE-related isoforms. DTU significance: Phase03 main threshold status `{MAIN_DTU}`. Expression descriptive transform: log2(size-factor normalized count + 1).

## 运行状态

SUCCESS

## 错误和修复

Complete Phase02 VST matrix was absent, so a descriptive size-factor normalized log2 expression matrix was generated for pattern mapping. No DTE/DTU model was rerun.

## 下一步

Manual review of HIGH/MEDIUM events before a separate ORF/NMD/domain annotation phase.
"""
    (OUT / "reports/phase04_worklog.md").write_text(worklog)
    provenance = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "phase": "Phase04_NAT_dynamic_patterns",
        "inputs": [str(p) for p in inputs],
        "n_expression_isoforms": len(expression_rows),
        "n_usage_isoforms": len(usage_rows),
        "n_integrated_isoforms": len(integrated_rows),
        "n_next_stage_annotation_input": len(next_rows),
        "n_high_events": len(high_rows),
        "expression_transform": "log2(size-factor normalized count + 1); Phase02 complete VST matrix unavailable",
    }
    (OUT / "provenance/phase04_provenance.json").write_text(json.dumps(provenance, indent=2))

    note = ROOT / "notebook/2026-06-25_phase04_NAT_dynamic_patterns"
    (note / "logs").mkdir(parents=True, exist_ok=True)
    (note / "figures").mkdir(parents=True, exist_ok=True)
    (note / "README.md").write_text(report)
    (note / "commands.sh").write_text(
        f"#!/usr/bin/env bash\nset -euo pipefail\ncd {OUT}\npython3 code/phase04_nat_dynamic_patterns.py 2>&1 | tee logs/phase04_nat_dynamic_patterns.log\n"
    )
    (note / "inputs.tsv").write_text("role\tpath\n" + "\n".join(f"input\t{p}" for p in inputs) + "\n")
    (note / "outputs.tsv").write_text(
        "role\tpath\n"
        f"integrated\t{OUT / 'tables/CRC_isoform_integrated_NAT_DTE_DTU_patterns.tsv'}\n"
        f"next_stage\t{OUT / 'tables/NAT_next_stage_annotation_input.tsv'}\n"
        f"summary\t{OUT / 'tables/NAT_dynamic_pattern_summary.tsv'}\n"
        f"report\t{OUT / 'reports/phase04_NAT_dynamic_report.md'}\n"
    )
    (note / "summary.tsv").write_text((OUT / "tables/NAT_dynamic_pattern_summary.tsv").read_text())

    print("PHASE04_NAT_DYNAMIC_COMPLETE")
    print(f"EXPRESSION_ROWS={len(expression_rows)}")
    print(f"USAGE_ROWS={len(usage_rows)}")
    print(f"INTEGRATED_ROWS={len(integrated_rows)}")
    print(f"NEXT_STAGE_ROWS={len(next_rows)}")
    print(f"HIGH_PLOTS={len(high_plot_manifest)}")


if __name__ == "__main__":
    main()
