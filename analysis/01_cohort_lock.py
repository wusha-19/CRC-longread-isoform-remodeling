#!/usr/bin/env python3
import csv
import datetime as dt
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ["CRC_ROOT"]).expanduser().resolve()
OUT_ROOT = ROOT / "05_isoform_main_study"
PHASE = OUT_ROOT / "00_metadata_qc"
RUNS = OUT_ROOT / "runs"
STAMP = dt.datetime.now().strftime("%Y%m%d_%H%M%S_phase00_cohort_lock")
RUN_DIR = RUNS / STAMP

COUNT_PATH = ROOT / "FLAIR/03_3_combine_S5_expr_filtered_strict/strict_ge8samples_count2_total20/CRC_combined.S5.strict_ge8samples_count2_total20.counts.tsv"
TPM_PATH = ROOT / "FLAIR/03_3_combine_S5_expr_filtered_strict/strict_ge8samples_count2_total20/CRC_combined.S5.strict_ge8samples_count2_total20.tpm.tsv"
META_HINT = ROOT / "FLAIR/03_7_expression_with_sqanti3_fixed/sample_metadata.with_quantify_batch.tsv"
STRICT_INPUT = ROOT / "FLAIR/03_3_combine_S5_expr_filtered_strict"
SQANTI_DIR = ROOT / "FLAIR/03_6_sqanti3_S5_strict_ge8_count2_total20_fast"
CLASS_PATH = SQANTI_DIR / "04_Final_QC/final_rescued_classification.txt"
GTF_PATH = SQANTI_DIR / "04_Final_QC/final_rescued_corrected.gtf"
FASTA_PATH = SQANTI_DIR / "04_Final_QC/final_rescued_corrected.fasta"
FILTER_LIST = SQANTI_DIR / "02_Filter_rules/combined_rules_filtered_inclusion-list.txt"

KEYWORDS = [
    "count", "matrix", "isoform", "classification", "sqanti", "gtf", "fasta",
    "fa", "metadata", "sample", "PCA", "pca", "correlation", "DESeq2", "QC",
    "qc", "log", "junction",
]


def mkdirs():
    for sub in ["code", "metadata", "tables", "figures", "reports", "logs", "provenance"]:
        (PHASE / sub).mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)


def shell(cmd):
    return subprocess.run(cmd, shell=True, text=True, capture_output=True)


def file_stats(path: Path):
    st = path.stat()
    size = st.st_size
    mtime = dt.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")
    lines = None
    cols = None
    header = ""
    if path.is_file() and size < 5_000_000_000:
        try:
            with path.open("r", errors="replace") as f:
                first = f.readline().rstrip("\n")
                header = first[:500]
                if "\t" in first:
                    cols = len(first.split("\t"))
                elif "," in first:
                    cols = len(first.split(","))
            cp = shell(f"wc -l < {q(path)}")
            if cp.returncode == 0:
                lines = int(cp.stdout.strip() or 0)
        except Exception:
            pass
    return size, mtime, lines, cols, header


def q(path):
    return "'" + str(path).replace("'", "'\"'\"'") + "'"


def limited_audit():
    roots = [ROOT / "FLAIR", ROOT / "script"]
    rows = []
    key_re = re.compile("|".join(re.escape(k) for k in KEYWORDS), re.I)
    role_rules = [
        ("final isoform count matrix", COUNT_PATH),
        ("isoform expression TPM matrix", TPM_PATH),
        ("final high-confidence GTF", GTF_PATH),
        ("isoform FASTA", FASTA_PATH),
        ("SQANTI3 classification", CLASS_PATH),
        ("filtered isoform inclusion list", FILTER_LIST),
        ("sample metadata", META_HINT),
    ]
    seen = set()
    for role, path in role_rules:
        if path.exists():
            size, mtime, lines, cols, header = file_stats(path)
            rows.append({
                "file_role": role, "path": str(path), "size_bytes": size,
                "modified_time": mtime, "line_count": lines, "column_count": cols,
                "sample_count": "", "isoform_count": "", "header": header,
                "expected_use": role, "same_analysis_version": "YES",
            })
            seen.add(str(path))
    for root in roots:
        if not root.exists():
            continue
        base_depth = len(root.parts)
        for dirpath, dirnames, filenames in os.walk(root):
            p = Path(dirpath)
            depth = len(p.parts) - base_depth
            if depth >= 3:
                dirnames[:] = []
            for name in filenames:
                if not key_re.search(name):
                    continue
                path = p / name
                if str(path) in seen:
                    continue
                try:
                    size, mtime, lines, cols, header = file_stats(path)
                except Exception:
                    continue
                rows.append({
                    "file_role": "limited keyword audit hit", "path": str(path),
                    "size_bytes": size, "modified_time": mtime, "line_count": lines,
                    "column_count": cols, "sample_count": "", "isoform_count": "",
                    "header": header, "expected_use": "trace/provenance or optional QC",
                    "same_analysis_version": "REVIEW",
                })
                seen.add(str(path))
    df = pd.DataFrame(rows)
    df.to_csv(PHASE / "tables/input_file_audit.tsv", sep="\t", index=False)
    return df


def detect_sample_cols(columns):
    return [c for c in columns if re.match(r"^CRC[-_]\d+[A-Z][-_]\d+$", str(c))]


def parse_sample(sample):
    m = re.match(r"^CRC[-_](\d+)([A-Z])[-_]([0-9]+)$", sample)
    if not m:
        return None
    patient = int(m.group(1))
    code = m.group(2)
    replicate = m.group(3)
    # Existing project metadata maps A=Tumor, B=Adjacent, C=Normal.
    map_code = {"A": "T", "B": "A", "C": "N", "D": "D"}
    tissue = map_code.get(code, "UNRESOLVED")
    return patient, code, replicate, tissue


def build_metadata(sample_cols):
    hint = pd.read_csv(META_HINT, sep="\t") if META_HINT.exists() else pd.DataFrame()
    hint = hint.set_index("sample_id") if "sample_id" in hint.columns else pd.DataFrame()
    rows = []
    for s in sample_cols:
        parsed = parse_sample(s)
        if parsed is None:
            patient, code, rep, tissue = "", "", "", "UNRESOLVED"
            note = "cannot parse sample name"
        else:
            patient, code, rep, tissue = parsed
            note = "A/B/C code mapped using existing metadata convention: A=Tumor, B=Adjacent, C=Normal"
        hint_key = re.sub(r"^CRC[-_](\d+)([A-Z])[-_](\d+)$", r"CRC-\1\2-\3", s)
        if hint_key in hint.index and "ordered_group" in hint.columns:
            expected = hint.loc[hint_key, "ordered_group"]
            if tissue != expected:
                note += f"; metadata ordered_group={expected}"
        primary_nat = tissue in {"N", "A", "T"}
        rows.append({
            "raw_sample_id": s,
            "normalized_sample_id": re.sub(r"^CRC[-_](\d+)([A-Z])[-_](\d+)$", r"CRC-\1\2-\3", s),
            "patient_id": f"CRC_{patient:02d}" if isinstance(patient, int) else "UNRESOLVED",
            "tissue_group": tissue,
            "paired_set": f"CRC_{patient:02d}_NAT" if primary_nat and isinstance(patient, int) else "",
            "primary_NAT_candidate": "YES" if primary_nat else "NO",
            "adenoma_exploratory": "YES" if tissue == "D" else "NO",
            "CRC20_status": "CRC20" if patient == 20 else "not_CRC20",
            "source_file": str(COUNT_PATH),
            "notes": note,
        })
    meta = pd.DataFrame(rows)
    meta.to_csv(PHASE / "metadata/CRC_77_sample_metadata.tsv", sep="\t", index=False)
    all_samples = meta["raw_sample_id"].tolist()
    nat75 = meta.loc[meta["tissue_group"].isin(["N", "A", "T"]), "raw_sample_id"].tolist()
    crc20 = set(meta.loc[meta["CRC20_status"] == "CRC20", "raw_sample_id"])
    primary72 = [s for s in nat75 if s not in crc20]
    write_list(PHASE / "metadata/cohort_all_77_samples.txt", all_samples)
    write_list(PHASE / "metadata/cohort_sensitivity_NAT_75_samples.txt", nat75)
    # The decision can later mark this pending, but the list is useful when evidence confirms exclusion.
    write_list(PHASE / "metadata/cohort_primary_NAT_72_samples.txt", primary72)
    return meta, all_samples, nat75, primary72


def write_list(path, items):
    path.write_text("\n".join(items) + "\n")


def robust_z(x):
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad == 0 or np.isnan(mad):
        sd = np.nanstd(x)
        if sd == 0 or np.isnan(sd):
            return np.zeros_like(x)
        return (x - np.nanmean(x)) / sd
    return 0.6745 * (x - med) / mad


def median_ratio_norm(counts):
    arr = counts.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_arr = np.where(arr > 0, np.log(arr), np.nan)
    gm = np.exp(np.nanmean(log_arr, axis=1))
    valid = np.isfinite(gm) & (gm > 0)
    ratios = arr[valid, :] / gm[valid, None]
    ratios[ratios <= 0] = np.nan
    sf = np.nanmedian(ratios, axis=0)
    if np.any(~np.isfinite(sf)) or np.any(sf <= 0):
        lib = arr.sum(axis=0)
        sf = lib / np.nanmedian(lib)
    return arr / sf, sf


def pca_svd(mat, ncomp=10):
    x = np.asarray(mat, dtype=float)
    x = x - x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, ddof=1)
    std[std == 0] = 1
    x = x / std
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :ncomp] * s[:ncomp]
    var = (s ** 2) / max(1, x.shape[0] - 1)
    explained = var / var.sum()
    return scores, explained[:ncomp]


def compute_qc(count_df, meta, sample_cols, nat75, primary72):
    counts = count_df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    ann = count_df.drop(columns=sample_cols)
    if "sqanti_isoform_id" in count_df.columns:
        counts.index = count_df["sqanti_isoform_id"].astype(str)
    elif "isoform_id" in count_df.columns:
        counts.index = count_df["isoform_id"].astype(str)
    elif "merge_id" in count_df.columns:
        counts.index = count_df["merge_id"].astype(str)
    else:
        counts.index = count_df.iloc[:, 0].astype(str)
    lib = counts.sum(axis=0)
    nonzero = (counts > 0).sum(axis=0)
    med_nonzero = counts.mask(counts <= 0).median(axis=0).fillna(0)
    zero_frac = 1 - nonzero / counts.shape[0]
    gene_count = pd.Series(index=sample_cols, dtype=float)
    if "associated_gene" in count_df.columns:
        gene_arr = count_df["associated_gene"].astype(str)
        for s in sample_cols:
            gene_count[s] = gene_arr[counts[s].values > 0].replace({"nan": np.nan}).dropna().nunique()
    qc = pd.DataFrame({
        "sample_id": sample_cols,
        "library_size": lib.loc[sample_cols].values,
        "nonzero_isoform_count": nonzero.loc[sample_cols].values,
        "detected_gene_count": gene_count.loc[sample_cols].values,
        "median_nonzero_count": med_nonzero.loc[sample_cols].values,
        "zero_fraction": zero_frac.loc[sample_cols].values,
    }).merge(meta, left_on="sample_id", right_on="raw_sample_id", how="left")
    for c in ["library_size", "nonzero_isoform_count", "detected_gene_count", "median_nonzero_count", "zero_fraction"]:
        qc[c + "_robust_z_NAT75"] = np.nan
        mask = qc["sample_id"].isin(nat75)
        qc.loc[mask, c + "_robust_z_NAT75"] = robust_z(qc.loc[mask, c].values)
    qc.to_csv(PHASE / "tables/sample_expression_basic_qc.tsv", sep="\t", index=False)

    nat_counts = counts[nat75]
    norm, sf = median_ratio_norm(nat_counts.values)
    lognorm = np.log2(norm + 1)
    variances = np.var(lognorm, axis=1)
    top_n = min(10000, lognorm.shape[0])
    top_idx = np.argsort(variances)[-top_n:]
    data_samples_by_features = lognorm[top_idx, :].T
    scores, explained = pca_svd(data_samples_by_features, ncomp=min(10, len(nat75)))
    pca = pd.DataFrame(scores[:, :5], columns=[f"PC{i}" for i in range(1, 6)])
    pca["sample_id"] = nat75
    pca["size_factor"] = sf
    pca["PC1_var_explained"] = explained[0]
    pca["PC2_var_explained"] = explained[1] if len(explained) > 1 else np.nan
    pca = pca.merge(meta, left_on="sample_id", right_on="raw_sample_id", how="left")
    pca.to_csv(PHASE / "tables/PCA_NAT75_coordinates.tsv", sep="\t", index=False)
    pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(explained))], "variance_explained": explained}).to_csv(
        PHASE / "tables/PCA_NAT75_variance_explained.tsv", sep="\t", index=False
    )

    # Group center distance in PC1-PC5.
    pc_cols = [c for c in ["PC1", "PC2", "PC3", "PC4", "PC5"] if c in pca.columns]
    centers = pca.groupby("tissue_group")[pc_cols].mean()
    dist_rows = []
    for _, row in pca.iterrows():
        center = centers.loc[row["tissue_group"]].values
        vec = row[pc_cols].values.astype(float)
        dist_rows.append(np.linalg.norm(vec - center))
    pca["distance_to_tissue_center_PC1_5"] = dist_rows
    pca["distance_to_tissue_center_robust_z"] = pca.groupby("tissue_group")["distance_to_tissue_center_PC1_5"].transform(lambda x: robust_z(x.values))
    pca.to_csv(PHASE / "tables/PCA_NAT75_coordinates_with_distances.tsv", sep="\t", index=False)

    corr = np.corrcoef(data_samples_by_features)
    corr_df = pd.DataFrame(corr, index=nat75, columns=nat75)
    corr_df.to_csv(PHASE / "tables/sample_correlation_NAT75.tsv", sep="\t")
    nearest_rows = []
    for i, s in enumerate(nat75):
        vals = corr[i].copy()
        vals[i] = -np.inf
        j = int(np.argmax(vals))
        nearest_rows.append({
            "sample_id": s,
            "nearest_neighbor": nat75[j],
            "nearest_neighbor_correlation": corr[i, j],
            "nearest_neighbor_tissue": meta.set_index("raw_sample_id").loc[nat75[j], "tissue_group"],
            "median_correlation_all_NAT75": np.median(corr[i][np.arange(len(nat75)) != i]),
            "median_correlation_same_tissue": np.median([corr[i, k] for k, x in enumerate(nat75)
                                                         if k != i and meta.set_index("raw_sample_id").loc[x, "tissue_group"] == meta.set_index("raw_sample_id").loc[s, "tissue_group"]]),
        })
    nearest = pd.DataFrame(nearest_rows)
    nearest["median_correlation_all_robust_z"] = robust_z(nearest["median_correlation_all_NAT75"].values)
    nearest.to_csv(PHASE / "tables/sample_nearest_neighbors_and_correlation_summary.tsv", sep="\t", index=False)

    patient_dist = []
    pca_idx = pca.set_index("sample_id")
    meta_nat = meta[meta["raw_sample_id"].isin(nat75)]
    for patient, g in meta_nat.groupby("patient_id"):
        ids = {r["tissue_group"]: r["raw_sample_id"] for _, r in g.iterrows()}
        if all(t in ids for t in ["N", "A", "T"]):
            pairs = [("N", "A"), ("A", "T"), ("N", "T")]
            vals = []
            for a, b in pairs:
                va = pca_idx.loc[ids[a], pc_cols].values.astype(float)
                vb = pca_idx.loc[ids[b], pc_cols].values.astype(float)
                d = float(np.linalg.norm(va - vb))
                vals.append(d)
                patient_dist.append({"patient_id": patient, "pair": f"{a}-{b}", "distance_PC1_5": d})
            patient_dist.append({"patient_id": patient, "pair": "mean_pair_distance", "distance_PC1_5": float(np.mean(vals))})
    pdist = pd.DataFrame(patient_dist)
    if not pdist.empty:
        m = pdist["pair"] == "mean_pair_distance"
        pdist.loc[m, "mean_pair_distance_robust_z"] = robust_z(pdist.loc[m, "distance_PC1_5"].values)
    pdist.to_csv(PHASE / "tables/patient_within_pair_distance.tsv", sep="\t", index=False)

    # Class composition by sample based on expressed isoforms.
    if "structural_category" in count_df.columns:
        comp_rows = []
        cats = count_df["structural_category"].astype(str)
        for s in nat75:
            expressed = counts[s].values > 0
            c = Counter(cats[expressed])
            total = sum(c.values())
            for cat, val in c.items():
                comp_rows.append({"sample_id": s, "structural_category": cat, "expressed_isoforms": val, "fraction": val / total if total else 0})
        comp = pd.DataFrame(comp_rows).merge(meta, left_on="sample_id", right_on="raw_sample_id", how="left")
        comp.to_csv(PHASE / "tables/SQANTI3_class_composition_by_sample.tsv", sep="\t", index=False)

    return counts, ann, qc, pca, corr_df, nearest, pdist


def filter_impact(counts, meta, sample_cols, nat75, primary72, category_map=None):
    all77 = sample_cols
    def calc(cols):
        sub = counts[cols]
        support = (sub >= 2).sum(axis=1)
        total = sub.sum(axis=1)
        passed = (support >= 8) & (total >= 20)
        return support, total, passed
    s77, t77, p77 = calc(all77)
    s75, t75, p75 = calc(nat75)
    s72, t72, p72 = calc(primary72)
    out = pd.DataFrame({
        "isoform_id": counts.index,
        "support_ge2_all77": s77.values,
        "total_count_all77": t77.values,
        "pass_all77": p77.values,
        "support_ge2_NAT75": s75.values,
        "total_count_NAT75": t75.values,
        "pass_NAT75": p75.values,
        "support_ge2_primary72": s72.values,
        "total_count_primary72": t72.values,
        "pass_primary72": p72.values,
    })
    if category_map is not None:
        out.insert(1, "structural_category", [category_map.get(x, "NA") for x in out["isoform_id"]])
    out["pass_change_due_to_D"] = out["pass_all77"] & ~out["pass_NAT75"]
    out["pass_change_due_to_CRC20"] = out["pass_NAT75"] & ~out["pass_primary72"]
    out.to_csv(PHASE / "tables/isoform_filter_cohort_impact.tsv", sep="\t", index=False)
    if "structural_category" in out.columns:
        cat_summary = out.groupby("structural_category", dropna=False).agg(
            n_isoforms=("isoform_id", "count"),
            pass_all77=("pass_all77", "sum"),
            pass_NAT75=("pass_NAT75", "sum"),
            pass_primary72=("pass_primary72", "sum"),
            lost_after_removing_D=("pass_change_due_to_D", "sum"),
            lost_after_removing_CRC20=("pass_change_due_to_CRC20", "sum"),
        ).reset_index()
        cat_summary.to_csv(PHASE / "tables/isoform_filter_cohort_impact_by_sqanti3_category.tsv", sep="\t", index=False)
    summary = {
        "n_isoforms": int(len(out)),
        "pass_all77": int(out["pass_all77"].sum()),
        "pass_NAT75": int(out["pass_NAT75"].sum()),
        "pass_primary72": int(out["pass_primary72"].sum()),
        "lost_after_removing_D": int(out["pass_change_due_to_D"].sum()),
        "lost_after_removing_CRC20": int(out["pass_change_due_to_CRC20"].sum()),
    }
    pd.DataFrame([summary]).to_csv(PHASE / "tables/isoform_filter_cohort_impact_summary.tsv", sep="\t", index=False)
    return out, summary


def id_concordance(count_df):
    if "sqanti_isoform_id" in count_df.columns:
        count_ids = set(count_df["sqanti_isoform_id"].astype(str))
    elif "isoform_id" in count_df.columns:
        count_ids = set(count_df["isoform_id"].astype(str))
    elif "merge_id" in count_df.columns:
        count_ids = set(count_df["merge_id"].astype(str))
    else:
        count_ids = set(count_df.iloc[:, 0].astype(str))
    orig_ids = set(count_df["original_isoform_id"].astype(str)) if "original_isoform_id" in count_df.columns else set()
    class_ids = set()
    class_counts = Counter()
    if CLASS_PATH.exists():
        cls = pd.read_csv(CLASS_PATH, sep="\t", usecols=["isoform", "structural_category"])
        class_ids = set(cls["isoform"].astype(str))
        class_counts = Counter(cls["structural_category"].astype(str))
    gtf_ids = set()
    if GTF_PATH.exists():
        pat = re.compile(r'transcript_id "([^"]+)"')
        with GTF_PATH.open(errors="replace") as f:
            for line in f:
                if "\ttranscript\t" in line or "transcript_id" in line:
                    m = pat.search(line)
                    if m:
                        gtf_ids.add(m.group(1))
    fa_ids = set()
    if FASTA_PATH.exists():
        with FASTA_PATH.open(errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    fa_ids.add(line[1:].strip().split()[0])
    rows = []
    count_label = "count_sqanti_isoform_id" if "sqanti_isoform_id" in count_df.columns else ("count_isoform_id" if "isoform_id" in count_df.columns else "count_merge_id")
    for name, ids in [(count_label, count_ids), ("count_original_isoform_id", orig_ids), ("classification_isoform", class_ids), ("gtf_transcript_id", gtf_ids), ("fasta_header_id", fa_ids)]:
        rows.append({"id_source": name, "n_ids": len(ids)})
    base = count_ids
    for name, ids in [("classification_isoform", class_ids), ("gtf_transcript_id", gtf_ids), ("fasta_header_id", fa_ids)]:
        rows.append({"id_source": f"intersection_{count_label}_vs_{name}", "n_ids": len(base & ids)})
        rows.append({"id_source": f"{count_label}_missing_from_{name}", "n_ids": len(base - ids)})
        if name == "classification_isoform":
            rows.append({"id_source": f"{name}_missing_from_{count_label}", "n_ids": len(ids - base)})
    pd.DataFrame(rows).to_csv(PHASE / "tables/isoform_id_concordance.tsv", sep="\t", index=False)
    pd.DataFrame([{"structural_category": k, "n_isoforms": v} for k, v in class_counts.most_common()]).to_csv(
        PHASE / "tables/final_sqanti3_structural_category_counts.tsv", sep="\t", index=False
    )
    return rows, class_counts


def decide_crc20(qc, pca, nearest, pdist):
    evidence = []
    against = []
    crc20_samples = set(qc.loc[qc["CRC20_status"] == "CRC20", "sample_id"])
    qc20 = qc[qc["sample_id"].isin(crc20_samples)]
    pca20 = pca[pca["sample_id"].isin(crc20_samples)]
    near20 = nearest[nearest["sample_id"].isin(crc20_samples)]
    pd20 = pdist[(pdist["patient_id"] == "CRC_20") & (pdist["pair"] == "mean_pair_distance")]
    dims = []
    for c in ["library_size", "nonzero_isoform_count", "detected_gene_count", "median_nonzero_count", "zero_fraction"]:
        zc = c + "_robust_z_NAT75"
        if zc in qc20:
            mx = float(np.nanmax(np.abs(qc20[zc].values)))
            if mx >= 3:
                dims.append(c)
                evidence.append(f"TECHNICAL_OUTLIER_EVIDENCE: CRC_20 has |robust z| >= 3 for {c} (max {mx:.2f}).")
            else:
                against.append(f"NO_CLEAR_OUTLIER_EVIDENCE: CRC_20 {c} max |robust z| = {mx:.2f}.")
    if "distance_to_tissue_center_robust_z" in pca20:
        mx = float(np.nanmax(np.abs(pca20["distance_to_tissue_center_robust_z"].values)))
        if mx >= 3:
            dims.append("PCA_tissue_center_distance")
            evidence.append(f"TECHNICAL_OUTLIER_EVIDENCE: CRC_20 is distant from same-tissue PCA centers (max |robust z| {mx:.2f}).")
        else:
            against.append(f"NO_CLEAR_OUTLIER_EVIDENCE: CRC_20 same-tissue PCA center distance max |robust z| = {mx:.2f}.")
    if not near20.empty:
        mx = float(np.nanmax(np.abs(near20["median_correlation_all_robust_z"].values)))
        if mx >= 3:
            dims.append("sample_correlation")
            evidence.append(f"TECHNICAL_OUTLIER_EVIDENCE: CRC_20 median sample correlation is abnormal (max |robust z| {mx:.2f}).")
        else:
            against.append(f"NO_CLEAR_OUTLIER_EVIDENCE: CRC_20 median correlation max |robust z| = {mx:.2f}.")
    if not pd20.empty and "mean_pair_distance_robust_z" in pd20:
        z = float(pd20["mean_pair_distance_robust_z"].iloc[0])
        if abs(z) >= 3:
            dims.append("within_patient_pair_distance")
            evidence.append(f"LABEL_OR_PAIRING_EVIDENCE: CRC_20 has abnormal within-patient N/A/T distance (robust z {z:.2f}).")
        else:
            against.append(f"NO_CLEAR_OUTLIER_EVIDENCE: CRC_20 within-patient mean distance robust z = {z:.2f}.")
    # Label evidence from nearest neighbors: if all nearest neighbors are unexpected tissue.
    label_flags = []
    for _, r in near20.iterrows():
        own_t = qc.set_index("sample_id").loc[r["sample_id"], "tissue_group"]
        if r["nearest_neighbor_tissue"] != own_t:
            label_flags.append(f"{r['sample_id']} nearest neighbor tissue={r['nearest_neighbor_tissue']} not {own_t}")
    if len(label_flags) >= 2:
        dims.append("nearest_neighbor_tissue_mismatch")
        evidence.append("LABEL_OR_PAIRING_EVIDENCE: " + "; ".join(label_flags))
    elif label_flags:
        against.append("Limited nearest-neighbor tissue mismatch: " + "; ".join(label_flags))
    if len(set(dims)) >= 2:
        decision = "CONFIRMED_EXCLUDE_WHOLE_PATIENT_FROM_PRIMARY"
    else:
        decision = "PENDING_MANUAL_REVIEW"
    return decision, evidence, against, dims


def write_reports(meta, audit, id_rows, class_counts, qc, pca, nearest, pdist, filt_summary, decision, evidence, against, dims):
    nat = meta[meta["tissue_group"].isin(["N", "A", "T"])]
    d = meta[meta["tissue_group"] == "D"]
    complete = []
    for patient, g in nat.groupby("patient_id"):
        complete.append({"patient_id": patient, "tissues": ",".join(sorted(g["tissue_group"])), "n": len(g),
                         "complete_NAT": set(g["tissue_group"]) == {"N", "A", "T"}})
    complete_df = pd.DataFrame(complete)
    complete_df.to_csv(PHASE / "tables/patient_pair_completeness.tsv", sep="\t", index=False)
    dup = meta["raw_sample_id"][meta["raw_sample_id"].duplicated()].tolist()
    unresolved = meta.loc[meta["tissue_group"] == "UNRESOLVED", "raw_sample_id"].tolist()
    crc20 = meta.loc[meta["CRC20_status"] == "CRC20", "raw_sample_id"].tolist()

    sample_audit = f"""# Phase 00 sample audit

## Inputs

- Count matrix: `{COUNT_PATH}`
- TPM matrix: `{TPM_PATH}`
- Strict SQANTI3 directory: `{SQANTI_DIR}`
- Strict filtered transcriptome input: `{STRICT_INPUT}`

## Sample composition

- Total samples in count matrix: {len(meta)}
- N/A/T samples: {len(nat)}
- D adenoma exploratory samples: {len(d)}
- Complete N/A/T patients: {int(complete_df['complete_NAT'].sum())}
- D samples: {', '.join(d['raw_sample_id'])}
- CRC_20 samples present: {', '.join(crc20)}
- Duplicate sample names: {', '.join(dup) if dup else 'None'}
- Unresolved sample names: {', '.join(unresolved) if unresolved else 'None'}

## Notes

Existing project metadata maps terminal tissue code A/B/C to Tumor/Adjacent/Normal. Therefore standardized tissue groups are T/A/N, while terminal D is marked as adenoma exploratory and excluded from N/A/T main analysis.

Input file audit table: `tables/input_file_audit.tsv`
Isoform ID concordance table: `tables/isoform_id_concordance.tsv`
"""
    (PHASE / "reports/phase00_sample_audit.md").write_text(sample_audit)

    ev = "\n".join(f"- {x}" for x in evidence) if evidence else "- No strong outlier evidence reached the predefined thresholds."
    ag = "\n".join(f"- {x}" for x in against) if against else "- No major counter-evidence recorded."
    crc20_report = f"""# CRC_20 outlier assessment

## Research background

The CRC paired long-read transcriptome study contains 77 samples: 75 N/A/T samples from 25 CRC patients plus 2 D adenoma exploratory samples. CRC_20 had previously been suspected as an outlier. This Phase 00 assessment rechecks CRC_20 using objective QC evidence before locking the main cohort.

## Data used

- Strict isoform count matrix with SQANTI3 annotations: `{COUNT_PATH}`
- Strict SQANTI3 classification: `{CLASS_PATH}`
- N/A/T samples only for outlier modeling: 75 samples
- D samples excluded from PCA/correlation/outlier modeling.

## QC dimensions used

1. Expression matrix basic QC: library size, nonzero isoform count, detected gene count, median nonzero count, zero fraction.
2. Size-factor normalized log2(count + 1) PCA using the top variable isoforms.
3. Sample correlation and nearest-neighbor tissue consistency.
4. Within-patient N/A/T distance in PC space.
5. SQANTI3 class composition by expressed isoforms when count > 0.

## Evidence supporting exclusion or concern

{ev}

## Evidence against exclusion

{ag}

## Technical abnormality vs biological extreme

Evidence categories triggered: {', '.join(sorted(set(dims))) if dims else 'none above threshold'}.

This assessment does not use downstream biological significance as an exclusion criterion. If CRC_20 is retained only in sensitivity analysis, that decision is based on pre-specified QC and paired-design considerations, not on making differential results stronger.

## Final determination

`{decision}`

## Main analysis recommendation

{"Exclude the whole CRC_20 patient from the primary paired N/A/T analysis and keep all three CRC_20 N/A/T samples in the 75-sample sensitivity cohort." if decision == "CONFIRMED_EXCLUDE_WHOLE_PATIENT_FROM_PRIMARY" else "Do not treat the 72-sample list as final without manual review; keep the 75-sample N/A/T sensitivity cohort."}

## Sensitivity analysis recommendation

Always retain a 75-sample N/A/T sensitivity cohort including CRC_20. The two D samples remain excluded from N/A/T statistics and are reserved for separate exploratory adenoma analyses.

## Limitations

Read length, mapping rate, and per-sample sequencing-quality metrics were only assessed if present in existing metadata/logs. Missing technical metrics were not imputed.
"""
    (PHASE / "reports/CRC20_outlier_assessment.md").write_text(crc20_report)

    cohort_lock = f"""# Cohort lock report

## Cohort composition

- Total samples: {len(meta)}
- N/A/T samples: {len(nat)}
- D exploratory adenoma samples: {len(d)} ({', '.join(d['raw_sample_id'])})
- Complete N/A/T patient triplets: {int(complete_df['complete_NAT'].sum())}

## CRC_20 handling

Final Phase 00 determination: `{decision}`.

If exclusion is used for the main analysis, the whole CRC_20 patient triplet is excluded to preserve complete paired design.

## Locked cohorts

- All samples: `metadata/cohort_all_77_samples.txt`
- Sensitivity N/A/T cohort: `metadata/cohort_sensitivity_NAT_75_samples.txt`
- Primary N/A/T 72-sample cohort: `metadata/cohort_primary_NAT_72_samples.txt` if the CRC_20 exclusion determination is accepted.

Primary analysis sample size under exclusion: 24 patients / 72 N/A/T samples.
Sensitivity analysis sample size: 25 patients / 75 N/A/T samples.

## Downstream paired statistical design

Use patient-blocked or paired models for N vs A vs T comparisons. D samples must not enter formal N/A/T paired testing.

## Strict filter impact

- Isoforms passing all77 threshold: {filt_summary['pass_all77']}
- Passing NAT75 after excluding D samples: {filt_summary['pass_NAT75']}
- Passing primary72 after excluding D and CRC_20: {filt_summary['pass_primary72']}
- Lost after removing D samples: {filt_summary['lost_after_removing_D']}
- Lost after removing CRC_20: {filt_summary['lost_after_removing_CRC20']}

Recommendation: do not rerun FLAIR collapse or SQANTI3 in Phase 00. In formal expression analysis, re-apply the support filter at the analysis layer for the locked cohort and report how many isoforms are affected.

Category-level impact table: `tables/isoform_filter_cohort_impact_by_sqanti3_category.tsv`
"""
    (PHASE / "reports/cohort_lock_report.md").write_text(cohort_lock)

    methods = f"""# CRC_20 exclusion methods text

## English draft

Before differential isoform analyses, we performed a pre-specified cohort-locking quality-control assessment using the strict long-read isoform count matrix annotated by SQANTI3. The two adenoma samples were retained in the full sample inventory but excluded from the paired normal-adjacent-tumor analysis. CRC_20 was re-evaluated using multiple independent QC dimensions, including library size, detected isoform and gene counts, robust outlier scores, size-factor-normalized PCA, sample-sample correlation, nearest-neighbor tissue consistency, and within-patient paired distances. When CRC_20 showed reproducible outlier evidence across independent QC dimensions, all three samples from this patient were excluded from the primary paired analysis to preserve a complete paired design. The full 75-sample N/A/T cohort including CRC_20 was retained for sensitivity analyses. This exclusion rule was defined before downstream differential testing and was not based on improving statistical significance.

## 中文说明

正式差异 isoform 分析前，先基于严格过滤后的长读长 isoform count matrix 和 SQANTI3 注释进行队列锁定 QC。两个 D 样本保留在总样本清单中，但不进入 N/A/T 配对主分析。CRC_20 使用多个独立维度复核，包括测序深度、检测到的 isoform/gene 数、稳健异常值分数、归一化 PCA、样本相关性、最近邻组织一致性和患者内配对距离。若 CRC_20 在多个独立 QC 维度中稳定异常，则从主分析中排除该患者的 N/A/T 三个样本，以维持完整配对设计；同时保留包含 CRC_20 的 75 个 N/A/T 样本作为敏感性分析队列。该处理基于预设 QC 和配对设计原则，不基于差异分析结果显著性。
"""
    (PHASE / "reports/CRC20_exclusion_methods_text.md").write_text(methods)

    worklog = f"""# Phase 00 worklog

## Start/end

- Start: {dt.datetime.now().isoformat(timespec='seconds')}
- End: generated at script completion

## Inputs

- `{COUNT_PATH}`
- `{TPM_PATH}`
- `{CLASS_PATH}`
- `{GTF_PATH}`
- `{FASTA_PATH}`
- `{META_HINT}`

## Commands/scripts

- Main script: `code/phase00_crc_cohort_lock.py`
- Figure script: `code/phase00_make_figures.R`
- R version is recorded in `provenance/session_info_R.txt` after figure generation.
- Python version: {sys.version.split()[0]}
- pandas: {pd.__version__}
- numpy: {np.__version__}

## Parameters

- Support threshold: count >= 2 in at least 8 samples and total count >= 20.
- PCA input: size-factor normalized log2(count + 1), top 10,000 most variable isoforms among N/A/T samples.
- Robust z-score: median/MAD scaled z.

## Outputs

- Metadata: `metadata/CRC_77_sample_metadata.tsv`
- Cohorts: `metadata/cohort_all_77_samples.txt`, `metadata/cohort_sensitivity_NAT_75_samples.txt`, `metadata/cohort_primary_NAT_72_samples.txt`
- Tables: see `tables/`
- Figures: generated by R into `figures/`
- Reports: this folder.

## Main results

- CRC_20 decision: `{decision}`
- Filter impact summary: {json.dumps(filt_summary, ensure_ascii=False)}

## Errors and fixes

- No upstream directories were modified.
- If optional technical QC metrics were unavailable, they were recorded as a limitation rather than imputed.

## Unresolved issues

- Manual review is still recommended for any label-swap suspicion raised by nearest-neighbor or clustering outputs.

## Next step

Proceed to Phase 01 only after reviewing this cohort lock report and accepting the primary/sensitivity cohort definitions.
"""
    (PHASE / "reports/phase00_worklog.md").write_text(worklog)


def main():
    mkdirs()
    audit = limited_audit()
    raw_count_df = pd.read_csv(COUNT_PATH, sep="\t")
    cls_for_main = pd.read_csv(CLASS_PATH, sep="\t")
    keep_ids = set(cls_for_main["isoform"].astype(str))
    if "isoform_id" not in raw_count_df.columns:
        raise RuntimeError(f"Expected isoform_id column in strict count matrix: {COUNT_PATH}")
    raw_count_df["sqanti_isoform_id"] = raw_count_df["isoform_id"].astype(str).str.replace(
        r"_ENSG[0-9A-Za-z_.-]+$", "", regex=True
    )
    count_df = raw_count_df[raw_count_df["sqanti_isoform_id"].astype(str).isin(keep_ids)].copy()
    # Add SQANTI3 annotation columns for QC summaries while keeping the strict count matrix as the expression source.
    ann_cols = [c for c in ["isoform", "structural_category", "associated_gene", "associated_transcript", "subcategory"] if c in cls_for_main.columns]
    count_df = count_df.merge(cls_for_main[ann_cols], left_on="sqanti_isoform_id", right_on="isoform", how="left")
    sample_cols = detect_sample_cols(count_df.columns)
    meta, all77, nat75, primary72 = build_metadata(sample_cols)
    counts, ann, qc, pca, corr, nearest, pdist = compute_qc(count_df, meta, sample_cols, nat75, primary72)
    category_map = {}
    if "sqanti_isoform_id" in count_df.columns and "structural_category" in count_df.columns:
        category_map = dict(zip(count_df["sqanti_isoform_id"].astype(str), count_df["structural_category"].astype(str)))
    filt, filt_summary = filter_impact(counts, meta, sample_cols, nat75, primary72, category_map)
    id_rows, class_counts = id_concordance(count_df)
    decision, evidence, against, dims = decide_crc20(qc, pca, nearest, pdist)
    if decision != "CONFIRMED_EXCLUDE_WHOLE_PATIENT_FROM_PRIMARY":
        pending = PHASE / "metadata/cohort_primary_NAT_PENDING.txt"
        pending.write_text("CRC_20 exclusion not confirmed automatically; manual review required before locking primary cohort.\n")
    write_reports(meta, audit, id_rows, class_counts, qc, pca, nearest, pdist, filt_summary, decision, evidence, against, dims)
    provenance = {
        "timestamp": STAMP,
        "root": str(ROOT),
        "out_root": str(OUT_ROOT),
        "count_path": str(COUNT_PATH),
        "tpm_path": str(TPM_PATH),
        "class_path": str(CLASS_PATH),
        "gtf_path": str(GTF_PATH),
        "fasta_path": str(FASTA_PATH),
        "sample_columns": sample_cols,
        "decision": decision,
    }
    (PHASE / "provenance/phase00_provenance.json").write_text(json.dumps(provenance, indent=2))
    print(f"PHASE00_COMPLETE\t{PHASE}\t{decision}")


if __name__ == "__main__":
    main()
