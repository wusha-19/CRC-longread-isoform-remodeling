#!/usr/bin/env python3
"""Build publication-ready supplementary three-line tables from frozen Fig. 3 robustness outputs."""

from __future__ import annotations

import csv
import hashlib
import math
import os
from pathlib import Path
from typing import Iterable

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


SEED = 20260827
ANALYSIS_ROOT = Path(
    os.environ["CRC_FIG3_ROBUSTNESS_ROOT"]
).expanduser().resolve()
RESULTS = ANALYSIS_ROOT / "results"
OUT = ANALYSIS_ROOT / "supplementary_tables" / "concise"
SOURCE = OUT / "source_data"
DOCX_PATH = OUT / "Supplementary_Table_Fig3_measurement_robustness_concise.docx"
REPORT_PATH = OUT / "SUPPLEMENTARY_TABLE_PACKAGE_REPORT_zh.md"
MANIFEST_PATH = OUT / "Supplementary_Table_source_manifest.tsv"


def read_tsv(name: str) -> list[dict[str, str]]:
    path = RESULTS / name
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing frozen result: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def f(value: str | float | int | None) -> float | None:
    if value is None or value == "" or str(value).upper() == "NA":
        return None
    return float(value)


def fmt(value: float | None, digits: int = 3) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def fmt_ci(low: float | None, high: float | None) -> str:
    if low is None or high is None:
        return "NA"
    return f"{fmt(low)} to {fmt(high)}"


def fmt_p(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_table_a() -> list[dict[str, object]]:
    lr = {r["contrast"]: r for r in read_tsv("Fig3_LR_TVD_reliability_summary.tsv")}
    sr = {r["contrast"]: r for r in read_tsv("Fig3_SR_log2FC_reliability_summary.tsv")}
    adj = {r["contrast"]: r for r in read_tsv("Fig3_reliability_adjusted_coupling.tsv")}
    partial = {
        r["contrast"]: r
        for r in read_tsv("Fig3_partial_correlation.tsv")
        if r["model"] == "primary"
    }
    rows: list[dict[str, object]] = []
    for contrast in ("T_vs_A", "T_vs_N"):
        lrow, srow, arow, prow = lr[contrast], sr[contrast], adj[contrast], partial[contrast]
        rows.extend(
            [
                {
                    "contrast": contrast,
                    "measure": "Long-read TVD pseudo-replicate reliability",
                    "n": int(lrow["n_replicates"]),
                    "n_unit": "pseudo-replicates",
                    "estimate": f(lrow["median"]),
                    "IQR_low": f(lrow["Q1"]),
                    "IQR_high": f(lrow["Q3"]),
                    "interval_low": f(lrow["percentile_2.5"]),
                    "interval_high": f(lrow["percentile_97.5"]),
                    "interval_type": "2.5th-97.5th percentile",
                    "P": None,
                    "definition": lrow["replicate_definition"],
                },
                {
                    "contrast": contrast,
                    "measure": "Matched short-read |log2FC| pseudo-library reliability",
                    "n": int(srow["n_replicates"]),
                    "n_unit": "pseudo-libraries",
                    "estimate": f(srow["median"]),
                    "IQR_low": f(srow["Q1"]),
                    "IQR_high": f(srow["Q3"]),
                    "interval_low": f(srow["percentile_2.5"]),
                    "interval_high": f(srow["percentile_97.5"]),
                    "interval_type": "2.5th-97.5th percentile",
                    "P": None,
                    "definition": srow["replicate_definition"],
                },
                {
                    "contrast": contrast,
                    "measure": "Observed cross-layer coupling",
                    "n": int(arow["n_genes"]),
                    "n_unit": "genes",
                    "estimate": f(arow["rho_observed"]),
                    "IQR_low": None,
                    "IQR_high": None,
                    "interval_low": f(arow["observed_CI_low"]),
                    "interval_high": f(arow["observed_CI_high"]),
                    "interval_type": "bootstrap 95% CI",
                    "P": f(arow["observed_P"]),
                    "definition": "Spearman correlation between matched short-read |log2FC| and long-read gene-level TVD.",
                },
                {
                    "contrast": contrast,
                    "measure": "Approximate reliability-adjusted coupling",
                    "n": int(arow["n_genes"]),
                    "n_unit": "genes",
                    "estimate": f(arow["rho_reliability_adjusted"]),
                    "IQR_low": None,
                    "IQR_high": None,
                    "interval_low": f(arow["adjusted_CI_low"]),
                    "interval_high": f(arow["adjusted_CI_high"]),
                    "interval_type": "bootstrap 95% CI",
                    "P": None,
                    "definition": arow["interpretation_boundary"],
                },
                {
                    "contrast": contrast,
                    "measure": "Primary covariate-adjusted coupling",
                    "n": int(prow["n_genes"]),
                    "n_unit": "genes",
                    "estimate": f(prow["partial_rho"]),
                    "IQR_low": None,
                    "IQR_high": None,
                    "interval_low": f(prow["CI_low"]),
                    "interval_high": f(prow["CI_high"]),
                    "interval_type": "bootstrap 95% CI",
                    "P": f(prow["P"]),
                    "definition": "Rank-residual partial correlation adjusted for long-read count, short-read baseMean, isoform multiplicity, and zero fraction.",
                },
            ]
        )
    return rows


def build_table_b() -> list[dict[str, object]]:
    files = [
        ("Long-read abundance", "Fig3_coupling_by_LR_abundance.tsv"),
        ("Short-read abundance", "Fig3_coupling_by_SR_abundance.tsv"),
        ("Isoform multiplicity", "Fig3_coupling_by_isoform_multiplicity.tsv"),
    ]
    rows: list[dict[str, object]] = []
    for label, name in files:
        for row in read_tsv(name):
            rows.append(
                {
                    "stratification": label,
                    "contrast": row["contrast"],
                    "stratum": row["stratum"],
                    "n_genes": int(row["n_genes"]),
                    "rho": f(row["rho"]),
                    "CI_low": f(row["CI_low"]),
                    "CI_high": f(row["CI_high"]),
                    "P": f(row["P"]),
                    "bootstrap_replicates": int(row["bootstrap_replicates"]),
                    "stratum_definition": row["stratum_type"],
                }
            )
    return rows


def build_table_c() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv("Fig3_incremental_explanatory_value.tsv"):
        rows.append(
            {
                "contrast": row["contrast"],
                "model": row["model_set"],
                "n_genes": int(row["n_genes"]),
                "apparent_metric": row["apparent_metric"],
                "baseline_apparent": f(row["baseline_apparent"]),
                "extended_apparent": f(row["extended_apparent"]),
                "delta_apparent": f(row["delta_apparent"]),
                "added_abs_log2FC_P": f(row["added_abs_log2FC_P"]),
                "baseline_CV_median": f(row["baseline_CV_median"]),
                "baseline_CV_CI_low": f(row["baseline_CV_CI_low"]),
                "baseline_CV_CI_high": f(row["baseline_CV_CI_high"]),
                "extended_CV_median": f(row["extended_CV_median"]),
                "extended_CV_CI_low": f(row["extended_CV_CI_low"]),
                "extended_CV_CI_high": f(row["extended_CV_CI_high"]),
                "delta_CV_median": f(row["delta_CV_median"]),
                "delta_CV_CI_low": f(row["delta_CV_CI_low"]),
                "delta_CV_CI_high": f(row["delta_CV_CI_high"]),
                "CV_repeats": int(row["CV_repeats"]),
                "CV_folds": int(row["CV_folds"]),
            }
        )
    return rows


def build_table_d() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    high = read_tsv("Fig3_high_confidence_summary.tsv")
    for row in high:
        contrast = row["contrast"]
        n = int(row["n_high_confidence_genes"])
        rows.extend(
            [
                {
                    "analysis_set": "Pre-specified high-confidence subset",
                    "contrast": contrast,
                    "n_genes": n,
                    "metric": "Observed Spearman rho",
                    "estimate": f(row["observed_rho"]),
                    "CI_low": f(row["observed_CI_low"]),
                    "CI_high": f(row["observed_CI_high"]),
                    "P": f(row["observed_P"]),
                    "definition": "TVD outcome; subset defined before examining coupling estimates.",
                },
                {
                    "analysis_set": "Pre-specified high-confidence subset",
                    "contrast": contrast,
                    "n_genes": n,
                    "metric": "Primary partial rho",
                    "estimate": f(row["partial_rho"]),
                    "CI_low": f(row["partial_CI_low"]),
                    "CI_high": f(row["partial_CI_high"]),
                    "P": f(row["partial_P"]),
                    "definition": "Covariate-adjusted rank-residual correlation.",
                },
                {
                    "analysis_set": "Pre-specified high-confidence subset",
                    "contrast": contrast,
                    "n_genes": n,
                    "metric": "Incremental cross-validated explained variance",
                    "estimate": f(row["incremental_delta_CV_median"]),
                    "CI_low": None,
                    "CI_high": None,
                    "P": None,
                    "definition": "Gain after adding matched short-read |log2FC| to the covariate model.",
                },
            ]
        )
    max_rows = read_tsv("Fig3_maxdelta_sensitivity.tsv")
    allowed = {"raw correlation", "partial correlation", "incremental CV"}
    for row in max_rows:
        if row["analysis"] not in allowed:
            continue
        if row["analysis"] == "partial correlation" and row["stratum"] != "primary":
            continue
        if row["analysis"] == "incremental CV" and row["stratum"] != "rank_linear":
            continue
        metric = {
            "raw correlation": "Observed Spearman rho",
            "partial correlation": "Primary partial rho",
            "incremental CV": "Incremental cross-validated explained variance",
        }[row["analysis"]]
        rows.append(
            {
                "analysis_set": "Maximum-absolute-delta sensitivity",
                "contrast": row["contrast"],
                "n_genes": int(row["n_genes"]),
                "metric": metric,
                "estimate": f(row["estimate"]),
                "CI_low": f(row["CI_low"]),
                "CI_high": f(row["CI_high"]),
                "P": f(row["P"]),
                "definition": row["detail"],
            }
        )
    return rows


def build_concise_rows(
    table_a: list[dict[str, object]],
    table_c: list[dict[str, object]],
    table_d: list[dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    category_map = {
        "Long-read TVD pseudo-replicate reliability": "Measurement reliability",
        "Matched short-read |log2FC| pseudo-library reliability": "Measurement reliability",
        "Observed cross-layer coupling": "Primary coupling",
        "Approximate reliability-adjusted coupling": "Primary coupling",
        "Primary covariate-adjusted coupling": "Primary coupling",
    }
    for row in table_a:
        rows.append(
            {
                "category": category_map[str(row["measure"])],
                "contrast": row["contrast"],
                "metric": row["measure"],
                "n": row["n"],
                "n_unit": row["n_unit"],
                "estimate": row["estimate"],
                "IQR_low": row["IQR_low"],
                "IQR_high": row["IQR_high"],
                "interval_low": row["interval_low"],
                "interval_high": row["interval_high"],
                "interval_type": row["interval_type"],
                "P": row["P"],
                "key_note": row["definition"],
                "source_file": (
                    "Fig3_LR_TVD_reliability_summary.tsv"
                    if str(row["measure"]).startswith("Long-read")
                    else "Fig3_SR_log2FC_reliability_summary.tsv"
                    if str(row["measure"]).startswith("Matched short-read")
                    else "Fig3_reliability_adjusted_coupling.tsv"
                    if "reliability-adjusted" in str(row["measure"]) or "Observed" in str(row["measure"])
                    else "Fig3_partial_correlation.tsv"
                ),
            }
        )
    for row in table_c:
        if row["model"] != "rank_linear":
            continue
        rows.append(
            {
                "category": "Incremental explanatory value",
                "contrast": row["contrast"],
                "metric": "Incremental cross-validated explained variance",
                "n": row["n_genes"],
                "n_unit": "genes",
                "estimate": row["delta_CV_median"],
                "IQR_low": None,
                "IQR_high": None,
                "interval_low": row["delta_CV_CI_low"],
                "interval_high": row["delta_CV_CI_high"],
                "interval_type": "2.5th-97.5th percentile",
                "P": None,
                "key_note": (
                    f"Rank-linear primary model; baseline CV={fmt(row['baseline_CV_median'])}, "
                    f"extended CV={fmt(row['extended_CV_median'])}; 100 x 10-fold CV."
                ),
                "source_file": "Fig3_incremental_explanatory_value.tsv",
            }
        )
    for row in table_d:
        keep = (
            row["metric"] == "Primary partial rho"
            and row["analysis_set"] in {
                "Pre-specified high-confidence subset",
                "Maximum-absolute-delta sensitivity",
            }
        )
        if not keep:
            continue
        rows.append(
            {
                "category": "Key sensitivity analysis",
                "contrast": row["contrast"],
                "metric": f"{row['analysis_set']}: primary partial rho",
                "n": row["n_genes"],
                "n_unit": "genes",
                "estimate": row["estimate"],
                "IQR_low": None,
                "IQR_high": None,
                "interval_low": row["CI_low"],
                "interval_high": row["CI_high"],
                "interval_type": "bootstrap 95% CI",
                "P": row["P"],
                "key_note": row["definition"],
                "source_file": (
                    "Fig3_high_confidence_summary.tsv"
                    if row["analysis_set"] == "Pre-specified high-confidence subset"
                    else "Fig3_maxdelta_sensitivity.tsv"
                ),
            }
        )
    return rows


def set_cell_text(cell, text: str, *, bold: bool = False, align: str = "left", size: float = 8.2) -> None:
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.paragraph_format.space_before = Pt(0)
    paragraph.paragraph_format.space_after = Pt(0)
    paragraph.paragraph_format.line_spacing = 1.05
    paragraph.alignment = {
        "left": WD_ALIGN_PARAGRAPH.LEFT,
        "center": WD_ALIGN_PARAGRAPH.CENTER,
        "right": WD_ALIGN_PARAGRAPH.RIGHT,
    }[align]
    run = paragraph.add_run(text)
    run.bold = bold
    run.font.name = "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    run.font.size = Pt(size)
    run.font.color.rgb = RGBColor(0, 0, 0)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def set_cell_margins(cell, top: int = 55, start: int = 75, bottom: int = 55, end: int = 75) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for tag, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{tag}"))
        if node is None:
            node = OxmlElement(f"w:{tag}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_cell_border(cell, *, top: str | None = None, bottom: str | None = None) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = tc_pr.first_child_found_in("w:tcBorders")
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tc_pr.append(borders)
    for edge, style in (("top", top), ("bottom", bottom)):
        if style is None:
            continue
        node = borders.find(qn(f"w:{edge}"))
        if node is None:
            node = OxmlElement(f"w:{edge}")
            borders.append(node)
        node.set(qn("w:val"), style)
        node.set(qn("w:sz"), "10")
        node.set(qn("w:space"), "0")
        node.set(qn("w:color"), "000000")


def set_three_line_table(table) -> None:
    table_pr = table._tbl.tblPr
    borders = table_pr.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        table_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        node = borders.find(qn(f"w:{edge}"))
        if node is None:
            node = OxmlElement(f"w:{edge}")
            borders.append(node)
        node.set(qn("w:val"), "nil")
    for cell in table.rows[0].cells:
        set_cell_border(cell, top="single", bottom="single")
    for cell in table.rows[-1].cells:
        set_cell_border(cell, bottom="single")


def set_repeat_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    tr_pr.append(header)


def set_table_geometry(table, widths: list[int]) -> None:
    total = sum(widths)
    table.autofit = False
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    tbl_pr = table._tbl.tblPr
    tbl_w = tbl_pr.first_child_found_in("w:tblW")
    if tbl_w is None:
        tbl_w = OxmlElement("w:tblW")
        tbl_pr.append(tbl_w)
    tbl_w.set(qn("w:w"), str(total))
    tbl_w.set(qn("w:type"), "dxa")
    tbl_layout = tbl_pr.first_child_found_in("w:tblLayout")
    if tbl_layout is None:
        tbl_layout = OxmlElement("w:tblLayout")
        tbl_pr.append(tbl_layout)
    tbl_layout.set(qn("w:type"), "fixed")
    tbl_ind = tbl_pr.first_child_found_in("w:tblInd")
    if tbl_ind is None:
        tbl_ind = OxmlElement("w:tblInd")
        tbl_pr.append(tbl_ind)
    tbl_ind.set(qn("w:w"), "0")
    tbl_ind.set(qn("w:type"), "dxa")
    grid = table._tbl.tblGrid
    for child in list(grid):
        grid.remove(child)
    for width in widths:
        col = OxmlElement("w:gridCol")
        col.set(qn("w:w"), str(width))
        grid.append(col)
    for row in table.rows:
        for idx, cell in enumerate(row.cells):
            tc_pr = cell._tc.get_or_add_tcPr()
            tc_w = tc_pr.first_child_found_in("w:tcW")
            if tc_w is None:
                tc_w = OxmlElement("w:tcW")
                tc_pr.append(tc_w)
            tc_w.set(qn("w:w"), str(widths[idx]))
            tc_w.set(qn("w:type"), "dxa")
            set_cell_margins(cell)


def add_table_title(doc: Document, label: str, title: str) -> None:
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(7)
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.keep_with_next = True
    run = p.add_run(f"{label}. {title}")
    run.bold = True
    run.font.name = "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    run.font.size = Pt(10.2)


def add_note(doc: Document, text: str) -> None:
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(3)
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.line_spacing = 1.0
    run = p.add_run(text)
    run.font.name = "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    run.font.size = Pt(7.7)


def add_table(doc: Document, headers: list[str], body: list[list[str]], widths: list[int], numeric_cols: Iterable[int]) -> None:
    table = doc.add_table(rows=1, cols=len(headers))
    for idx, header in enumerate(headers):
        set_cell_text(table.rows[0].cells[idx], header, bold=True, align="center", size=8.0)
    numeric_set = set(numeric_cols)
    for row_values in body:
        cells = table.add_row().cells
        for idx, value in enumerate(row_values):
            align = "right" if idx in numeric_set else ("center" if idx in {0, 2, 3} else "left")
            set_cell_text(cells[idx], value, align=align, size=7.9)
    set_repeat_header(table.rows[0])
    set_table_geometry(table, widths)
    set_three_line_table(table)


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.extend([begin, instr, end])


def build_docx(a: list[dict[str, object]], b: list[dict[str, object]], c: list[dict[str, object]], d: list[dict[str, object]]) -> None:
    doc = Document()
    section = doc.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width = Inches(11)
    section.page_height = Inches(8.5)
    section.top_margin = Inches(0.55)
    section.bottom_margin = Inches(0.55)
    section.left_margin = Inches(0.55)
    section.right_margin = Inches(0.55)
    section.header_distance = Inches(0.28)
    section.footer_distance = Inches(0.28)

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    normal.font.size = Pt(9.5)
    normal.paragraph_format.space_after = Pt(4)
    normal.paragraph_format.line_spacing = 1.08

    header = section.header.paragraphs[0]
    header.text = "Supplementary Table Sx | Figure 3 measurement robustness"
    header.alignment = WD_ALIGN_PARAGRAPH.LEFT
    for run in header.runs:
        run.font.name = "Times New Roman"
        run.font.size = Pt(8)
        run.font.color.rgb = RGBColor(90, 90, 90)
    add_page_number(section.footer.paragraphs[0])

    title = doc.add_paragraph()
    title.paragraph_format.space_after = Pt(2)
    run = title.add_run("Supplementary Table Sx | Measurement robustness of matched short-read gene expression and long-read isoform-composition remodeling")
    run.bold = True
    run.font.name = "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    run.font.size = Pt(13.5)

    subtitle = doc.add_paragraph()
    subtitle.paragraph_format.space_after = Pt(7)
    run = subtitle.add_run("Paired primary72 cohort; frozen common testable universe of 1,503 genes")
    run.italic = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor(70, 70, 70)

    add_table_title(doc, "A", "Reliability and primary cross-layer coupling estimates")
    a_body: list[list[str]] = []
    for row in a:
        estimate = fmt(row["estimate"])
        if row["IQR_low"] is not None:
            estimate += f" (IQR {fmt(row['IQR_low'])}-{fmt(row['IQR_high'])})"
        a_body.append(
            [
                str(row["contrast"]).replace("_", " "),
                str(row["measure"]),
                str(row["n"]),
                str(row["n_unit"]),
                estimate,
                fmt_ci(row["interval_low"], row["interval_high"]),
                fmt_p(row["P"]),
                str(row["definition"]),
            ]
        )
    add_table(
        doc,
        ["Contrast", "Measure", "n", "Unit", "Estimate", "Interval", "P", "Definition / interpretation boundary"],
        a_body,
        [900, 2700, 650, 1000, 1550, 1500, 850, 5010],
        numeric_cols={2, 4, 6},
    )
    add_note(
        doc,
        "TVD, total variation distance of within-gene isoform usage. Reliability rows report medians across count-resampling pseudo-replicates or pseudo-libraries; these are not technical replicates. Correlation intervals are bootstrap 95% confidence intervals. The reliability-adjusted estimate is an approximate sensitivity estimate and must not be interpreted as a true deattenuated correlation.",
    )

    doc.add_page_break()
    add_table_title(doc, "B", "Pre-specified abundance and isoform-multiplicity stratification")
    b_body = [
        [
            str(row["stratification"]),
            str(row["contrast"]).replace("_", " "),
            str(row["stratum"]),
            str(row["n_genes"]),
            fmt(row["rho"]),
            fmt_ci(row["CI_low"], row["CI_high"]),
            fmt_p(row["P"]),
        ]
        for row in b
    ]
    add_table(
        doc,
        ["Stratification", "Contrast", "Stratum", "n genes", "Spearman rho", "Bootstrap 95% CI", "P"],
        b_body,
        [2450, 1050, 1100, 900, 1200, 1750, 950],
        numeric_cols={3, 4, 6},
    )
    add_note(
        doc,
        "Long-read abundance strata were defined by mean gene-level long-read count quartiles; short-read strata were defined by DESeq2 baseMean quartiles. Isoform multiplicity denotes the number of eligible long-read isoforms per gene. Strata were pre-specified and do not alter the frozen 1,503-gene universe.",
    )

    doc.add_page_break()
    add_table_title(doc, "C", "Incremental explanatory value of matched short-read gene expression")
    c_body: list[list[str]] = []
    for row in c:
        apparent = f"{fmt(row['delta_apparent'])} (P={fmt_p(row['added_abs_log2FC_P'])})"
        baseline = f"{fmt(row['baseline_CV_median'])} ({fmt(row['baseline_CV_CI_low'])}-{fmt(row['baseline_CV_CI_high'])})"
        extended = f"{fmt(row['extended_CV_median'])} ({fmt(row['extended_CV_CI_low'])}-{fmt(row['extended_CV_CI_high'])})"
        delta = f"{fmt(row['delta_CV_median'])} ({fmt(row['delta_CV_CI_low'])}-{fmt(row['delta_CV_CI_high'])})"
        c_body.append(
            [
                str(row["contrast"]).replace("_", " "),
                str(row["model"]).replace("_", " "),
                str(row["n_genes"]),
                str(row["apparent_metric"]),
                apparent,
                baseline,
                extended,
                delta,
                f"{row['CV_repeats']} x {row['CV_folds']}-fold",
            ]
        )
    add_table(
        doc,
        ["Contrast", "Model", "n genes", "Metric", "Apparent gain", "CV baseline", "CV extended", "CV gain", "Resampling"],
        c_body,
        [850, 1350, 700, 1500, 1450, 1550, 1550, 1550, 1200],
        numeric_cols={2},
    )
    add_note(
        doc,
        "The baseline model included long-read abundance, short-read abundance, isoform multiplicity, and zero fraction; the extended model additionally included matched short-read |log2FC|. Rank-linear explained variance was the primary model and quasibinomial explained deviance was a model-form sensitivity analysis. CV summaries are medians and 2.5th-97.5th percentiles across 100 repetitions of 10-fold cross-validation.",
    )

    add_table_title(doc, "D", "High-confidence and outcome-definition sensitivity analyses")
    d_body = [
        [
            str(row["analysis_set"]),
            str(row["contrast"]).replace("_", " "),
            str(row["n_genes"]),
            str(row["metric"]),
            fmt(row["estimate"]),
            fmt_ci(row["CI_low"], row["CI_high"]),
            fmt_p(row["P"]),
            str(row["definition"]),
        ]
        for row in d
    ]
    add_table(
        doc,
        ["Analysis set", "Contrast", "n genes", "Metric", "Estimate", "95% interval", "P", "Definition"],
        d_body,
        [2300, 950, 700, 2200, 950, 1450, 900, 3810],
        numeric_cols={2, 4, 6},
    )
    add_note(
        doc,
        "The high-confidence subset required long-read mean gene count at or above the common-universe median, eligible isoform multiplicity <=5, and zero fraction at or below the common-universe median. The maximum-absolute-delta analysis replaced TVD with the largest absolute isoform-usage change per gene. These analyses test robustness and do not redefine the primary universe or thresholds.",
    )

    add_note(
        doc,
        "Overall interpretation: matched short-read gene-abundance change and long-read isoform-composition remodeling showed limited but non-zero coupling after reliability, abundance, isoform multiplicity, and zero-count structure were considered. The matched short-read data were generated from aliquots of the same RNA extracts and therefore provide complementary measurement context, not independent external validation. Formal Figure 3 and its frozen source data were not modified.",
    )

    doc.core_properties.title = "Supplementary Table Sx - Figure 3 measurement robustness"
    doc.core_properties.subject = "CRC paired N/A/T long-read and matched short-read robustness analysis"
    doc.core_properties.author = "CRC N/A/T isoform remodeling project"
    doc.core_properties.keywords = "CRC; isoform usage; long-read; matched short-read; robustness"
    DOCX_PATH.parent.mkdir(parents=True, exist_ok=True)
    doc.save(DOCX_PATH)


def build_concise_docx(rows: list[dict[str, object]]) -> None:
    doc = Document()
    section = doc.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width = Inches(11)
    section.page_height = Inches(8.5)
    section.top_margin = Inches(0.55)
    section.bottom_margin = Inches(0.55)
    section.left_margin = Inches(0.55)
    section.right_margin = Inches(0.55)
    section.header_distance = Inches(0.28)
    section.footer_distance = Inches(0.28)

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    normal.font.size = Pt(9.2)
    normal.paragraph_format.space_after = Pt(4)
    normal.paragraph_format.line_spacing = 1.05

    header = section.header.paragraphs[0]
    header.text = "Supplementary Table Sx | Figure 3 key robustness results"
    for run in header.runs:
        run.font.name = "Times New Roman"
        run.font.size = Pt(8)
        run.font.color.rgb = RGBColor(90, 90, 90)
    add_page_number(section.footer.paragraphs[0])

    title = doc.add_paragraph()
    title.paragraph_format.space_after = Pt(2)
    run = title.add_run(
        "Supplementary Table Sx | Key measurement-robustness results for matched short-read gene expression and long-read isoform remodeling"
    )
    run.bold = True
    run.font.name = "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    run.font.size = Pt(13.2)

    subtitle = doc.add_paragraph()
    subtitle.paragraph_format.space_after = Pt(6)
    run = subtitle.add_run("Paired primary72 cohort; frozen common testable universe of 1,503 genes")
    run.italic = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor(70, 70, 70)

    body: list[list[str]] = []
    for row in rows:
        estimate = fmt(row["estimate"])
        if row["IQR_low"] is not None:
            estimate += f" (IQR {fmt(row['IQR_low'])}-{fmt(row['IQR_high'])})"
        body.append(
            [
                str(row["category"]),
                str(row["contrast"]).replace("_", " "),
                str(row["metric"]),
                f"{row['n']} {row['n_unit']}",
                estimate,
                fmt_ci(row["interval_low"], row["interval_high"]),
                fmt_p(row["P"]),
                str(row["key_note"]),
            ]
        )
    add_table(
        doc,
        ["Category", "Contrast", "Metric", "Analysis units", "Estimate", "Interval", "P", "Key definition / boundary"],
        body,
        [1700, 850, 3000, 1150, 1350, 1450, 850, 3810],
        numeric_cols={4, 6},
    )
    add_note(
        doc,
        "TVD, total variation distance of within-gene isoform usage; CV, cross-validation. Reliability values are medians across count-resampling pseudo-replicates or pseudo-libraries and are not technical-replicate estimates. Reliability intervals are 2.5th-97.5th percentiles; correlation intervals are bootstrap 95% confidence intervals. The reliability-adjusted rho is an approximate sensitivity estimate, not a true deattenuated correlation.",
    )
    add_note(
        doc,
        "The primary partial correlation adjusted normal-score ranks for long-read count, short-read baseMean, isoform multiplicity, and zero fraction. Incremental explained variance used 100 repetitions of 10-fold cross-validation. The high-confidence subset and maximum-absolute-delta analysis were pre-specified sensitivity checks and did not redefine the frozen primary universe or thresholds.",
    )
    add_note(
        doc,
        "Interpretation boundary: the results support limited but non-zero cross-layer coupling. Matched short-read data were obtained from aliquots of the same RNA extracts and provide complementary measurement context, not independent external validation. Formal Figure 3 and its source data were not modified.",
    )

    doc.core_properties.title = "Supplementary Table Sx - Figure 3 key measurement robustness results"
    doc.core_properties.subject = "CRC paired N/A/T long-read and matched short-read robustness analysis"
    doc.core_properties.author = "CRC N/A/T isoform remodeling project"
    doc.core_properties.keywords = "CRC; isoform usage; long-read; matched short-read; robustness"
    DOCX_PATH.parent.mkdir(parents=True, exist_ok=True)
    doc.save(DOCX_PATH)


def build_report(source_files: list[str], table_paths: list[Path]) -> None:
    lines = [
        "# Figure 3 robustness 补充三线表导出报告",
        "",
        "## 定位",
        "",
        "- 本轮未重跑任何统计模型，仅将已冻结的 Figure 3 measurement robustness 结果重排为发表用补充表。",
        "- 正式 Figure 3 及其 source data 未修改。",
        "- 原 robustness candidate/supplementary figures 仅保留为内部 QA 和分析追溯材料，**不纳入 Figure 3 补充图，不用于投稿**。",
        "- 根据最终审阅要求，发表版已压缩为 **1 张核心三线表**；不呈现 24 行分层明细、500/200 次逐次模拟结果或内部 LEVEL 判定。",
        "",
        "## 发表表格结构",
        "",
        "- 主表仅保留：LR/SR 内部可靠性、observed/reliability-adjusted/primary partial coupling、rank-linear CV 增量解释度、high-confidence partial rho 和 max-delta partial rho。",
        "- 丰度四分位、isoform multiplicity 分层和 quasibinomial 完整明细仍保留在冻结 results 中，不进入发表表。",
        "",
        "## 格式",
        "",
        "- DOCX 采用 Times New Roman、横向 Letter 页面和无竖线的真三线表（表顶线、表头底线、表底线）。",
        "- XLSX 保留未舍入数值、来源字段和结论边界，用于编辑、审稿和最终编号。",
        "- 500/200 次逐次 pseudo-replicate/pseudo-library 明细未重复放入发表表中；原始明细仍保留在冻结 results 目录。",
        "",
        "## 科学边界",
        "",
        "- pseudo-replicates/pseudo-libraries 不是 technical replicates。",
        "- reliability-adjusted rho 是近似敏感性估计，不是真实或去噪相关。",
        "- 低相关只支持 limited/incomplete coupling，不支持两个生物学层面完全独立。",
        "- matched short-read 来自同一批 RNA 提取物分装，是 complementary measurement context，不是独立外部验证。",
        "",
        "## 输出",
        "",
    ]
    for path in table_paths:
        lines.append(f"- `{path}`")
    lines.extend(["", "## 冻结来源", ""])
    for name in source_files:
        lines.append(f"- `{RESULTS / name}`")
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    source_files = [
        "Fig3_LR_TVD_reliability_summary.tsv",
        "Fig3_SR_log2FC_reliability_summary.tsv",
        "Fig3_reliability_adjusted_coupling.tsv",
        "Fig3_partial_correlation.tsv",
        "Fig3_incremental_explanatory_value.tsv",
        "Fig3_high_confidence_summary.tsv",
        "Fig3_maxdelta_sensitivity.tsv",
    ]
    for name in source_files:
        read_tsv(name)

    table_a = build_table_a()
    table_c = build_table_c()
    table_d = build_table_d()
    concise_rows = build_concise_rows(table_a, table_c, table_d)

    outputs = {
        "Supplementary_Table_Sx_key_measurement_robustness.tsv": (
            concise_rows,
            [
                "category", "contrast", "metric", "n", "n_unit", "estimate", "IQR_low", "IQR_high",
                "interval_low", "interval_high", "interval_type", "P", "key_note", "source_file",
            ],
        ),
    }
    table_paths: list[Path] = []
    for name, (rows, columns) in outputs.items():
        path = SOURCE / name
        write_tsv(path, rows, columns)
        table_paths.append(path)

    build_concise_docx(concise_rows)
    table_paths.insert(0, DOCX_PATH)
    build_report(source_files, table_paths)
    table_paths.append(REPORT_PATH)

    manifest_rows = []
    for path in [*(RESULTS / name for name in source_files), *table_paths]:
        manifest_rows.append(
            {
                "file_role": "frozen_input" if path.parent == RESULTS else "publication_table_output",
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    write_tsv(MANIFEST_PATH, manifest_rows, ["file_role", "path", "size_bytes", "sha256"])
    print(f"DOCX={DOCX_PATH}")
    print(f"REPORT={REPORT_PATH}")
    print(f"SOURCE_TABLES={len(outputs)}")


if __name__ == "__main__":
    main()
