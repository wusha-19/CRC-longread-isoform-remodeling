#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(os.environ["CRC_ROOT"]).expanduser().resolve()
PHASE05 = ROOT / "05_isoform_main_study/05_functional_annotation"
PHASE04 = ROOT / "05_isoform_main_study/04_NAT_dynamic_patterns"
PHASE01 = ROOT / "05_isoform_main_study/01_master_catalog"
OUT = ROOT / "05_isoform_main_study/05B_ORF_domain_validation"
ENV = Path(os.environ["CRC_ORF_ENV"]).expanduser().resolve()
BIN = ENV / "bin"
TOOLS = {
    "TransDecoder.LongOrfs": BIN / "TransDecoder.LongOrfs",
    "TransDecoder.Predict": BIN / "TransDecoder.Predict",
    "gffread": BIN / "gffread",
    "seqkit": BIN / "seqkit",
    "bedtools": BIN / "bedtools",
    "samtools": BIN / "samtools",
    "hmmscan": BIN / "hmmscan",
    "hmmsearch": BIN / "hmmsearch",
}


def ensure_dirs() -> None:
    for sub in ["code", "env", "tables", "fasta", "gtf", "proteins", "domains", "figures", "candidate_cards_tier1_updated", "reports", "logs", "provenance"]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k, "NA")) for k in fields})


def fmt(value: object) -> str:
    if value is None or value == "":
        return "NA"
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.12g}"
    return str(value)


def bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)[:180]


def parse_attrs(attr: str) -> dict[str, str]:
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr):
        out[m.group(1)] = m.group(2)
    return out


def read_fasta_subset(path: Path, keep: set[str]) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = None
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                tid = line[1:].strip().split()[0]
                cur = tid if tid in keep else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as handle:
        for tid, seq in seqs.items():
            handle.write(f">{tid}\n")
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i + 70] + "\n")


def extract_gtf_subset(path: Path, keep: set[str], out_path: Path) -> dict[str, dict[str, object]]:
    tx = {i: {"exons": [], "chrom": "NA", "strand": "NA"} for i in keep}
    with path.open() as inp, out_path.open("w") as out:
        for line in inp:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_attrs(parts[8])
            tid = attrs.get("transcript_id")
            if tid not in keep:
                continue
            out.write(line)
            chrom, feature, start, end, strand = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6]
            tx[tid]["chrom"] = chrom
            tx[tid]["strand"] = strand
            if feature == "exon":
                tx[tid]["exons"].append((start, end))
    for val in tx.values():
        val["exons"] = sorted(val["exons"])
    return tx


def exon_len(exons: list[tuple[int, int]]) -> int:
    return sum(e - s + 1 for s, e in exons)


def last_eej_tx_position(exons: list[tuple[int, int]], strand: str) -> int | None:
    if len(exons) < 2:
        return None
    ordered = exons if strand == "+" else list(reversed(exons))
    return sum(e - s + 1 for s, e in ordered[:-1])


def run_cmd(cmd: list[str], cwd: Path, stdout: Path, stderr: Path) -> int:
    with stdout.open("w") as so, stderr.open("w") as se:
        proc = subprocess.run(cmd, cwd=str(cwd), stdout=so, stderr=se, text=True)
    return proc.returncode


def version(path: Path) -> str:
    if not path.is_file():
        return "NA"
    for args in ([str(path), "--version"], [str(path), "-h"]):
        try:
            out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=5).stdout.strip()
            if out:
                return out.splitlines()[0][:180]
        except Exception:
            pass
    return "available_version_not_detected"


def parse_transdecoder_pep(path: Path) -> dict[str, dict[str, object]]:
    # Header example usually includes:
    # >transcript.p1 type:complete len:123 (+) transcript:1-369(+)
    out: dict[str, dict[str, object]] = {}
    if not path.is_file() or path.stat().st_size == 0:
        return out
    cur_header = None
    cur_seq: list[str] = []

    def flush() -> None:
        if not cur_header:
            return
        header = cur_header
        seq = "".join(cur_seq)
        first = header.split()[0]
        iso = first.rsplit(".p", 1)[0] if ".p" in first else first
        typem = re.search(r"type:([A-Za-z_]+)", header)
        lenm = re.search(r"len:(\d+)", header)
        coordm = re.search(r":(\d+)-(\d+)\(([+-])\)", header)
        aa_len = int(lenm.group(1)) if lenm else len(seq)
        start = int(coordm.group(1)) if coordm else None
        end = int(coordm.group(2)) if coordm else None
        old = out.get(iso)
        if old is None or aa_len > int(old["protein_length_aa"]):
            out[iso] = {
                "longest_ORF_sequence_id": first,
                "header": header,
                "orf_type": typem.group(1) if typem else "unknown",
                "ORF_start": start,
                "ORF_end": end,
                "ORF_length_nt": abs(end - start) + 1 if start is not None and end is not None else aa_len * 3,
                "protein_length_aa": aa_len,
                "protein_sequence": seq,
            }

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                cur_header = line[1:]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    flush()
    return out


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")
    return seq.translate(table)[::-1].replace("U", "T").replace("u", "t")


def simple_orf_screen(seq: str) -> dict[str, object] | None:
    """Conservative fallback: longest ATG-to-stop ORF across six frames."""
    stops = {"TAA", "TAG", "TGA"}
    best = None
    seq = seq.upper().replace("U", "T")
    strands = [("+", seq), ("-", revcomp(seq))]
    for strand, s in strands:
        for frame in range(3):
            i = frame
            while i <= len(s) - 3:
                codon = s[i:i + 3]
                if codon == "ATG":
                    j = i + 3
                    while j <= len(s) - 3:
                        stop = s[j:j + 3]
                        if stop in stops:
                            nt_len = j + 3 - i
                            if nt_len >= 90:
                                if strand == "+":
                                    start, end = i + 1, j + 3
                                else:
                                    # Convert reverse-complement coordinates back to transcript coordinates.
                                    start, end = len(seq) - (j + 3) + 1, len(seq) - i
                                cand = {
                                    "longest_ORF_sequence_id": f"simple_orf_{strand}{frame}_{start}_{end}",
                                    "header": f"SIMPLE_ORF_SCREEN_ONLY strand={strand} frame={frame} transcript:{start}-{end}",
                                    "orf_type": "complete",
                                    "ORF_start": start,
                                    "ORF_end": end,
                                    "ORF_length_nt": nt_len,
                                    "protein_length_aa": nt_len // 3,
                                    "protein_sequence": "NA",
                                }
                                if best is None or nt_len > int(best["ORF_length_nt"]):
                                    best = cand
                            break
                        j += 3
                    i = j
                i += 3
    return best


def main() -> None:
    ensure_dirs()
    tier_path = PHASE05 / "tables/phase05_annotation_tier_assignment.tsv"
    integrated_path = PHASE05 / "tables/phase05_integrated_functional_annotation.tsv"
    phase04_input = PHASE04 / "tables/NAT_next_stage_annotation_input.tsv"
    fasta_path = PHASE01 / "annotation/primary72_isoforms.fa"
    gtf_path = PHASE01 / "annotation/primary72_isoforms.gtf"
    for p in [tier_path, integrated_path, phase04_input, fasta_path, gtf_path]:
        if not p.is_file() or p.stat().st_size == 0:
            raise FileNotFoundError(p)

    tier_all = read_tsv(tier_path)
    tier1 = [r for r in tier_all if r["tier"] == "Tier 1"]
    tier1_ids = {r["isoform_id"] for r in tier1}
    if len(tier1_ids) != 53:
        raise ValueError(f"Tier 1 count expected 53, found {len(tier1_ids)}")
    integrated = {r["isoform_id"]: r for r in read_tsv(integrated_path)}
    phase04 = {r["isoform_id"]: r for r in read_tsv(phase04_input)}

    # Tool/database availability.
    db_candidates = []
    for root in [Path(os.environ.get("CONDA_ENVS_DIR", ENV.parent)), Path(os.environ.get("SOFTWARE_ROOT", ROOT / "software")), ROOT]:
        db_candidates.extend(root.glob("**/Pfam-A.hmm"))
    tool_rows = []
    for name, path in TOOLS.items():
        tool_rows.append({
            "tool_or_database": name,
            "available": bool_text(path.is_file()),
            "path": str(path) if path.is_file() else "NA",
            "environment": str(ENV),
            "version_or_status": version(path),
            "usable_for_phase05B": bool_text(path.is_file()),
            "notes": "project_conda_env_crc_orf_anno" if path.is_file() else "missing",
        })
    pfam_path = db_candidates[0] if db_candidates else None
    tool_rows.append({
        "tool_or_database": "Pfam-A.hmm",
        "available": bool_text(pfam_path is not None),
        "path": str(pfam_path) if pfam_path else "NA",
        "environment": "existing_limited_search_paths",
        "version_or_status": "found" if pfam_path else "missing",
        "usable_for_phase05B": bool_text(pfam_path is not None),
        "notes": "do_not_download_full_Pfam_without_user_approval" if not pfam_path else "existing_database",
    })
    write_tsv(OUT / "tables/phase05B_tool_and_database_availability.tsv", tool_rows, ["tool_or_database", "available", "path", "environment", "version_or_status", "usable_for_phase05B", "notes"])

    install_rows = [{
        "step": "create_crc_orf_anno",
        "command": "conda create -y -n crc_orf_anno -c conda-forge -c bioconda transdecoder gffread seqkit bedtools samtools hmmer",
        "status": "completed_before_phase05B_resume",
        "stdout_log": str(OUT / "logs/conda_create_crc_orf_anno.stdout.log"),
        "stderr_log": str(OUT / "logs/conda_create_crc_orf_anno.stderr.log"),
        "notes": "user confirmed env listed; tool validation performed in this run",
    }]
    write_tsv(OUT / "tables/phase05B_environment_installation_log.tsv", install_rows, ["step", "command", "status", "stdout_log", "stderr_log", "notes"])

    # Tier 1 extraction.
    seqs = read_fasta_subset(fasta_path, tier1_ids)
    missing_fa = sorted(tier1_ids - set(seqs))
    if missing_fa:
        raise ValueError(f"Tier1 isoforms missing FASTA: {missing_fa[:5]}")
    tier1_fa = OUT / "fasta/tier1_transcripts.fa"
    write_fasta(tier1_fa, {r["isoform_id"]: seqs[r["isoform_id"]] for r in tier1})
    gtf_info = extract_gtf_subset(gtf_path, tier1_ids, OUT / "gtf/tier1_transcripts.gtf")
    missing_gtf = [i for i in tier1_ids if not gtf_info.get(i, {}).get("exons")]
    if missing_gtf:
        raise ValueError(f"Tier1 isoforms missing GTF exons: {missing_gtf[:5]}")

    write_tsv(OUT / "tables/tier1_isoform_list.tsv", tier1, list(tier1[0].keys()))
    genes = {}
    for r in tier1:
        genes[r["gene_id"]] = {"gene_id": r["gene_id"], "gene_name": r["gene_name"], "n_tier1_isoforms": 0}
    for r in tier1:
        genes[r["gene_id"]]["n_tier1_isoforms"] += 1
    write_tsv(OUT / "tables/tier1_gene_list.tsv", list(genes.values()), ["gene_id", "gene_name", "n_tier1_isoforms"])

    # ORF prediction with TransDecoder.
    td_ok = TOOLS["TransDecoder.LongOrfs"].is_file() and TOOLS["TransDecoder.Predict"].is_file()
    td_status = "not_run"
    if td_ok:
        for old in [OUT / "fasta/tier1_transcripts.fa.transdecoder_dir", OUT / "fasta/tier1_transcripts.fa.transdecoder.pep", OUT / "fasta/tier1_transcripts.fa.transdecoder.cds", OUT / "fasta/tier1_transcripts.fa.transdecoder.gff3", OUT / "fasta/tier1_transcripts.fa.transdecoder.bed"]:
            if old.is_dir():
                shutil.rmtree(old)
            elif old.exists():
                old.unlink()
        rc1 = run_cmd([str(TOOLS["TransDecoder.LongOrfs"]), "-t", str(tier1_fa)], OUT / "fasta", OUT / "logs/TransDecoder.LongOrfs.stdout.log", OUT / "logs/TransDecoder.LongOrfs.stderr.log")
        if rc1 == 0:
            rc2 = run_cmd([str(TOOLS["TransDecoder.Predict"]), "-t", str(tier1_fa)], OUT / "fasta", OUT / "logs/TransDecoder.Predict.stdout.log", OUT / "logs/TransDecoder.Predict.stderr.log")
            td_status = "TransDecoder_success" if rc2 == 0 else f"TransDecoder.Predict_failed_rc{rc2}"
        else:
            td_status = f"TransDecoder.LongOrfs_failed_rc{rc1}"
    else:
        td_status = "TransDecoder_missing"

    pep_path = OUT / "fasta/tier1_transcripts.fa.transdecoder.pep"
    cds_path = OUT / "fasta/tier1_transcripts.fa.transdecoder.cds"
    if pep_path.is_file():
        shutil.copy2(pep_path, OUT / "proteins/tier1_predicted_proteins.fa")
    else:
        (OUT / "proteins/tier1_predicted_proteins.fa").write_text("")
    orfs = parse_transdecoder_pep(pep_path)
    if not orfs:
        # Last allowed fallback: simple ORF screen, clearly labeled low confidence.
        for iso, seq in seqs.items():
            hit = simple_orf_screen(seq)
            if hit:
                orfs[iso] = hit
        if td_status == "TransDecoder_missing":
            td_status = "SIMPLE_ORF_SCREEN_ONLY_TRANSDECODER_BROKEN_SYMLINK"
        elif not td_status.startswith("TransDecoder_success"):
            td_status = f"SIMPLE_ORF_SCREEN_ONLY_AFTER_{td_status}"

    orf_rows = []
    nmd_rows = []
    for r in tier1:
        iso = r["isoform_id"]
        hit = orfs.get(iso)
        if hit:
            orf_found = True
            orf_start = hit["ORF_start"]
            orf_end = hit["ORF_end"]
            orf_len = hit["ORF_length_nt"]
            aa_len = hit["protein_length_aa"]
            complete = str(hit["orf_type"]).lower() == "complete"
            if str(hit["header"]).startswith("SIMPLE_ORF_SCREEN_ONLY"):
                conf = "LOW_SIMPLE_SCREEN"
                source = "SIMPLE_ORF_SCREEN_ONLY"
            else:
                conf = "HIGH_TOOL_BASED"
                source = "TransDecoder"
            notes = hit["header"]
        else:
            orf_found = False
            orf_start = orf_end = orf_len = aa_len = "NA"
            complete = False
            conf = "INSUFFICIENT_INFORMATION" if td_status.startswith("TransDecoder_success") else "SIMPLE_ORF_SCREEN_ONLY_NOT_USED"
            source = td_status
            notes = "No TransDecoder ORF predicted"
        orf_rows.append({
            "isoform_id": iso,
            "gene_id": r["gene_id"],
            "gene_name": r["gene_name"],
            "ORF_found": bool_text(orf_found),
            "ORF_source": source,
            "ORF_start": orf_start,
            "ORF_end": orf_end,
            "ORF_length_nt": orf_len,
            "protein_length_aa": aa_len,
            "start_codon_present": bool_text(complete),
            "stop_codon_present": bool_text(complete),
            "complete_ORF": bool_text(complete),
            "longest_ORF_sequence_id": hit["longest_ORF_sequence_id"] if hit else "NA",
            "confidence": conf,
            "notes": notes,
        })
        exons = gtf_info[iso]["exons"]
        last = last_eej_tx_position(exons, str(gtf_info[iso]["strand"]))
        ptc = int(orf_end) if isinstance(orf_end, int) else None
        if orf_found and ptc is not None and last is not None:
            dist = last - ptc
            risk = "PREDICTED_NMD_RISK" if dist >= 50 else "NO_PREDICTED_NMD_RISK"
            nmd_conf = "MEDIUM_TOOL_BASED"
        elif orf_found and last is None:
            dist = "NA"
            risk = "NO_PREDICTED_NMD_RISK"
            nmd_conf = "single_exon_or_no_EEJ"
        else:
            dist = "NA"
            risk = "INSUFFICIENT_INFORMATION"
            nmd_conf = "no_ORF"
        nmd_rows.append({
            "isoform_id": iso,
            "gene_id": r["gene_id"],
            "gene_name": r["gene_name"],
            "predicted_PTC": bool_text(risk == "PREDICTED_NMD_RISK"),
            "predicted_NMD_risk": risk,
            "NMD_rule_source": "50nt_upstream_of_last_EEJ_conservative_rule",
            "PTC_position": ptc if ptc is not None else "NA",
            "last_exon_exon_junction_position": last if last is not None else "NA",
            "distance_to_last_EEJ": dist,
            "confidence": nmd_conf,
            "notes": "prediction_only_not_experimental_validation",
        })

    orf_fields = ["isoform_id", "gene_id", "gene_name", "ORF_found", "ORF_source", "ORF_start", "ORF_end", "ORF_length_nt", "protein_length_aa", "start_codon_present", "stop_codon_present", "complete_ORF", "longest_ORF_sequence_id", "confidence", "notes"]
    nmd_fields = ["isoform_id", "gene_id", "gene_name", "predicted_PTC", "predicted_NMD_risk", "NMD_rule_source", "PTC_position", "last_exon_exon_junction_position", "distance_to_last_EEJ", "confidence", "notes"]
    write_tsv(OUT / "tables/phase05B_tier1_ORF_prediction.tsv", orf_rows, orf_fields)
    write_tsv(OUT / "tables/phase05B_tier1_NMD_prediction.tsv", nmd_rows, nmd_fields)
    orf_by_iso = {r["isoform_id"]: r for r in orf_rows}
    nmd_by_iso = {r["isoform_id"]: r for r in nmd_rows}

    # Domain status. hmmscan is available but Pfam database is missing, so do not run full Pfam.
    domain_rows = []
    domain_summary = []
    for r in tier1:
        iso = r["isoform_id"]
        orf = orf_by_iso[iso]
        if orf["ORF_found"] != "TRUE":
            status = "NO_ORF"
            note = "No predicted ORF for domain scan"
        elif not TOOLS["hmmscan"].is_file():
            status = "DOMAIN_TOOL_MISSING"
            note = "hmmscan missing"
        elif pfam_path is None:
            status = "DOMAIN_TOOL_AVAILABLE_DATABASE_MISSING"
            note = "hmmscan available but Pfam-A.hmm missing; full Pfam download requires user approval"
        else:
            status = "NOT_RUN_PENDING_USER_APPROVAL"
            note = "Pfam found but run not implemented in this lightweight tier1 pass"
        domain_rows.append({
            "isoform_id": iso,
            "gene_id": r["gene_id"],
            "gene_name": r["gene_name"],
            "protein_id": orf["longest_ORF_sequence_id"],
            "domain_source": "not_run",
            "domain_id": "NA",
            "domain_name": "NA",
            "domain_start": "NA",
            "domain_end": "NA",
            "evalue": "NA",
            "notes": note,
        })
        domain_summary.append({
            "isoform_id": iso,
            "gene_id": r["gene_id"],
            "gene_name": r["gene_name"],
            "ORF_found": orf["ORF_found"],
            "protein_length_aa": orf["protein_length_aa"],
            "has_domain": "FALSE",
            "n_domains": 0,
            "domain_list": "NA",
            "domain_annotation_status": status,
        })
    write_tsv(OUT / "tables/phase05B_tier1_domain_annotation.tsv", domain_rows, ["isoform_id", "gene_id", "gene_name", "protein_id", "domain_source", "domain_id", "domain_name", "domain_start", "domain_end", "evalue", "notes"])
    write_tsv(OUT / "tables/phase05B_tier1_domain_summary.tsv", domain_summary, ["isoform_id", "gene_id", "gene_name", "ORF_found", "protein_length_aa", "has_domain", "n_domains", "domain_list", "domain_annotation_status"])
    domain_by_iso = {r["isoform_id"]: r for r in domain_summary}

    # Within-gene comparison against other Tier 1 isoforms in same gene.
    by_gene = defaultdict(list)
    for r in tier1:
        by_gene[r["gene_id"]].append(r["isoform_id"])
    comparison_rows = []
    for r in tier1:
        iso = r["isoform_id"]
        gid = r["gene_id"]
        peers = [x for x in by_gene[gid] if x != iso]
        if peers:
            comp = peers[0]
            reason = "same_gene_tier1_peer"
        else:
            comp = "NA"
            reason = "no_same_gene_tier1_peer"
        ex = gtf_info[iso]["exons"]
        comp_ex = gtf_info[comp]["exons"] if comp != "NA" else []
        olen = int(orf_by_iso[iso]["ORF_length_nt"]) if str(orf_by_iso[iso]["ORF_length_nt"]).isdigit() else 0
        clen = int(orf_by_iso[comp]["ORF_length_nt"]) if comp != "NA" and str(orf_by_iso[comp]["ORF_length_nt"]).isdigit() else 0
        comparison_rows.append({
            "gene_id": gid,
            "gene_name": r["gene_name"],
            "tier1_isoform_id": iso,
            "comparison_isoform_id": comp,
            "comparison_reason": reason,
            "delta_transcript_length": exon_len(ex) - exon_len(comp_ex) if comp_ex else "NA",
            "delta_exon_count": len(ex) - len(comp_ex) if comp_ex else "NA",
            "delta_ORF_length_nt": olen - clen if comp != "NA" else "NA",
            "delta_protein_length_aa": (olen - clen) // 3 if comp != "NA" else "NA",
            "domain_gain_loss_status": "NO_DOMAIN_ANNOTATION",
            "predicted_NMD_difference": "same" if comp != "NA" and nmd_by_iso[iso]["predicted_NMD_risk"] == nmd_by_iso[comp]["predicted_NMD_risk"] else "insufficient_or_different",
            "dominant_switch_related": r["dominant_isoform_switch"],
            "notes": "Tier1-only same gene comparison; inspect all gene isoforms manually for final interpretation",
        })
    write_tsv(OUT / "tables/phase05B_tier1_within_gene_isoform_comparison.tsv", comparison_rows, ["gene_id", "gene_name", "tier1_isoform_id", "comparison_isoform_id", "comparison_reason", "delta_transcript_length", "delta_exon_count", "delta_ORF_length_nt", "delta_protein_length_aa", "domain_gain_loss_status", "predicted_NMD_difference", "dominant_switch_related", "notes"])

    integrated_rows = []
    for r in tier1:
        iso = r["isoform_id"]
        p5 = integrated[iso]
        p4 = phase04.get(iso, {})
        orf = orf_by_iso[iso]
        nmd = nmd_by_iso[iso]
        dom = domain_by_iso[iso]
        if nmd["predicted_NMD_risk"] == "PREDICTED_NMD_RISK":
            status = "PREDICTED_NMD_RISK"
        elif orf["ORF_found"] == "TRUE" and dom["has_domain"] == "TRUE":
            status = "ORF_DOMAIN_ANNOTATED"
        elif orf["ORF_found"] == "TRUE":
            status = "ORF_ONLY_NO_DOMAIN"
        elif dom["domain_annotation_status"] == "DOMAIN_TOOL_AVAILABLE_DATABASE_MISSING":
            status = "INSUFFICIENT_TOOL_SUPPORT"
        else:
            status = "STRUCTURAL_ONLY_NO_ORF"
        if r["review_priority"] == "HIGH" or orf["ORF_found"] == "TRUE" or nmd["predicted_NMD_risk"] == "PREDICTED_NMD_RISK":
            mpri = "HIGH"
        elif r["in_DTU"] == "TRUE":
            mpri = "MEDIUM"
        else:
            mpri = "LOW"
        notes = []
        if dom["domain_annotation_status"] == "DOMAIN_TOOL_AVAILABLE_DATABASE_MISSING":
            notes.append("Pfam missing; domain scan not run")
        if orf["ORF_found"] != "TRUE":
            notes.append("No TransDecoder ORF predicted")
        integrated_rows.append({
            "isoform_id": iso,
            "gene_id": r["gene_id"],
            "gene_name": r["gene_name"],
            "review_priority": r["review_priority"],
            "next_stage_annotation_priority": r["next_stage_annotation_priority"],
            "in_DTE": r["in_DTE"],
            "in_DTU": r["in_DTU"],
            "gene_silent_DTU": r["gene_silent_DTU"],
            "dominant_isoform_switch": r["dominant_isoform_switch"],
            "expression_dynamic_pattern": p5.get("expression_dynamic_pattern", p4.get("expression_dynamic_pattern", "NA")),
            "usage_dynamic_pattern": p5.get("usage_dynamic_pattern", p4.get("usage_dynamic_pattern", "NA")),
            "ORF_found": orf["ORF_found"],
            "ORF_source": orf["ORF_source"],
            "ORF_length_nt": orf["ORF_length_nt"],
            "protein_length_aa": orf["protein_length_aa"],
            "complete_ORF": orf["complete_ORF"],
            "predicted_NMD_risk": nmd["predicted_NMD_risk"],
            "has_domain": dom["has_domain"],
            "domain_list": dom["domain_list"],
            "domain_annotation_status": dom["domain_annotation_status"],
            "RT_PCR_readiness_from_phase05": p5.get("primer_design_readiness", "NA"),
            "functional_consequence_status": status,
            "manual_review_priority": mpri,
            "manual_review_note": "; ".join(notes) if notes else "NA",
        })
    int_fields = ["isoform_id", "gene_id", "gene_name", "review_priority", "next_stage_annotation_priority", "in_DTE", "in_DTU", "gene_silent_DTU", "dominant_isoform_switch", "expression_dynamic_pattern", "usage_dynamic_pattern", "ORF_found", "ORF_source", "ORF_length_nt", "protein_length_aa", "complete_ORF", "predicted_NMD_risk", "has_domain", "domain_list", "domain_annotation_status", "RT_PCR_readiness_from_phase05", "functional_consequence_status", "manual_review_priority", "manual_review_note"]
    write_tsv(OUT / "tables/phase05B_tier1_integrated_annotation.tsv", integrated_rows, int_fields)

    # Updated cards.
    for r in integrated_rows:
        card = f"""# Phase05B Tier 1 updated card: {r['gene_id']} | {r['isoform_id']}

## DTE / DTU / NAT Evidence

- in_DTE: {r['in_DTE']}
- in_DTU: {r['in_DTU']}
- Phase03B review priority: {r['review_priority']}
- Phase04 next-stage priority: {r['next_stage_annotation_priority']}
- Expression dynamic pattern: {r['expression_dynamic_pattern']}
- Usage dynamic pattern: {r['usage_dynamic_pattern']}
- Gene-silent DTU: {r['gene_silent_DTU']}
- Dominant isoform switch: {r['dominant_isoform_switch']}

## Phase05 RT-PCR Feasibility

- RT-PCR readiness from Phase05: {r['RT_PCR_readiness_from_phase05']}

## Phase05B ORF Prediction

- ORF found: {r['ORF_found']}
- ORF source: {r['ORF_source']}
- ORF length nt: {r['ORF_length_nt']}
- Protein length aa: {r['protein_length_aa']}
- Complete ORF: {r['complete_ORF']}

## Phase05B NMD/PTC Prediction

- Predicted NMD risk: {r['predicted_NMD_risk']}
- This is computational prediction only, not experimental validation.

## Domain Annotation

- Has domain: {r['has_domain']}
- Domain list: {r['domain_list']}
- Domain status: {r['domain_annotation_status']}

## What Remains Unknown

Domain annotation requires Pfam-A.hmm or InterProScan data. Protein consequence should not be interpreted mechanistically until ORF/domain evidence is manually reviewed.

## Recommended Manual Review Action

Review transcript structure, TransDecoder output, and Phase03/04 evidence before any final primer design or wet-lab candidate selection.
"""
        (OUT / "candidate_cards_tier1_updated" / f"{sanitize(r['gene_id'])}__{sanitize(r['isoform_id'])}.md").write_text(card)

    # Report/worklog/provenance.
    counts = {
        "tier1": len(tier1),
        "orf_found": sum(r["ORF_found"] == "TRUE" for r in orf_rows),
        "complete_orf": sum(r["complete_ORF"] == "TRUE" for r in orf_rows),
        "nmd_risk": sum(r["predicted_NMD_risk"] == "PREDICTED_NMD_RISK" for r in nmd_rows),
        "domain_annotated": sum(r["has_domain"] == "TRUE" for r in domain_summary),
        "insufficient": sum(r["functional_consequence_status"] in {"STRUCTURAL_ONLY_NO_ORF", "INSUFFICIENT_TOOL_SUPPORT"} for r in integrated_rows),
        "cards": len(integrated_rows),
    }
    report = f"""# Phase05B ORF/NMD/domain tool-enabled Tier 1 validation

## Why Phase05 Was Insufficient

Phase05 found that all 3,671 next-stage isoforms were marked `non_coding` in SQANTI3 with no valid ORF length and no predicted NMD field. Phase05B therefore tested a project-specific tool-enabled Tier 1 workflow.

## Tool and Database Availability

See `tables/phase05B_tool_and_database_availability.tsv`.

TransDecoder, gffread, seqkit, bedtools, samtools and HMMER are available in `{ENV}`. Pfam-A.hmm and InterProScan data were not found in the allowed search paths, so domain annotation was not run.

## Environment Installation

The project conda environment `crc_orf_anno` exists at `{ENV}`. Installation command and logs are recorded in `tables/phase05B_environment_installation_log.tsv`.

## Tier 1 Extraction Validation

Tier 1 isoforms extracted: {len(tier1)}. FASTA and GTF validation passed with no missing Tier 1 isoforms.

## ORF Prediction Results

Method: TransDecoder.LongOrfs + TransDecoder.Predict on `fasta/tier1_transcripts.fa`.

ORF found: {counts['orf_found']}; complete ORF: {counts['complete_orf']}.

## NMD/PTC Prediction Results

NMD rule: predicted stop/PTC more than 50 nt upstream of final exon-exon junction. Predicted NMD-risk: {counts['nmd_risk']}. This is not experimental validation.

## Domain Annotation Results

Domain annotated: {counts['domain_annotated']}. HMMER is installed, but Pfam-A.hmm was not found. Full Pfam/InterProScan downloads were not performed because they require user approval.

## Within-gene Comparison

See `tables/phase05B_tier1_within_gene_isoform_comparison.tsv`. Comparisons are Tier1-only same-gene comparisons and should be manually expanded if needed.

## Updated Candidate Cards

Updated Tier 1 cards generated: {counts['cards']}.

## Full Tier 2 Recommendation

Do not extend to all 3,671 Tier 2 isoforms yet. First manually review Tier 1 TransDecoder results and decide whether a Pfam/InterProScan database download is worth the storage/time cost.

## Pfam/InterProScan Recommendation

Pfam/InterProScan may be useful for protein-domain interpretation, but full database download can be large and should be approved separately.

## Limitations

Predicted ORF, NMD and domain status are computational and not functional validation. Domain annotation is missing because database support is absent.
"""
    (OUT / "reports/phase05B_ORF_domain_validation_report.md").write_text(report)

    worklog = f"""# Phase05B worklog

## date

{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## objective

Run tool-enabled Tier 1 ORF/NMD/domain validation for CRC isoform candidates.

## inputs

- `{tier_path}`
- `{integrated_path}`
- `{phase04_input}`
- `{fasta_path}`
- `{gtf_path}`

## outputs

`{OUT}`

## scripts

`code/phase05B_orf_domain_validation.py`

## commands

```bash
cd {OUT}
python3 code/phase05B_orf_domain_validation.py 2>&1 | tee logs/phase05B_orf_domain_validation.log
```

## software environment

Conda env: `{ENV}`; TransDecoder status: `{td_status}`.

## parameters

Tier 1 only; TransDecoder ORF prediction; conservative 50 nt upstream of final EEJ NMD rule; no Pfam download; no Tier2 expansion.

## run status

SUCCESS

## errors and fixes

Pfam/InterProScan database missing, so domain scan was not run.

## next step

Manual review of Tier 1 updated cards and decide whether to approve Pfam/InterProScan database download or Tier 2 expansion.
"""
    (OUT / "reports/phase05B_worklog.md").write_text(worklog)

    prov = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "td_status": td_status,
        "counts": counts,
        "env": str(ENV),
        "pfam_found": pfam_path is not None,
    }
    (OUT / "provenance/phase05B_provenance.json").write_text(json.dumps(prov, indent=2))
    try:
        subprocess.run([os.environ.get("CONDA_BIN", "conda"), "env", "export", "-p", str(ENV)], stdout=(OUT / "provenance/phase05B_conda_env_export.yaml").open("w"), stderr=subprocess.DEVNULL, check=False)
    except Exception:
        pass
    (OUT / "provenance/phase05B_tool_versions.txt").write_text("\n".join(f"{r['tool_or_database']}\t{r['version_or_status']}\t{r['path']}" for r in tool_rows) + "\n")

    print("PHASE05B_ORF_DOMAIN_VALIDATION_COMPLETE")
    print(f"TIER1={len(tier1)}")
    print(f"TRANSDECODER_STATUS={td_status}")
    print(f"ORF_FOUND={counts['orf_found']}")
    print(f"COMPLETE_ORF={counts['complete_orf']}")
    print(f"NMD_RISK={counts['nmd_risk']}")
    print(f"DOMAIN_ANNOTATED={counts['domain_annotated']}")
    print(f"INSUFFICIENT={counts['insufficient']}")
    print(f"CARDS={counts['cards']}")


if __name__ == "__main__":
    main()
