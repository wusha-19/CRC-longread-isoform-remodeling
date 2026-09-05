#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


PROJECT = Path(os.environ["CRC_STUDY_ROOT"]).expanduser().resolve()
OUT = PROJECT / "05E_comparison_isoform_domain"
PHASE05D = PROJECT / "05D_Pfam_domain_annotation"
TIER1_UPDATE = PHASE05D / "tables/phase05D_tier1_integrated_domain_update.tsv"
TIER1_PFAM_SUMMARY = PHASE05D / "tables/phase05D_tier1_pfam_domain_summary.tsv"
WITHIN_GENE = PROJECT / "05B_ORF_domain_validation/tables/phase05B_tier1_within_gene_isoform_comparison.tsv"
PRIMARY_GTF = PROJECT / "01_master_catalog/annotation/primary72_isoforms.gtf"
PRIMARY_FASTA = PROJECT / "01_master_catalog/annotation/primary72_isoforms.fa"
MASTER = PROJECT / "01_master_catalog/tables/CRC_isoform_master_catalog.tsv"
PFAM = PHASE05D / "db/Pfam-A.hmm"
DOMTBLOUT = OUT / "domains/comparison_vs_pfam.domtblout"

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def ensure_dirs() -> None:
    for sub in ["code", "tables", "fasta", "gtf", "proteins", "domains", "figures", "reports", "logs", "provenance", "candidate_cards_domain_comparison"]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "NA") for k in fields})


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    header = None
    seq: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    out[header] = "".join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        out[header] = "".join(seq)
    return out


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w") as fh:
        for header, seq in records.items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def valid(value: str | None) -> bool:
    return value not in {None, "", "NA", "nan", "None"}


def to_int(value: str, default: int = -1) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTUacgtu", "TGCAAtgcaa"))[::-1].upper().replace("U", "T")


def translate(nt: str) -> str:
    nt = nt.upper().replace("U", "T")
    return "".join(GENETIC_CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def longest_simple_orf(seq: str) -> dict[str, object] | None:
    """Same fallback as Phase05B: longest ATG-to-stop ORF across six frames, minimum 90 nt."""
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
                                    start, end = len(seq) - (j + 3) + 1, len(seq) - i
                                nt = s[i:j + 3]
                                cand = {
                                    "strand": strand,
                                    "frame": frame,
                                    "ORF_start": start,
                                    "ORF_end": end,
                                    "ORF_length_nt": nt_len,
                                    "protein_sequence_raw": translate(nt),
                                    "longest_ORF_sequence_id": f"simple_orf_{strand}{frame}_{start}_{end}",
                                }
                                if best is None or nt_len > int(best["ORF_length_nt"]):
                                    best = cand
                            break
                        j += 3
                    i = j
                i += 3
    return best


def parse_attrs(line: str) -> dict[str, str]:
    attrs = {}
    for key, value in re.findall(r'([A-Za-z0-9_]+) "([^"]+)"', line):
        attrs[key] = value
    return attrs


def select_comparisons() -> list[dict[str, object]]:
    tier1 = read_tsv(TIER1_UPDATE)
    within = {r["tier1_isoform_id"]: r for r in read_tsv(WITHIN_GENE)}
    master_rows = read_tsv(MASTER)
    fasta_ids = set(read_fasta(PRIMARY_FASTA))
    tier1_ids = {r["isoform_id"] for r in tier1}
    by_gene = defaultdict(list)
    for r in master_rows:
        if r.get("analysis_set_primary72") == "True" or r.get("analysis_set_primary72") == "TRUE":
            by_gene[r["gene_id"]].append(r)

    selected = []
    for r in tier1:
        iso = r["isoform_id"]
        gene = r["gene_id"]
        comp = within.get(iso, {}).get("comparison_isoform_id", "NA")
        rule = "NO_COMPARISON_ISOFORM_AVAILABLE"
        reason = "No same-gene non-candidate isoform was available in primary72 FASTA."
        notes = "NA"
        if valid(comp) and comp != iso:
            rule = "PRIORITY1_PHASE05B_WITHIN_GENE_TABLE"
            reason = "Explicit comparison_isoform_id provided by Phase05B within-gene comparison table."
        else:
            candidates = [x for x in by_gene.get(gene, []) if x["isoform_id"] != iso and x["isoform_id"] not in tier1_ids and x["isoform_id"] in fasta_ids]
            candidates.sort(key=lambda x: (to_int(x.get("strict_count_total_primary72", "0")), to_int(x.get("support_ge2_primary72", "0")), to_int(x.get("length", "0"))), reverse=True)
            if candidates:
                comp = candidates[0]["isoform_id"]
                if str(r.get("dominant_isoform_switch", "")).upper() == "TRUE":
                    rule = "PRIORITY3_DOMINANT_SWITCH_PROXY_PRIMARY72_SUPPORT"
                    reason = "Candidate marked dominant_isoform_switch; explicit switch partner unavailable, selected highest primary72 supported non-candidate same-gene isoform."
                else:
                    rule = "PRIORITY2_HIGHEST_PRIMARY72_SUPPORT_PROXY"
                    reason = "No explicit comparison_isoform_id; selected highest primary72 strict_count_total non-candidate same-gene isoform with transcript FASTA."
                notes = f"strict_count_total_primary72={candidates[0].get('strict_count_total_primary72', 'NA')}; support_ge2_primary72={candidates[0].get('support_ge2_primary72', 'NA')}"
            else:
                comp = "NA"
        selected.append({
            "gene_id": gene,
            "gene_name": r.get("gene_name", "NA"),
            "tier1_isoform_id": iso,
            "comparison_isoform_id": comp,
            "selection_rule": rule,
            "selection_reason": reason,
            "comparison_available": str(valid(comp) and comp in fasta_ids).upper(),
            "notes": notes,
        })
    fields = ["gene_id", "gene_name", "tier1_isoform_id", "comparison_isoform_id", "selection_rule", "selection_reason", "comparison_available", "notes"]
    write_tsv(OUT / "tables/phase05E_comparison_isoform_selection.tsv", selected, fields)
    return selected


def extract_fasta_gtf(selection_rows: list[dict[str, object]]) -> dict[str, str]:
    fasta = read_fasta(PRIMARY_FASTA)
    requested = {str(r["comparison_isoform_id"]) for r in selection_rows if r["comparison_available"] == "TRUE"}
    extracted = {iso: fasta[iso] for iso in sorted(requested) if iso in fasta}
    write_fasta(OUT / "fasta/comparison_isoforms.fa", extracted)

    gtf_ids = set()
    with PRIMARY_GTF.open() as inp, (OUT / "gtf/comparison_isoforms.gtf").open("w") as out:
        for line in inp:
            attrs = parse_attrs(line)
            tid = attrs.get("transcript_id")
            if tid in requested:
                out.write(line)
                gtf_ids.add(tid)
    rows = []
    for iso in sorted(requested):
        rows.append({
            "comparison_isoform_id": iso,
            "present_in_primary72_fasta": str(iso in fasta).upper(),
            "present_in_primary72_gtf": str(iso in gtf_ids).upper(),
            "extracted": str(iso in extracted and iso in gtf_ids).upper(),
            "notes": "OK" if iso in extracted and iso in gtf_ids else "Missing FASTA or GTF; downstream ORF/Pfam skipped if FASTA absent.",
        })
    write_tsv(OUT / "tables/comparison_isoform_extraction_status.tsv", rows, ["comparison_isoform_id", "present_in_primary72_fasta", "present_in_primary72_gtf", "extracted", "notes"])
    return extracted


def predict_orfs(comparison_fasta: dict[str, str], selection_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    gene_for = {}
    name_for = {}
    for r in selection_rows:
        if r["comparison_available"] == "TRUE":
            gene_for[str(r["comparison_isoform_id"])] = str(r["gene_id"])
            name_for[str(r["comparison_isoform_id"])] = str(r["gene_name"])
    rows = []
    proteins = {}
    clean = {}
    for iso, seq in comparison_fasta.items():
        hit = longest_simple_orf(seq)
        if hit is None:
            rows.append({
                "comparison_isoform_id": iso, "gene_id": gene_for.get(iso, "NA"), "gene_name": name_for.get(iso, "NA"),
                "ORF_found": "FALSE", "ORF_source": "SIMPLE_ORF_SCREEN_ONLY", "ORF_start": "NA", "ORF_end": "NA",
                "ORF_length_nt": "NA", "protein_length_aa": 0, "complete_ORF": "FALSE", "has_internal_stop": "FALSE",
                "n_internal_stop": 0, "ORF_confidence": "INSUFFICIENT_NO_SIMPLE_ORF", "notes": "No ATG-to-stop ORF >=90 nt found by Phase05B simple ORF fallback.",
            })
            continue
        prot = str(hit["protein_sequence_raw"])
        terminal_stop_removed = prot.endswith("*")
        if terminal_stop_removed:
            prot = prot[:-1]
        n_internal = prot.count("*")
        clean_prot = prot.replace("*", "X")
        header = f"{iso}|gene_id={gene_for.get(iso, 'NA')}|ORF_start={hit['ORF_start']}|ORF_end={hit['ORF_end']}|protein_length_aa={len(clean_prot)}|ORF_source=SIMPLE_ORF_SCREEN_ONLY"
        proteins[header] = prot
        clean[header] = clean_prot
        has_internal = n_internal > 0
        rows.append({
            "comparison_isoform_id": iso,
            "gene_id": gene_for.get(iso, "NA"),
            "gene_name": name_for.get(iso, "NA"),
            "ORF_found": "TRUE",
            "ORF_source": "SIMPLE_ORF_SCREEN_ONLY",
            "ORF_start": hit["ORF_start"],
            "ORF_end": hit["ORF_end"],
            "ORF_length_nt": hit["ORF_length_nt"],
            "protein_length_aa": len(clean_prot),
            "complete_ORF": "TRUE",
            "has_internal_stop": str(has_internal).upper(),
            "n_internal_stop": n_internal,
            "ORF_confidence": "LOW_SIMPLE_ORF_INTERNAL_STOP" if has_internal else "MODERATE_SIMPLE_ORF_NO_INTERNAL_STOP",
            "notes": f"strand={hit['strand']};frame={hit['frame']};terminal_stop_removed={str(terminal_stop_removed).upper()}",
        })
    fields = ["comparison_isoform_id", "gene_id", "gene_name", "ORF_found", "ORF_source", "ORF_start", "ORF_end", "ORF_length_nt", "protein_length_aa", "complete_ORF", "has_internal_stop", "n_internal_stop", "ORF_confidence", "notes"]
    write_tsv(OUT / "tables/phase05E_comparison_ORF_prediction.tsv", rows, fields)
    write_fasta(OUT / "proteins/comparison_predicted_proteins.fa", proteins)
    write_fasta(OUT / "proteins/comparison_predicted_proteins.hmmscan_clean.fa", clean)
    return rows


def run_hmmscan_if_needed() -> None:
    protein_fa = OUT / "proteins/comparison_predicted_proteins.hmmscan_clean.fa"
    n = 0
    if protein_fa.exists():
        with protein_fa.open() as fh:
            n = sum(1 for line in fh if line.startswith(">"))
    if n == 0:
        (OUT / "domains/comparison_vs_pfam.hmmscan.txt").write_text("NO_COMPARISON_PROTEINS\n")
        DOMTBLOUT.write_text("")
        return
    cmd = [
        "hmmscan",
        "--cpu", "4",
        "--domtblout", str(DOMTBLOUT),
        str(PFAM),
        str(protein_fa),
    ]
    with (OUT / "domains/comparison_vs_pfam.hmmscan.txt").open("w") as out:
        subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.STDOUT)


def parse_domtblout() -> list[dict[str, object]]:
    orf = {r["comparison_isoform_id"]: r for r in read_tsv(OUT / "tables/phase05E_comparison_ORF_prediction.tsv")}
    rows = []
    if not DOMTBLOUT.exists():
        return rows
    with DOMTBLOUT.open() as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            parts = raw.rstrip("\n").split(maxsplit=22)
            if len(parts) < 23:
                continue
            pfam_name, pfam_acc, query = parts[0], parts[1], parts[3]
            full_e = float(parts[6])
            i_eval = float(parts[12])
            hit_pass = i_eval <= 1e-5 or full_e <= 1e-5
            iso = query.split("|", 1)[0]
            st = orf.get(iso, {})
            rows.append({
                "comparison_isoform_id": iso,
                "protein_id": query,
                "gene_id": st.get("gene_id", "NA"),
                "gene_name": st.get("gene_name", "NA"),
                "pfam_acc": pfam_acc,
                "pfam_name": pfam_name,
                "domain_i_evalue": i_eval,
                "domain_score": parts[13],
                "hmm_start": parts[15],
                "hmm_end": parts[16],
                "ali_start": parts[17],
                "ali_end": parts[18],
                "env_start": parts[19],
                "env_end": parts[20],
                "description": parts[22],
                "hit_pass": str(hit_pass).upper(),
                "has_internal_stop": st.get("has_internal_stop", "NA"),
                "ORF_confidence": st.get("ORF_confidence", "NA"),
                "domain_interpretation_confidence": "LOW_INTERNAL_STOP_WARNING" if st.get("has_internal_stop") == "TRUE" else "HIGH_NO_INTERNAL_STOP",
                "full_sequence_evalue": full_e,
            })
    fields = ["comparison_isoform_id", "protein_id", "gene_id", "gene_name", "pfam_acc", "pfam_name", "domain_i_evalue", "domain_score", "hmm_start", "hmm_end", "ali_start", "ali_end", "env_start", "env_end", "description", "hit_pass", "has_internal_stop", "ORF_confidence", "domain_interpretation_confidence"]
    write_tsv(OUT / "tables/phase05E_comparison_pfam_domain_hits.tsv", rows, fields)
    return rows


def domain_set(domain_list: str) -> set[str]:
    if not valid(domain_list):
        return set()
    return {x.strip() for x in domain_list.split(";") if x.strip() and x.strip() != "NA"}


def summarize_and_assess(selection_rows: list[dict[str, object]]) -> None:
    orf = {r["comparison_isoform_id"]: r for r in read_tsv(OUT / "tables/phase05E_comparison_ORF_prediction.tsv")}
    hits = parse_domtblout()
    pass_hits = [h for h in hits if h["hit_pass"] == "TRUE"]
    by_iso = defaultdict(list)
    for h in pass_hits:
        by_iso[h["comparison_isoform_id"]].append(h)
    comp_summary = []
    for iso, st in orf.items():
        phits = sorted(by_iso.get(iso, []), key=lambda h: (float(h["domain_i_evalue"]), -float(h["domain_score"])))
        dlist = ";".join(dict.fromkeys(f"{h['pfam_name']}({h['pfam_acc']})" for h in phits)) if phits else "NA"
        has_internal = st["has_internal_stop"]
        if has_internal == "TRUE":
            status = "LOW_CONFIDENCE_ORF_INTERNAL_STOP"
            conf = "LOW_INTERNAL_STOP_WARNING"
        elif phits:
            status = "PFAM_DOMAIN_FOUND"
            conf = "HIGH_NO_INTERNAL_STOP"
        else:
            status = "NO_PFAM_DOMAIN_FOUND"
            conf = "INSUFFICIENT"
        comp_summary.append({
            "comparison_isoform_id": iso,
            "gene_id": st["gene_id"],
            "gene_name": st["gene_name"],
            "protein_length_aa": st["protein_length_aa"],
            "has_internal_stop": has_internal,
            "ORF_confidence": st["ORF_confidence"],
            "has_pfam_domain": str(bool(phits)).upper(),
            "n_pfam_domains": len(phits),
            "pfam_domain_list": dlist,
            "best_domain_i_evalue": min((float(h["domain_i_evalue"]) for h in phits), default="NA"),
            "domain_annotation_status": status,
            "domain_interpretation_confidence": conf,
        })
    write_tsv(OUT / "tables/phase05E_comparison_pfam_domain_summary.tsv", comp_summary, ["comparison_isoform_id", "gene_id", "gene_name", "protein_length_aa", "has_internal_stop", "ORF_confidence", "has_pfam_domain", "n_pfam_domains", "pfam_domain_list", "best_domain_i_evalue", "domain_annotation_status", "domain_interpretation_confidence"])

    tier1_summary = {r["isoform_id"]: r for r in read_tsv(TIER1_PFAM_SUMMARY)}
    tier1_update = {r["isoform_id"]: r for r in read_tsv(TIER1_UPDATE)}
    comp_by_iso = {r["comparison_isoform_id"]: r for r in comp_summary}
    assess_rows = []
    integrated = []
    for sel in selection_rows:
        tiso = str(sel["tier1_isoform_id"])
        ciso = str(sel["comparison_isoform_id"])
        t = tier1_summary.get(tiso, {})
        tu = tier1_update.get(tiso, {})
        c = comp_by_iso.get(ciso, {})
        if sel["comparison_available"] != "TRUE" or not c:
            status = "INSUFFICIENT_COMPARISON_DATA"
            conf = "INSUFFICIENT"
            gain = loss = shared = "NA"
            notes = "No usable comparison isoform protein/Pfam annotation."
        elif t.get("has_internal_stop") == "TRUE" or c.get("has_internal_stop") == "TRUE":
            status = "LOW_CONFIDENCE_INTERNAL_STOP"
            conf = "LOW"
            td = domain_set(t.get("pfam_domain_list", "NA"))
            cd = domain_set(c.get("pfam_domain_list", "NA"))
            gain = ";".join(sorted(td - cd)) if td - cd else "NA"
            loss = ";".join(sorted(cd - td)) if cd - td else "NA"
            shared = ";".join(sorted(td & cd)) if td & cd else "NA"
            notes = "Either candidate or comparison protein has internal stop; domain differences are low-confidence hints only."
        else:
            td = domain_set(t.get("pfam_domain_list", "NA"))
            cd = domain_set(c.get("pfam_domain_list", "NA"))
            gain_set = td - cd
            loss_set = cd - td
            shared_set = td & cd
            gain = ";".join(sorted(gain_set)) if gain_set else "NA"
            loss = ";".join(sorted(loss_set)) if loss_set else "NA"
            shared = ";".join(sorted(shared_set)) if shared_set else "NA"
            if not td and not cd:
                status = "NO_DOMAIN_IN_EITHER"
            elif gain_set and loss_set:
                status = "DOMAIN_GAIN_AND_LOSS"
            elif gain_set:
                status = "DOMAIN_GAIN"
            elif loss_set:
                status = "DOMAIN_LOSS"
            else:
                status = "DOMAIN_RETAINED"
            conf = "HIGH"
            notes = "Comparable simple ORF/Pfam annotation without internal stop warnings; computational annotation only."
        assess_rows.append({
            "gene_id": sel["gene_id"],
            "gene_name": sel["gene_name"],
            "tier1_isoform_id": tiso,
            "comparison_isoform_id": ciso,
            "tier1_domain_list": t.get("pfam_domain_list", "NA"),
            "comparison_domain_list": c.get("pfam_domain_list", "NA"),
            "tier1_has_internal_stop": t.get("has_internal_stop", "NA"),
            "comparison_has_internal_stop": c.get("has_internal_stop", "NA"),
            "tier1_ORF_confidence": t.get("ORF_confidence", "NA"),
            "comparison_ORF_confidence": c.get("ORF_confidence", "NA"),
            "domain_gain_loss_status": status,
            "domain_gain_list": gain,
            "domain_loss_list": loss,
            "domain_shared_list": shared,
            "assessment_confidence": conf,
            "notes": notes,
        })
        integrated.append({
            "gene_id": sel["gene_id"],
            "gene_name": sel["gene_name"],
            "tier1_isoform_id": tiso,
            "comparison_isoform_id": ciso,
            "review_priority": tu.get("review_priority", "NA"),
            "in_DTE": tu.get("in_DTE", "NA"),
            "in_DTU": tu.get("in_DTU", "NA"),
            "gene_silent_DTU": tu.get("gene_silent_DTU", "NA"),
            "dominant_isoform_switch": tu.get("dominant_isoform_switch", "NA"),
            "expression_dynamic_pattern": tu.get("expression_dynamic_pattern", "NA"),
            "usage_dynamic_pattern": tu.get("usage_dynamic_pattern", "NA"),
            "tier1_has_pfam_domain": t.get("has_pfam_domain", "NA"),
            "comparison_has_pfam_domain": c.get("has_pfam_domain", "NA"),
            "tier1_domain_list": t.get("pfam_domain_list", "NA"),
            "comparison_domain_list": c.get("pfam_domain_list", "NA"),
            "domain_gain_loss_status": status,
            "assessment_confidence": conf,
            "tier1_has_internal_stop": t.get("has_internal_stop", "NA"),
            "comparison_has_internal_stop": c.get("has_internal_stop", "NA"),
            "RT_PCR_readiness": tu.get("RT_PCR_readiness", "NA"),
            "manual_review_note": tu.get("manual_review_note", "NA"),
        })
    assess_fields = ["gene_id", "gene_name", "tier1_isoform_id", "comparison_isoform_id", "tier1_domain_list", "comparison_domain_list", "tier1_has_internal_stop", "comparison_has_internal_stop", "tier1_ORF_confidence", "comparison_ORF_confidence", "domain_gain_loss_status", "domain_gain_list", "domain_loss_list", "domain_shared_list", "assessment_confidence", "notes"]
    write_tsv(OUT / "tables/phase05E_domain_gain_loss_assessment.tsv", assess_rows, assess_fields)
    int_fields = ["gene_id", "gene_name", "tier1_isoform_id", "comparison_isoform_id", "review_priority", "in_DTE", "in_DTU", "gene_silent_DTU", "dominant_isoform_switch", "expression_dynamic_pattern", "usage_dynamic_pattern", "tier1_has_pfam_domain", "comparison_has_pfam_domain", "tier1_domain_list", "comparison_domain_list", "domain_gain_loss_status", "assessment_confidence", "tier1_has_internal_stop", "comparison_has_internal_stop", "RT_PCR_readiness", "manual_review_note"]
    write_tsv(OUT / "tables/phase05E_integrated_tier1_domain_comparison.tsv", integrated, int_fields)
    make_cards(integrated, assess_rows)
    make_summary_figures(assess_rows, comp_summary)
    make_reports(selection_rows, orf, comp_summary, assess_rows)


def make_cards(integrated: list[dict[str, object]], assess_rows: list[dict[str, object]]) -> None:
    outdir = OUT / "candidate_cards_domain_comparison"
    for old in outdir.glob("*.md"):
        old.unlink()
    assess = {r["tier1_isoform_id"]: r for r in assess_rows}
    for row in integrated:
        a = assess[row["tier1_isoform_id"]]
        safe = str(row["tier1_isoform_id"]).replace("/", "_").replace("|", "_")
        text = f"""# Phase05E domain comparison card: {row['tier1_isoform_id']}

## Candidate domain list

{row['tier1_domain_list']}

## Comparison isoform

- comparison_isoform_id: {row['comparison_isoform_id']}
- comparison_domain_list: {row['comparison_domain_list']}

## Domain gain/loss assessment

- status: {row['domain_gain_loss_status']}
- confidence: {row['assessment_confidence']}
- gain: {a['domain_gain_list']}
- loss: {a['domain_loss_list']}
- shared: {a['domain_shared_list']}

## Internal stop caveat

- tier1_has_internal_stop: {row['tier1_has_internal_stop']}
- comparison_has_internal_stop: {row['comparison_has_internal_stop']}

## Interpretation boundary

This is a computational simple ORF/Pfam comparison. It is not experimental mechanism evidence and must not be used as final wet-lab candidate selection by itself.
"""
        (outdir / f"{safe}.phase05E_domain_comparison_card.md").write_text(text)


def make_summary_figures(assess_rows: list[dict[str, object]], comp_summary: list[dict[str, object]]) -> None:
    status_counts = Counter(r["domain_gain_loss_status"] for r in assess_rows)
    conf_counts = Counter(r["assessment_confidence"] for r in assess_rows)
    comp_domain_counts = Counter("with_pfam_domain" if r["has_pfam_domain"] == "TRUE" else "without_pfam_domain" for r in comp_summary)
    summary_rows = [{"metric": "tier1_isoforms_processed", "value": len(assess_rows)}]
    summary_rows += [{"metric": f"domain_status_{k}", "value": v} for k, v in sorted(status_counts.items())]
    summary_rows += [{"metric": f"assessment_confidence_{k}", "value": v} for k, v in sorted(conf_counts.items())]
    summary_rows += [{"metric": f"comparison_{k}", "value": v} for k, v in sorted(comp_domain_counts.items())]
    write_tsv(OUT / "tables/phase05E_domain_comparison_summary.tsv", summary_rows, ["metric", "value"])
    write_tsv(OUT / "figures/domain_gain_loss_status_summary.data.tsv", [{"category": k, "count": v} for k, v in sorted(status_counts.items())], ["category", "count"])
    write_tsv(OUT / "figures/domain_assessment_confidence_summary.data.tsv", [{"category": k, "count": v} for k, v in sorted(conf_counts.items())], ["category", "count"])
    r_code = r'''
plot_bar <- function(infile, outfile, main) {
  d <- read.delim(infile, check.names=FALSE)
  pdf(outfile, width=7.5, height=4.8)
  op <- par(mar=c(8,4,3,1))
  barplot(d$count, names.arg=d$category, las=2, col="#4C78A8", ylab="count", main=main)
  par(op)
  dev.off()
}
plot_bar("figures/domain_gain_loss_status_summary.data.tsv", "figures/domain_gain_loss_status_summary.pdf", "Phase05E domain gain/loss status")
plot_bar("figures/domain_assessment_confidence_summary.data.tsv", "figures/domain_assessment_confidence_summary.pdf", "Phase05E assessment confidence")
'''
    (OUT / "code/phase05E_make_figures.R").write_text(r_code)
    subprocess.run(["Rscript", "code/phase05E_make_figures.R"], cwd=OUT, check=True)


def make_reports(selection_rows, orf, comp_summary, assess_rows) -> None:
    selected = sum(1 for r in selection_rows if r["comparison_available"] == "TRUE")
    proteins = sum(1 for r in orf.values() if r["ORF_found"] == "TRUE")
    with_domains = sum(1 for r in comp_summary if r["has_pfam_domain"] == "TRUE")
    status_counts = Counter(r["domain_gain_loss_status"] for r in assess_rows)
    high_assessable = sum(1 for r in assess_rows if r["assessment_confidence"] == "HIGH" and r["domain_gain_loss_status"] not in {"INSUFFICIENT_COMPARISON_DATA", "LOW_CONFIDENCE_INTERNAL_STOP"})
    report = f"""# Phase05E same-gene comparison isoform ORF/Pfam annotation

## Why Phase05E was needed

Phase05D annotated Pfam domains only for Tier 1 candidate proteins, so domain gain/loss could not be assessed. Phase05E adds same-gene comparison isoform simple ORF/Pfam annotation with the same conservative fallback used in Phase05B.

## Comparison isoform selection method

Priority 1 used explicit `comparison_isoform_id` from Phase05B. If absent, the pipeline selected the highest primary72-supported non-candidate same-gene isoform with transcript FASTA, using `strict_count_total_primary72` as a conservative support proxy.

## Counts

- Tier 1 isoforms processed: {len(selection_rows)}
- Comparison isoforms selected: {selected}
- Comparison proteins generated: {proteins}
- Comparison proteins with Pfam domains: {with_domains}
- High-confidence assessable comparisons: {high_assessable}

## Domain gain/loss assessment summary

{chr(10).join(f'- {k}: {v}' for k, v in sorted(status_counts.items()))}

## Internal stop caveat

If either candidate or comparison protein has an internal stop warning, the assessment is marked `LOW_CONFIDENCE_INTERNAL_STOP`; domain differences are hints only.

## Candidates with assessable domain consequence

See `tables/phase05E_domain_gain_loss_assessment.tsv`. High-confidence rows are those with no internal stop and usable comparison ORF/Pfam annotation.

## Limitations

ORF prediction remains `SIMPLE_ORF_SCREEN_ONLY`, not validated coding potential. Pfam hits are computational annotations, not experimental mechanism evidence. Fallback comparison selection uses primary72 support when explicit biological switch partners are unavailable.

## Recommendation for Phase05C revised prioritization

Use Phase05E as a conservative evidence layer for prioritization, combined with DTE/DTU/NAT evidence and RT-PCR feasibility. Do not make final wet-lab candidate selection solely from domain gain/loss status.
"""
    (OUT / "reports/phase05E_domain_comparison_report.md").write_text(report)
    worklog = f"""# Phase05E worklog

- Time: {now()}
- No DTE/DTU/SUPPA2 rerun, InterProScan, new database download, Tier2 expansion, primer design, or final wet-lab candidate selection was performed.
- Pfam database reused from Phase05D: `{PFAM}`
- hmmscan threshold: domain i-Evalue <= 1e-5 or full sequence E-value <= 1e-5
- Tier 1 isoforms processed: {len(selection_rows)}
- Comparison isoforms selected: {selected}
- Comparison proteins generated: {proteins}
- Comparison proteins with Pfam domains: {with_domains}
- High-confidence assessable comparisons: {high_assessable}
"""
    (OUT / "reports/phase05E_worklog.md").write_text(worklog)
    prov = {
        "time": now(),
        "output": str(OUT),
        "pfam": str(PFAM),
        "threshold": "domain i-Evalue <= 1e-5 or full sequence E-value <= 1e-5",
        "forbidden_steps_not_run": ["DTE", "DTU", "SUPPA2", "InterProScan", "new_database_download", "Tier2_expansion", "primer_design", "final_wet_lab_candidate_selection"],
    }
    (OUT / "provenance/phase05E_provenance.json").write_text(json.dumps(prov, indent=2))


def main() -> None:
    ensure_dirs()
    for path in [TIER1_UPDATE, TIER1_PFAM_SUMMARY, WITHIN_GENE, PRIMARY_GTF, PRIMARY_FASTA, MASTER, PFAM]:
        if not path.exists() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required input: {path}")
    selection = select_comparisons()
    comp_fasta = extract_fasta_gtf(selection)
    predict_orfs(comp_fasta, selection)
    run_hmmscan_if_needed()
    summarize_and_assess(selection)
    print("PHASE05E_COMPLETE")


if __name__ == "__main__":
    main()
