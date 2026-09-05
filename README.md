# CRC paired N/A/T long-read isoform remodeling

This repository contains the analysis code accompanying the manuscript on
patient-matched normal mucosa (N), tumor-adjacent tissue (A), and colorectal
tumor tissue (T) transcriptomes.

## Scope

The code is organized by the analysis layer described in the manuscript:

- `workflow/`: short-read STAR alignment and long-read Pychopper/FLAIR/SQANTI3
  reconstruction and annotation workflows.
- `analysis/`: cohort locking, catalog construction, DTE, satuRn-stageR DTU,
  N/A/T tissue-state classification, molecular annotation, and cross-platform
  gene/isoform comparisons.
- `figures/`: final computational Figure 1-3 generation scripts.
- `config/`: example path configuration. No study data are included.
- `provenance/`: method-to-script mapping and release QA.

## Data and privacy boundary

This code release contains no FASTQ, BAM, CRAM, count/usage matrix, clinical
spreadsheet, patient-level metadata, phone number, personal name, credential,
token, or server key. De-identified sample identifiers can appear only where
they are part of an analysis rule or a documented technical exception.

Patient-level clinical processing and wet-laboratory figure assembly are
intentionally excluded from this public computational release. They require
separate privacy and primary-source audits and are not part of the STAR/FLAIR
analysis workflow released here.

Raw long- and short-read sequencing data are associated with GSA-Human
accession `HRA018980`; access is governed by the repository record. Processed
matrices and Source Data are released separately from this code repository.

## Reproducibility boundary

The scripts preserve the statistical definitions used in the manuscript.
They do not rerun or alter frozen results during repository preparation.
Site-specific absolute paths are replaced by explicit environment variables in
the public copies. The original-source hashes and every release-layer change
are recorded in `provenance/`.

Copy `config/paths.example.env`, set the paths for the local installation, and
source it before running a script. `CRC_ROOT`, `CRC_STUDY_ROOT`, and
`CRC_V2_ROOT` locate the frozen project inputs. `CRC_FIG3_ROBUSTNESS_ROOT` is
needed only for the Figure 3 supplementary-table renderer, and `CRC_ORF_ENV`
is needed only for the ORF/domain annotation workflow.

The upstream files in `workflow/` are portable publication wrappers assembled
from the retained final commands and output provenance. They preserve the
analysis parameters but omit server scheduling, absolute paths, and restart
bookkeeping. They are not presented as byte-identical copies of every original
site-specific launcher. The FLAIR align interface is reconstructed from the
retained manuscript method and the versioned FLAIR interface because its exact
site launcher was not recovered. See `provenance/WORKFLOW_PROVENANCE.md`.

Basecalling software, model, and minimum read-quality settings were not
recoverable from the retained analysis records at the time of this release and
are therefore not invented here. The repository starts from the retained
Pychopper input stage unless those records are restored.

## Analysis overview

1. Match short reads with STAR splice-junction evidence.
2. Process long reads and reconstruct transcripts with FLAIR.
3. Combine and quantify the cohort-specific transcript atlas.
4. Classify and curate transcript structures with SQANTI3.
5. Freeze the N/A/T cohort and build the master isoform catalog.
6. Run paired DTE and satuRn-stageR DTU analyses.
7. Classify cross-sectional N/A/T tissue-state patterns.
8. Audit coding, NMD, ORF, protein-length, and Pfam consequences.
9. Compare long-read isoform remodeling with short-read gene abundance.
10. Render the final computational Figures 1-3 from frozen Source Data.

See `provenance/METHODS_TO_SCRIPT_MAP.tsv` for the exact code mapping.
