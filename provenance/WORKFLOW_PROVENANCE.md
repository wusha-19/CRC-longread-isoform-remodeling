# Upstream workflow provenance

The files in `workflow/` are publication-facing wrappers. They retain the
scientific commands and fixed parameters found in the final server-side output
provenance, while replacing installation-specific paths and orchestration with
arguments and environment variables. They are not byte-identical copies of all
historical launchers.

| Stage | Retained source evidence | Public file | Preserved method |
|---|---|---|---|
| Pychopper | `script/02_pychopper.sh` | `workflow/01_pychopper.sh` | 20 threads; report, unclassified, rescued, and full-length outputs |
| FLAIR align | manuscript method and FLAIR 3.0.0b1 interface; exact site launcher not recovered | `workflow/02_flair_align.sh` | splice-aware FLAIR/minimap2 alignment to the declared genome; aligned BAM, BAM index, and BED outputs |
| STAR | accepted paired-end STAR BAM and junction provenance; recovery launcher audit | `workflow/02_star_alignment.sh` | 8 threads; gzipped paired reads; coordinate-sorted BAM; `SJ.out.tab` |
| FLAIR transcriptome | `script/Flair/03_2_flair_transcriptome.sh` | `workflow/03_flair_transcriptome.sh` | GRCh38.p14/GENCODE v48; STAR junction evidence; `-s 5`; `byregion`; `nosubset`; stringent splice checking |
| FLAIR combine | `script/Flair/03_3_flair_combine.sh` | `workflow/04_flair_combine.sh` | manifest-based cohort combine with `-c` |
| FLAIR quantify | `script/Flair/03_4_flair_quantify_nosubset.sh` | `workflow/05_flair_quantify.sh` | four-column reads manifest; TPM output; stringent splice checking |
| Expression support | `script/Flair/filter_S5_transcriptome_by_strict_expression_support.py` | `workflow/06_filter_expression_support.py` | count >= 2 in >= 8 samples and total count >= 20 |
| SQANTI3 | final expression-supported QC/filter/rescue/final-QC launcher | `workflow/07_sqanti3_annotation.sh` | default rules filter; automatic rescue; final QC; `--force_id_ignore`; `--skipORF` |

The portable STAR wrapper excludes sample-specific recovery checks used to
repair a small number of damaged BAM artifacts. Those checks concern recovery
provenance rather than the alignment method. The public FLAIR wrappers likewise
exclude scheduler settings, concurrent batch splitting, completion markers,
and restart bookkeeping.

Basecalling software, model, minimum read-quality settings, and the byte-exact
site-specific FLAIR align launcher were not recoverable from the retained
records. No values are inferred. The documented workflow therefore begins at
the Pychopper input stage, and the alignment wrapper is explicitly identified
as a reconstruction from the retained method and versioned software interface.

The retained runtime audit identified minimap2 2.30-r1287, whereas an older
manuscript draft named version 2.28. The discrepancy is exposed in
`environment/SOFTWARE_VERSIONS.tsv` and should be reconciled against an archived
environment export before the manuscript states an exact minimap2 version.
