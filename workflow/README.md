# Upstream workflow

These scripts expose the retained final command lines as portable one-sample or
one-stage wrappers. Set the reference paths in `config/paths.example.env` before
running them.

## Order

1. `01_pychopper.sh INPUT_FASTQ OUTPUT_DIR SAMPLE_ID`
2. Long-read branch: `02_flair_align.sh FULL_LENGTH_FASTQ OUTPUT_DIR SAMPLE_ID`
3. Matched short-read branch: `02_star_alignment.sh R1_FASTQ_GZ R2_FASTQ_GZ OUTPUT_DIR SAMPLE_ID`
4. `03_flair_transcriptome.sh LONGREAD_BAM STAR_SJ OUTPUT_DIR SAMPLE_ID`
5. `04_flair_combine.sh COMBINE_MANIFEST OUTPUT_PREFIX`
6. `05_flair_quantify.sh READS_MANIFEST COMBINED_FASTA COMBINED_BED OUTPUT_PREFIX TEMP_DIR`
7. `06_filter_expression_support.py ...`
8. `07_sqanti3_annotation.sh FILTERED_GTF FILTERED_FASTA FILTERED_BED OUTPUT_DIR`

The FLAIR reads manifest has four tab-separated columns: sample identifier,
condition, batch identifier, and full-length read path. The combine manifest is
the FLAIR combine manifest used by the installed FLAIR version.

The retained catalog support rule is fixed at count >= 2 in at least 8 samples
and total count >= 20. The scripts do not include raw data or sample metadata.

The exact site-specific FLAIR align launcher was not recovered during the
release audit. `02_flair_align.sh` reconstructs that interface from the retained
manuscript method and the FLAIR 3.0.0b1 command contract; this provenance status
is intentionally distinct from the directly recovered downstream commands.
