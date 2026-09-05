#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(20260827)

required_packages <- c(
  "DESeq2", "ggplot2", "patchwork", "dplyr", "tidyr", "readr",
  "scales", "ragg", "systemfonts"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

v2_base <- Sys.getenv(
  "V2_BASE",
  Sys.getenv("CRC_V2_ROOT")
)
study_root <- dirname(v2_base)
framework_root <- file.path(study_root, "16_paired_NAT_multi_platform_framework_20260827")
phase03_root <- file.path(study_root, "03_DTU_paired", "tables")

dir_tables <- file.path(v2_base, "tables")
dir_source <- file.path(v2_base, "source_data")
dir_reports <- file.path(v2_base, "reports")
dir_logs <- file.path(v2_base, "logs")
dir_main_pdf <- file.path(v2_base, "figures", "main", "pdf")
dir_main_svg <- file.path(v2_base, "figures", "main", "svg")
dir_main_png <- file.path(v2_base, "figures", "main", "png")
dir_supp_pdf <- file.path(v2_base, "figures", "supplementary", "pdf")
dir_supp_svg <- file.path(v2_base, "figures", "supplementary", "svg")
dir_supp_png <- file.path(v2_base, "figures", "supplementary", "png")
for (d in c(
  dir_tables, dir_source, dir_reports, dir_logs,
  dir_main_pdf, dir_main_svg, dir_main_png,
  dir_supp_pdf, dir_supp_svg, dir_supp_png
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

progress_path <- file.path(dir_logs, "25_matched_shortread_Fig3_progress.txt")
record_progress <- function(stage, detail) {
  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\t", stage, "\t", detail, "\n",
    file = progress_path, append = TRUE, sep = ""
  )
  message("[", stage, "] ", detail)
}

stop_if <- function(condition, message_text) {
  if (isTRUE(condition)) {
    record_progress("FAIL", message_text)
    stop(message_text, call. = FALSE)
  }
}

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x) | x == "." | toupper(x) == "NA"] <- NA_character_
  x
}

strip_ensembl_version <- function(x) {
  sub("[.][0-9]+$", "", clean_chr(x))
}

sha256_file <- function(path) {
  stop_if(!file.exists(path), paste0("Missing file for SHA256: ", path))
  out <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  stop_if(length(out) == 0L || !grepl("^[0-9a-f]{64}", out[[1]]), paste0("SHA256 failed: ", path))
  sub(" .*", "", out[[1]])
}

write_tsv2 <- function(x, path) {
  readr::write_tsv(x, path, na = "NA")
  path
}

write_lines2 <- function(x, path) {
  writeLines(x, path, useBytes = TRUE)
  path
}

assert_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  stop_if(length(missing) > 0L, paste0(label, " missing columns: ", paste(missing, collapse = ", ")))
}

file_snapshot <- function(paths) {
  paths <- sort(unique(paths[file.exists(paths) & !dir.exists(paths)]))
  tibble(
    path = paths,
    size_bytes = file.info(paths)$size,
    mtime = format(file.info(paths)$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    sha256 = vapply(paths, sha256_file, character(1))
  )
}

FONT_FAMILY <- "Times New Roman"
font_match <- systemfonts::match_fonts(FONT_FAMILY)[1, , drop = FALSE]
stop_if(
  is.null(font_match$path) || !nzchar(font_match$path) || !file.exists(font_match$path),
  "Times New Roman is unavailable in the remote R graphics environment."
)

theme_manuscript <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY, colour = "#202020"),
      axis.line = element_line(linewidth = 0.28, colour = "#333333"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#333333"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.6),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.9, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#555555", margin = margin(b = 2)),
      plot.caption = element_text(size = base_size - 1.0, colour = "#555555", hjust = 0),
      plot.margin = margin(3, 3, 3, 3),
      panel.grid.major.y = element_line(linewidth = 0.16, colour = "#E8E8E8"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}
theme_set(theme_manuscript())
update_geom_defaults("text", list(family = FONT_FAMILY))
update_geom_defaults("label", list(family = FONT_FAMILY))

save_bundle <- function(plot, basename_no_ext, width_mm, height_mm, supplementary = FALSE, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  out_dirs <- if (supplementary) {
    c(pdf = dir_supp_pdf, svg = dir_supp_svg, png = dir_supp_png)
  } else {
    c(pdf = dir_main_pdf, svg = dir_main_svg, png = dir_main_png)
  }
  paths <- c(
    pdf = file.path(out_dirs[["pdf"]], paste0(basename_no_ext, ".pdf")),
    svg = file.path(out_dirs[["svg"]], paste0(basename_no_ext, ".svg")),
    png = file.path(out_dirs[["png"]], paste0(basename_no_ext, ".png"))
  )
  grDevices::cairo_pdf(paths[["pdf"]], width = width_in, height = height_in,
                       family = FONT_FAMILY, bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::svg(paths[["svg"]], width = width_in, height = height_in,
                 family = FONT_FAMILY, bg = "white", onefile = TRUE)
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(paths[["png"]], width = width_in, height = height_in,
                units = "in", res = dpi, background = "white")
  print(plot)
  grDevices::dev.off()
  tibble(format = names(paths), path = unname(paths))
}

contrast_levels <- c("A_vs_N", "T_vs_A", "T_vs_N")
contrast_labels <- c(A_vs_N = "A vs N", T_vs_A = "T vs A", T_vs_N = "T vs N")
contrast_pairs <- list(
  A_vs_N = c("A", "N"),
  T_vs_A = c("T", "A"),
  T_vs_N = c("T", "N")
)
expected_strict_units <- c(A_vs_N = 0L, T_vs_A = 155L, T_vs_N = 118L)
expected_strict_genes <- c(A_vs_N = 0L, T_vs_A = 99L, T_vs_N = 74L)
expected_dte_units <- c(A_vs_N = 0L, T_vs_A = 4329L, T_vs_N = 4044L)
expected_source_hashes <- c(
  A_vs_N = "7ac71baab936bc589bcd66ed36e1d956528dcd2c944ed7f685cdeb337325a31c",
  T_vs_A = "60074de59d6ed821a4a8be472405179beabe3aa80dc3fa7a9eaa15eff1f0e0f6",
  T_vs_N = "a3d28b4631b04b2fd0b9eef311a1af4a42d3b537c97f82fc0bb06ccb0a91d6c0"
)

input_paths <- c(
  contract = file.path(framework_root, "validated_data_objects.rds"),
  dtu_lock = file.path(v2_base, "tables", "DTU_primary72_strict_source_lock.tsv"),
  dte_counts = file.path(v2_base, "source_data", "Fig2A_DTE_DTU_counts_strict_primary_v5.tsv"),
  dtu_A_vs_N = file.path(phase03_root, "primary72_DTU_isoform_A_vs_N.tsv"),
  dtu_T_vs_A = file.path(phase03_root, "primary72_DTU_isoform_T_vs_A.tsv"),
  dtu_T_vs_N = file.path(phase03_root, "primary72_DTU_isoform_T_vs_N.tsv"),
  old_overlap = file.path(v2_base, "source_data", "Fig3B_detection_overlap_strict_primary_v8.tsv"),
  old_fraction = file.path(v2_base, "source_data", "Fig3C_DTU_only_fraction_strict_primary_v8.tsv"),
  old_gene_effect = file.path(v2_base, "source_data", "Fig3D_gene_log2FC_among_DTU_strict_primary_v8.tsv"),
  old_isoform_effect = file.path(v2_base, "source_data", "Fig3E_isoform_usage_vs_gene_log2FC_strict_primary_v8.tsv")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
stop_if(length(missing_inputs) > 0L, paste0("Missing required input(s): ", paste(missing_inputs, collapse = "; ")))

# Freeze unrelated main figures before writing any new Figure 3 files.
guard_paths <- unlist(lapply(
  c(dir_main_pdf, dir_main_svg, dir_main_png),
  function(d) list.files(d, pattern = "^Fig(1|2|4|5|6|7|8)", full.names = TRUE)
))
guard_before <- file_snapshot(guard_paths)

input_manifest <- tibble(
  input_role = names(input_paths),
  path = unname(input_paths),
  size_bytes = file.info(unname(input_paths))$size,
  mtime = format(file.info(unname(input_paths))$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  sha256 = vapply(unname(input_paths), sha256_file, character(1)),
  access_mode = "read-only"
)

record_progress("INPUT_AUDIT", "Loading validated cross-platform data contract")
contract <- readRDS(input_paths[["contract"]])
required_objects <- c("lr_counts", "sr_gene_counts", "tx2gene", "meta", "sqanti")
stop_if(!all(required_objects %in% names(contract)), "Validated contract is missing one or more required objects.")

lr_counts <- contract[["lr_counts"]]
sr_gene_counts <- contract[["sr_gene_counts"]]
tx2gene <- as.data.frame(contract[["tx2gene"]])
meta <- as.data.frame(contract[["meta"]])
assert_columns(meta, c("sample_id", "patient_id", "tissue"), "meta")
assert_columns(tx2gene, c("isoform_id", "gene_id", "gene_symbol"), "tx2gene")
stop_if(!is.matrix(sr_gene_counts) || !is.numeric(sr_gene_counts), "sr_gene_counts must be a numeric matrix.")
stop_if(anyNA(sr_gene_counts) || any(!is.finite(sr_gene_counts)) || any(sr_gene_counts < 0), "sr_gene_counts contains invalid values.")
stop_if(!identical(colnames(sr_gene_counts), meta$sample_id), "sr_gene_counts sample order differs from meta.")
stop_if(!identical(colnames(lr_counts), meta$sample_id), "lr_counts sample order differs from meta.")
stop_if(anyDuplicated(meta$sample_id) > 0L, "meta contains duplicated sample IDs.")
stop_if(anyDuplicated(rownames(sr_gene_counts)) > 0L, "sr_gene_counts contains duplicated source row IDs.")

primary_keep <- meta$patient_id != "CRC_20"
primary_meta <- meta[primary_keep, , drop = FALSE]
primary_samples <- primary_meta$sample_id
primary_sr_counts_raw <- sr_gene_counts[, primary_samples, drop = FALSE]
stop_if(nrow(primary_meta) != 72L, paste0("primary72 sample count is ", nrow(primary_meta), ", expected 72."))
stop_if(length(unique(primary_meta$patient_id)) != 24L, "primary72 does not contain exactly 24 patients.")
observed_tissue_counts <- table(factor(primary_meta$tissue, levels = c("N", "A", "T")))
stop_if(!identical(as.integer(observed_tissue_counts), c(24L, 24L, 24L)), "primary72 tissue counts are not N/A/T = 24/24/24.")
pairing_table <- table(primary_meta$patient_id, factor(primary_meta$tissue, levels = c("N", "A", "T")))
stop_if(!all(pairing_table == 1L), "One or more primary72 patients lack an exact N/A/T triplet.")
stop_if(!identical(primary_samples, colnames(primary_sr_counts_raw)), "primary72 short-read sample order mismatch.")
stop_if(!identical(primary_samples, colnames(lr_counts)[primary_keep]), "primary72 long-read/short-read sample order mismatch.")

sample_audit <- primary_meta %>%
  mutate(
    primary72_included = TRUE,
    longread_column_present = sample_id %in% colnames(lr_counts),
    shortread_column_present = sample_id %in% colnames(sr_gene_counts),
    exact_cross_platform_match = longread_column_present & shortread_column_present,
    patient_complete_NAT = patient_id %in% rownames(pairing_table)[apply(pairing_table == 1L, 1, all)]
  ) %>%
  arrange(patient_id, factor(tissue, levels = c("N", "A", "T")))

fractional_mask <- abs(primary_sr_counts_raw - round(primary_sr_counts_raw)) > 1e-8
primary_sr_counts <- round(primary_sr_counts_raw)
stop_if(max(primary_sr_counts) > .Machine$integer.max, "Rounded short-read count exceeds R integer range.")
storage.mode(primary_sr_counts) <- "integer"
rounding_audit <- tibble(
  metric = c(
    "input_rows", "input_samples", "fractional_entries", "fractional_entry_fraction",
    "maximum_absolute_rounding_delta", "raw_total_count", "rounded_total_count",
    "maximum_sample_library_relative_delta"
  ),
  value = c(
    nrow(primary_sr_counts_raw), ncol(primary_sr_counts_raw), sum(fractional_mask),
    mean(fractional_mask), max(abs(primary_sr_counts_raw - primary_sr_counts)),
    sum(primary_sr_counts_raw), sum(primary_sr_counts),
    max(abs(colSums(primary_sr_counts_raw) - colSums(primary_sr_counts)) /
          pmax(colSums(primary_sr_counts_raw), 1))
  ),
  interpretation = c(
    "Existing Kallisto/tximport gene rows; no upstream rerun",
    "Exact primary72 sample columns",
    "Estimated-count cells requiring deterministic rounding for DESeq2",
    "Fraction of matrix entries that were not integer-valued",
    "Round-to-nearest was applied once before DESeqDataSetFromMatrix",
    "Sum before rounding",
    "Sum after rounding",
    "Largest relative library-total change caused by rounding"
  )
)

primary_meta <- primary_meta %>%
  mutate(
    patient_id = factor(patient_id),
    tissue = factor(tissue, levels = c("N", "A", "T"))
  )
rownames(primary_meta) <- primary_meta$sample_id

record_progress("MAPPING_AUDIT", "Building stable Ensembl gene-ID mapping without symbol-only matches")
sr_id_map <- tibble(
  sr_gene_id_original = rownames(primary_sr_counts),
  stable_gene_id = strip_ensembl_version(rownames(primary_sr_counts))
) %>%
  group_by(stable_gene_id) %>%
  summarise(
    n_shortread_rows = n(),
    shortread_source_ids = paste(sort(unique(sr_gene_id_original)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    valid_ensembl_gene_id = grepl("^ENSG[0-9]+$", stable_gene_id),
    shortread_mapping_status = case_when(
      !valid_ensembl_gene_id ~ "excluded_non_Ensembl_gene_id",
      n_shortread_rows > 1L ~ "excluded_ambiguous_shortread_rows_after_version_strip",
      TRUE ~ "shortread_unique_stable_Ensembl"
    )
  )

tx_clean <- tx2gene %>%
  transmute(
    isoform_id = clean_chr(isoform_id),
    longread_gene_id_original = clean_chr(gene_id),
    stable_gene_id = strip_ensembl_version(gene_id),
    gene_symbol = clean_chr(gene_symbol)
  )
stop_if(anyDuplicated(tx_clean$isoform_id) > 0L, "tx2gene contains duplicated isoform IDs; ambiguous mapping cannot be resolved automatically.")

lr_gene_map <- tx_clean %>%
  group_by(stable_gene_id) %>%
  summarise(
    n_longread_gene_versions = n_distinct(longread_gene_id_original),
    longread_gene_ids = paste(sort(unique(longread_gene_id_original)), collapse = ";"),
    n_longread_isoforms = n_distinct(isoform_id),
    n_gene_symbols = n_distinct(gene_symbol[!is.na(gene_symbol)]),
    gene_symbol = if (n_gene_symbols == 1L) sort(unique(gene_symbol[!is.na(gene_symbol)]))[[1]] else NA_character_,
    .groups = "drop"
  ) %>%
  mutate(
    valid_ensembl_gene_id = grepl("^ENSG[0-9]+$", stable_gene_id),
    longread_mapping_status = case_when(
      !valid_ensembl_gene_id ~ "excluded_non_Ensembl_gene_id",
      n_gene_symbols > 1L ~ "excluded_ambiguous_gene_symbol_for_stable_id",
      TRUE ~ "longread_unique_stable_Ensembl"
    )
  )

gene_mapping_audit <- full_join(
  sr_id_map %>% select(-valid_ensembl_gene_id),
  lr_gene_map %>% select(-valid_ensembl_gene_id),
  by = "stable_gene_id"
) %>%
  mutate(
    n_shortread_rows = replace_na(n_shortread_rows, 0L),
    n_longread_gene_versions = replace_na(n_longread_gene_versions, 0L),
    n_longread_isoforms = replace_na(n_longread_isoforms, 0L),
    n_gene_symbols = replace_na(n_gene_symbols, 0L),
    valid_ensembl_gene_id = grepl("^ENSG[0-9]+$", stable_gene_id),
    mapping_status = case_when(
      !valid_ensembl_gene_id ~ "excluded_non_Ensembl_gene_id",
      n_shortread_rows > 1L ~ "excluded_ambiguous_shortread_rows_after_version_strip",
      n_gene_symbols > 1L ~ "excluded_ambiguous_longread_annotation",
      n_shortread_rows == 1L & n_longread_isoforms > 0L ~ "shared_unambiguous_Ensembl",
      n_shortread_rows == 1L & n_longread_isoforms == 0L ~ "shortread_only",
      n_shortread_rows == 0L & n_longread_isoforms > 0L ~ "longread_only",
      TRUE ~ "excluded_unresolved"
    ),
    overlap_eligible = mapping_status == "shared_unambiguous_Ensembl",
    symbol_used_for_matching = FALSE
  ) %>%
  arrange(mapping_status, stable_gene_id)

shared_gene_universe <- gene_mapping_audit$stable_gene_id[gene_mapping_audit$overlap_eligible]
stop_if(anyDuplicated(shared_gene_universe) > 0L, "Shared Ensembl gene universe is not unique.")

record_progress("STRICT_DTU_GATE", "Checking satuRn + stageR strict primary call set and immutable counts")
dtu_source_paths <- c(
  A_vs_N = input_paths[["dtu_A_vs_N"]],
  T_vs_A = input_paths[["dtu_T_vs_A"]],
  T_vs_N = input_paths[["dtu_T_vs_N"]]
)
observed_dtu_hashes <- vapply(dtu_source_paths, sha256_file, character(1))
stop_if(!identical(unname(observed_dtu_hashes), unname(expected_source_hashes)), "Strict DTU source SHA256 mismatch; analysis stopped.")

dtu_all <- bind_rows(lapply(contrast_levels, function(contrast_name) {
  x <- readr::read_tsv(dtu_source_paths[[contrast_name]], show_col_types = FALSE)
  assert_columns(
    x,
    c("gene_id", "gene_name", "isoform_id", "contrast", "analysis_cohort", "delta_usage", "padj"),
    paste0("DTU ", contrast_name)
  )
  stop_if(!all(x$contrast == contrast_name), paste0("DTU contrast label mismatch: ", contrast_name))
  stop_if(!all(x$analysis_cohort == "primary72"), paste0("DTU cohort is not primary72: ", contrast_name))
  x %>%
    mutate(
      contrast = contrast_name,
      gene_id = clean_chr(gene_id),
      stable_gene_id = strip_ensembl_version(gene_id),
      source_gene_symbol = clean_chr(gene_name),
      isoform_id = clean_chr(isoform_id),
      delta_usage = as.numeric(delta_usage),
      transcript_adjusted_FDR = as.numeric(padj),
      strict_primary_DTU = is.finite(transcript_adjusted_FDR) &
        transcript_adjusted_FDR < 0.10 & is.finite(delta_usage) & abs(delta_usage) >= 0.05
    )
}))
dtu_strict <- dtu_all %>% filter(strict_primary_DTU)
observed_strict_units <- table(factor(dtu_strict$contrast, levels = contrast_levels))
observed_strict_genes <- vapply(
  contrast_levels,
  function(z) n_distinct(dtu_strict$stable_gene_id[dtu_strict$contrast == z]),
  integer(1)
)
stop_if(!identical(as.integer(observed_strict_units), as.integer(expected_strict_units)), paste0(
  "Strict DTU unit gate failed: observed ", paste(as.integer(observed_strict_units), collapse = "/"),
  "; expected 0/155/118."
))
stop_if(!identical(as.integer(observed_strict_genes), as.integer(expected_strict_genes)), paste0(
  "Strict DTU gene gate failed: observed ", paste(as.integer(observed_strict_genes), collapse = "/"),
  "; expected 0/99/74."
))
stop_if(any(!dtu_strict$stable_gene_id %in% shared_gene_universe), "At least one strict DTU gene is outside the shared unambiguous Ensembl universe.")

dte_source <- readr::read_tsv(input_paths[["dte_counts"]], show_col_types = FALSE)
assert_columns(dte_source, c("analysis_type", "statistical_unit", "contrast", "count", "threshold"), "DTE count source")
dte_counts <- dte_source %>%
  filter(analysis_type == "DTE events", contrast %in% contrast_levels) %>%
  mutate(contrast = as.character(contrast), count = as.integer(count)) %>%
  arrange(match(contrast, contrast_levels))
stop_if(nrow(dte_counts) != 3L, "Frozen long-read DTE source does not contain exactly three contrasts.")
stop_if(!identical(dte_counts$count, as.integer(expected_dte_units)), "Frozen long-read DTE counts differ from 0/4329/4044.")

input_manifest <- input_manifest %>%
  mutate(
    source_lock_status = case_when(
      input_role == "dtu_A_vs_N" ~ if_else(sha256 == expected_source_hashes[["A_vs_N"]], "PASS", "FAIL"),
      input_role == "dtu_T_vs_A" ~ if_else(sha256 == expected_source_hashes[["T_vs_A"]], "PASS", "FAIL"),
      input_role == "dtu_T_vs_N" ~ if_else(sha256 == expected_source_hashes[["T_vs_N"]], "PASS", "FAIL"),
      TRUE ~ "not_applicable"
    )
  )

record_progress("DESEQ2", "Running one paired DESeq2 model on rounded existing tximport estimated counts")
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = primary_sr_counts,
  colData = primary_meta,
  design = ~ patient_id + tissue
)
dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = FALSE)

extract_deseq_result <- function(contrast_name) {
  pair <- contrast_pairs[[contrast_name]]
  res <- DESeq2::results(
    dds,
    contrast = c("tissue", pair[[1]], pair[[2]]),
    alpha = 0.05,
    independentFiltering = TRUE,
    cooksCutoff = TRUE
  )
  out <- as.data.frame(res) %>%
    tibble::rownames_to_column("shortread_gene_id") %>%
    mutate(
      stable_gene_id = strip_ensembl_version(shortread_gene_id),
      contrast = contrast_name,
      numerator = pair[[1]],
      denominator = pair[[2]],
      shortread_DEG = !is.na(padj) & padj < 0.05 & !is.na(log2FoldChange) & abs(log2FoldChange) > 1,
      DE_threshold = "DESeq2 adjusted P < 0.05 and abs(log2FC) > 1",
      model = "DESeq2 paired fixed-effect design: ~ patient_id + tissue",
      count_input = "existing Kallisto/tximport estimated counts rounded once to nearest integer"
    ) %>%
    left_join(
      gene_mapping_audit %>%
        select(stable_gene_id, gene_symbol, mapping_status, overlap_eligible),
      by = "stable_gene_id"
    ) %>%
    mutate(
      mapping_status = replace_na(mapping_status, "excluded_unresolved"),
      overlap_eligible = replace_na(overlap_eligible, FALSE)
    )
  out
}

shortread_de_results <- bind_rows(lapply(contrast_levels, extract_deseq_result)) %>%
  select(
    contrast, numerator, denominator, shortread_gene_id, stable_gene_id, gene_symbol,
    baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, shortread_DEG,
    mapping_status, overlap_eligible, DE_threshold, model, count_input
  ) %>%
  arrange(factor(contrast, levels = contrast_levels), stable_gene_id)

stop_if(nrow(shortread_de_results) != nrow(sr_gene_counts) * 3L, "DESeq2 result row count mismatch.")
stop_if(anyDuplicated(shortread_de_results[c("contrast", "stable_gene_id")]) > 0L, "DESeq2 result contains duplicate contrast/stable-gene rows.")

shortread_de_shared <- shortread_de_results %>% filter(overlap_eligible)

record_progress("INTEGRATION", "Recomputing strict DTU versus matched short-read DE overlap")
overlap_results <- bind_rows(lapply(c("T_vs_A", "T_vs_N"), function(contrast_name) {
  dtu_genes <- sort(unique(dtu_strict$stable_gene_id[dtu_strict$contrast == contrast_name]))
  sr_de_total <- sort(unique(shortread_de_results$stable_gene_id[
    shortread_de_results$contrast == contrast_name & shortread_de_results$shortread_DEG
  ]))
  sr_de_shared <- sort(unique(shortread_de_shared$stable_gene_id[
    shortread_de_shared$contrast == contrast_name & shortread_de_shared$shortread_DEG
  ]))
  both <- intersect(dtu_genes, sr_de_shared)
  dtu_only <- setdiff(dtu_genes, sr_de_shared)
  sr_only <- setdiff(sr_de_shared, dtu_genes)
  tibble(
    contrast = contrast_name,
    n_shared_unambiguous_Ensembl_genes = length(shared_gene_universe),
    n_strict_longread_DTU_units = sum(dtu_strict$contrast == contrast_name),
    n_strict_longread_DTU_genes = length(dtu_genes),
    n_shortread_DE_genes_all_Ensembl = length(sr_de_total),
    n_shortread_DE_genes_shared_unambiguous = length(sr_de_shared),
    n_DTU_only_genes = length(dtu_only),
    n_DTU_plus_shortread_DE_genes = length(both),
    n_shortread_DE_only_genes = length(sr_only),
    fraction_DTU_not_shortread_DE = length(dtu_only) / length(dtu_genes),
    percent_DTU_not_shortread_DE = 100 * fraction_DTU_not_shortread_DE,
    DTU_threshold = "satuRn + stageR transcript adjusted FDR < 0.10 and abs(delta usage) >= 0.05",
    shortread_DE_threshold = "DESeq2 adjusted P < 0.05 and abs(log2FC) > 1",
    overlap_universe = "shared unambiguous version-stripped Ensembl gene IDs; no gene-symbol matching"
  )
}))

dtu_shortread_gene_effect <- dtu_strict %>%
  select(
    contrast, gene_id, stable_gene_id, source_gene_symbol, isoform_id,
    delta_usage, transcript_adjusted_FDR, analysis_cohort
  ) %>%
  left_join(
    lr_gene_map %>% select(stable_gene_id, tx2gene_symbol = gene_symbol),
    by = "stable_gene_id"
  ) %>%
  mutate(gene_symbol = coalesce(tx2gene_symbol, source_gene_symbol)) %>%
  select(-tx2gene_symbol) %>%
  left_join(
    shortread_de_results %>%
      select(
        contrast, stable_gene_id, shortread_baseMean = baseMean,
        shortread_log2FC = log2FoldChange, shortread_lfcSE = lfcSE,
        shortread_pvalue = pvalue, shortread_padj = padj,
        shortread_DEG, mapping_status, overlap_eligible
      ),
    by = c("contrast", "stable_gene_id")
  ) %>%
  mutate(
    mapping_level = if_else(overlap_eligible, "shared_unambiguous_Ensembl", mapping_status),
    direction_relation = case_when(
      !is.finite(shortread_log2FC) ~ "shortread_effect_unavailable",
      sign(delta_usage) == sign(shortread_log2FC) ~ "same_sign",
      sign(delta_usage) != sign(shortread_log2FC) ~ "opposite_sign",
      TRUE ~ "zero_effect"
    ),
    DTU_threshold = "transcript adjusted FDR < 0.10 and abs(delta usage) >= 0.05",
    shortread_DE_threshold = "adjusted P < 0.05 and abs(log2FC) > 1"
  ) %>%
  select(
    contrast, gene_id, stable_gene_id, gene_symbol, isoform_id,
    delta_usage, transcript_adjusted_FDR, analysis_cohort,
    shortread_baseMean, shortread_log2FC, shortread_lfcSE,
    shortread_pvalue, shortread_padj, shortread_DEG,
    mapping_level, direction_relation, DTU_threshold, shortread_DE_threshold
  ) %>%
  arrange(factor(contrast, levels = contrast_levels), stable_gene_id, isoform_id)

stop_if(nrow(dtu_shortread_gene_effect) != sum(expected_strict_units), "Strict DTU-to-short-read effect table row count is not 273.")
stop_if(any(is.na(dtu_shortread_gene_effect$shortread_log2FC)), "At least one strict DTU row lacks a matched short-read gene log2FC.")

panel_a_source <- bind_rows(
  shortread_de_results %>%
    group_by(contrast) %>%
    summarise(count = n_distinct(stable_gene_id[shortread_DEG]), .groups = "drop") %>%
    mutate(
      analysis_type = "Short-read DEG genes",
      statistical_unit = "genes; all Ensembl rows in the matched short-read analysis",
      threshold = "DESeq2 adjusted P < 0.05 and abs(log2FC) > 1"
    ),
  dte_counts %>%
    transmute(
      contrast, count,
      analysis_type = "Long-read DTE",
      statistical_unit = statistical_unit,
      threshold = threshold
    ),
  tibble(
    contrast = contrast_levels,
    count = as.integer(expected_strict_units),
    analysis_type = "Long-read strict DTU",
    statistical_unit = "isoform-level DTU calls",
    threshold = "satuRn + stageR transcript adjusted FDR < 0.10 and abs(delta usage) >= 0.05"
  )
) %>%
  mutate(
    contrast = factor(contrast, levels = contrast_levels),
    comparison = factor(unname(contrast_labels[as.character(contrast)]), levels = unname(contrast_labels)),
    analysis_type = factor(
      analysis_type,
      levels = c("Short-read DEG genes", "Long-read DTE", "Long-read strict DTU")
    )
  )

panel_b_source <- overlap_results %>%
  select(
    contrast, n_DTU_only_genes, n_DTU_plus_shortread_DE_genes,
    n_shortread_DE_only_genes, n_shared_unambiguous_Ensembl_genes
  ) %>%
  pivot_longer(
    cols = c(n_DTU_only_genes, n_DTU_plus_shortread_DE_genes, n_shortread_DE_only_genes),
    names_to = "set_name", values_to = "n_genes"
  ) %>%
  mutate(
    set_label = recode(
      set_name,
      n_DTU_only_genes = "DTU only",
      n_DTU_plus_shortread_DE_genes = "DTU + short-read DE",
      n_shortread_DE_only_genes = "Short-read DE only"
    ),
    set_label = factor(set_label, levels = c("DTU only", "DTU + short-read DE", "Short-read DE only")),
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N"))
  )

panel_c_source <- overlap_results %>%
  transmute(
    contrast,
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    n_DTU_genes = n_strict_longread_DTU_genes,
    n_DTU_not_shortread_DE = n_DTU_only_genes,
    fraction = fraction_DTU_not_shortread_DE,
    percent = percent_DTU_not_shortread_DE,
    label = sprintf("%.1f%% (%d/%d)", percent, n_DTU_not_shortread_DE, n_DTU_genes)
  )

panel_c_plot <- panel_c_source %>%
  mutate(label_display = sprintf("%.1f%%\n%d/%d", percent, n_DTU_not_shortread_DE, n_DTU_genes))

panel_d_source <- dtu_shortread_gene_effect %>%
  filter(contrast %in% c("T_vs_A", "T_vs_N")) %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    shortread_status = if_else(shortread_DEG, "Short-read DE", "Not short-read DE"),
    shortread_status = factor(shortread_status, levels = c("Not short-read DE", "Short-read DE"))
  )

record_progress("OUTPUT_TABLES", "Writing complete DE results, overlap, mapping and panel source data")
output_paths <- c(
  matched_de = write_tsv2(shortread_de_results, file.path(dir_tables, "matched_shortread_primary72_DE_results.tsv")),
  overlap = write_tsv2(overlap_results, file.path(dir_tables, "DTU_shortread_DE_overlap.tsv")),
  gene_effect = write_tsv2(dtu_shortread_gene_effect, file.path(dir_tables, "DTU_shortread_gene_effect.tsv")),
  sample_audit = write_tsv2(sample_audit, file.path(dir_tables, "matched_shortread_primary72_sample_audit.tsv")),
  mapping_audit = write_tsv2(gene_mapping_audit, file.path(dir_tables, "matched_shortread_gene_mapping_audit.tsv")),
  rounding_audit = write_tsv2(rounding_audit, file.path(dir_tables, "matched_shortread_count_rounding_audit.tsv")),
  input_manifest = write_tsv2(input_manifest, file.path(dir_tables, "matched_shortread_input_manifest.tsv")),
  panel_A = write_tsv2(panel_a_source, file.path(dir_source, "Figure3A_matched_shortread_longread_counts.tsv")),
  panel_B = write_tsv2(panel_b_source, file.path(dir_source, "Figure3B_DTU_shortread_overlap.tsv")),
  panel_C = write_tsv2(panel_c_source, file.path(dir_source, "Figure3C_DTU_not_shortread_DE_fraction.tsv")),
  panel_D = write_tsv2(panel_d_source, file.path(dir_source, "Figure3D_delta_usage_vs_shortread_gene_log2FC.tsv"))
)

COL <- c(
  "Short-read DEG genes" = "#2F6B8A",
  "Long-read DTE" = "#D97928",
  "Long-read strict DTU" = "#7B4E9D",
  "DTU only" = "#D97928",
  "DTU + short-read DE" = "#7B4E9D",
  "Short-read DE only" = "#7E8B91",
  "Not short-read DE" = "#B9BEC1",
  "Short-read DE" = "#B53A3A"
)

panel_a <- ggplot(panel_a_source, aes(comparison, count, fill = analysis_type)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(y = pmax(count, ave(count, analysis_type, FUN = max) * 0.045), label = comma(count)),
    vjust = -0.35, size = 2.55
  ) +
  facet_wrap(~ analysis_type, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = COL, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "A  Cross-analysis detection counts",
    subtitle = "Short-read counts are genes; long-read DTE/DTU counts are isoform-level calls",
    x = NULL, y = "Count"
  ) +
  theme_manuscript(7.1) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    strip.text = element_text(size = 6.8, face = "bold"),
    panel.spacing.x = grid::unit(6, "mm")
  )

panel_b <- ggplot(panel_b_source, aes(set_label, comparison)) +
  geom_tile(aes(fill = set_label), width = 0.88, height = 0.68, colour = "white", linewidth = 0.7) +
  geom_text(aes(label = comma(n_genes)), size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = COL, guide = "none") +
  labs(
    title = "B  DTU gene and short-read DE overlap",
    subtitle = "Exact disjoint-set counts; shared unambiguous Ensembl-gene universe",
    x = NULL, y = NULL
  ) +
  theme_manuscript(7.2) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.line = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

panel_c <- ggplot(panel_c_plot, aes(comparison, percent, fill = comparison)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label_display), vjust = -0.35, size = 2.7, lineheight = 0.9) +
  scale_fill_manual(values = c("T vs A" = "#D97928", "T vs N" = "#7B4E9D"), guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
  labs(
    title = "C  DTU genes below short-read DE threshold",
    subtitle = "Numerator/denominator shown above each bar",
    x = NULL, y = "DTU genes not called as short-read DE"
  ) +
  theme_manuscript(7.2)

panel_d <- ggplot(panel_d_source, aes(shortread_log2FC, delta_usage)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, colour = "#7A7A7A") +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", linewidth = 0.3, colour = "#7A7A7A") +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#B0B0B0") +
  geom_point(aes(fill = shortread_status), shape = 21, colour = "white", stroke = 0.18,
             size = 1.8, alpha = 0.82) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = COL, drop = FALSE) +
  labs(
    title = "D  Isoform usage change versus matched gene expression",
    subtitle = "Each point is one strict long-read DTU isoform call",
    x = "Matched short-read gene log2FC",
    y = "Long-read isoform delta usage",
    fill = NULL
  ) +
  theme_manuscript(7.2) +
  theme(legend.position = "top", panel.grid.major.x = element_line(linewidth = 0.16, colour = "#E8E8E8"))

main_figure <- panel_a / (panel_b | panel_c) / panel_d +
  plot_layout(heights = c(0.86, 0.92, 1.20)) +
  plot_annotation(
    title = "Matched short-read gene expression contextualizes long-read isoform remodeling",
    subtitle = "Paired primary72 cohort; matched complementary analysis of aliquots from the same RNA extracts",
    theme = theme(
      plot.title = element_text(family = FONT_FAMILY, face = "bold", size = 11.2, colour = "#202020"),
      plot.subtitle = element_text(family = FONT_FAMILY, size = 8.0, colour = "#555555"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

record_progress("FIGURE3", "Rendering matched short-read/long-read Figure 3")
main_figure_paths <- save_bundle(
  main_figure, "Figure3_matched_shortread_longread",
  width_mm = 183, height_mm = 174, supplementary = FALSE, dpi = 600
)

record_progress("SUPPLEMENT", "Moving the internal long-read-derived gene screen to a supplementary rendering")
old_overlap <- readr::read_tsv(input_paths[["old_overlap"]], show_col_types = FALSE)
old_fraction <- readr::read_tsv(input_paths[["old_fraction"]], show_col_types = FALSE)
old_gene_effect <- readr::read_tsv(input_paths[["old_gene_effect"]], show_col_types = FALSE)
old_isoform_effect <- readr::read_tsv(input_paths[["old_isoform_effect"]], show_col_types = FALSE)
assert_columns(old_overlap, c("contrast", "detection_class", "n_genes"), "old internal overlap")
assert_columns(old_fraction, c("contrast", "n_DTU_only", "n_total_DTU_genes", "DTU_only_fraction", "estimable"), "old internal fraction")
assert_columns(old_gene_effect, c("contrast", "gene_log2FC", "detection_class"), "old internal gene effect")
assert_columns(old_isoform_effect, c("contrast", "delta_usage", "gene_log2FC", "detection_class"), "old internal isoform effect")

supp_a_data <- old_overlap %>%
  filter(contrast %in% c("T_vs_A", "T_vs_N")) %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    detection_class = factor(
      detection_class,
      levels = c("DTU only", "DTU + derived gene-level", "Derived gene-level only")
    )
  )
supp_a <- ggplot(supp_a_data, aes(detection_class, log10(n_genes + 1), fill = detection_class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = comma(n_genes)), vjust = -0.3, size = 2.5) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = c(
    "DTU only" = "#D97928", "DTU + derived gene-level" = "#7B4E9D",
    "Derived gene-level only" = "#9A9A9A"
  ), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "A  Internal derived-gene overlap", x = NULL, y = "log10(gene count + 1)") +
  theme_manuscript(7.1) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

supp_b_data <- old_fraction %>%
  filter(contrast %in% c("T_vs_A", "T_vs_N"), estimable) %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    percent = 100 * DTU_only_fraction,
    label = sprintf("%.1f%%\n%d/%d", percent, n_DTU_only, n_total_DTU_genes)
  )
supp_b <- ggplot(supp_b_data, aes(comparison, percent, fill = comparison)) +
  geom_col(width = 0.58) +
  geom_text(aes(label = label), vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = c("T vs A" = "#D97928", "T vs N" = "#7B4E9D"), guide = "none") +
  scale_y_continuous(limits = c(0, 108), labels = function(x) paste0(x, "%")) +
  labs(title = "B  DTU-only fraction in internal screen", x = NULL, y = "DTU genes without derived gene-level call") +
  theme_manuscript(7.1)

supp_c <- old_gene_effect %>%
  filter(contrast %in% c("T_vs_A", "T_vs_N")) %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    detection_class = factor(detection_class, levels = c("DTU only", "DTU + derived gene-level"))
  ) %>%
  ggplot(aes(detection_class, gene_log2FC, fill = detection_class)) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.45) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = c("DTU only" = "#D97928", "DTU + derived gene-level" = "#7B4E9D"), guide = "none") +
  labs(title = "C  Internal derived gene log2FC", x = NULL, y = "Long-read-derived gene log2FC") +
  theme_manuscript(7.1) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

supp_d <- old_isoform_effect %>%
  filter(contrast %in% c("T_vs_A", "T_vs_N")) %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    detection_class = factor(detection_class, levels = c("DTU only", "DTU + derived gene-level"))
  ) %>%
  ggplot(aes(gene_log2FC, delta_usage, fill = detection_class)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", linewidth = 0.3, colour = "#777777") +
  geom_point(shape = 21, colour = "white", stroke = 0.15, size = 1.6, alpha = 0.8) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = c("DTU only" = "#D97928", "DTU + derived gene-level" = "#7B4E9D")) +
  labs(
    title = "D  Usage versus internally derived gene effect",
    x = "Long-read-derived gene log2FC", y = "Long-read delta usage", fill = NULL
  ) +
  theme_manuscript(7.1) +
  theme(legend.position = "top")

supp_figure <- (supp_a | supp_b) / (supp_c | supp_d) +
  plot_annotation(
    title = "Supplementary: internal long-read-derived gene-level screen",
    subtitle = "Aggregated long-read counts; paired log2CPM Wilcoxon + BH. This is not an independent gene-expression atlas.",
    theme = theme(
      plot.title = element_text(family = FONT_FAMILY, face = "bold", size = 10.8),
      plot.subtitle = element_text(family = FONT_FAMILY, size = 7.8, colour = "#555555")
    )
  )
supp_figure_paths <- save_bundle(
  supp_figure, "Supplementary_internal_longread_gene_screen",
  width_mm = 183, height_mm = 145, supplementary = TRUE, dpi = 600
)

shortread_summary <- shortread_de_results %>%
  group_by(contrast) %>%
  summarise(
    n_rows = n(),
    n_tested_with_finite_p = sum(is.finite(pvalue)),
    n_padj_available = sum(!is.na(padj)),
    n_DEG_all_Ensembl = sum(shortread_DEG),
    n_DEG_shared_unambiguous = sum(shortread_DEG & overlap_eligible),
    n_up_all_Ensembl = sum(shortread_DEG & log2FoldChange > 1),
    n_down_all_Ensembl = sum(shortread_DEG & log2FoldChange < -1),
    .groups = "drop"
  )

key_statistics <- bind_rows(
  tibble(section = "cohort", metric = c("patients", "samples", "N", "A", "T"), value = c(24, 72, 24, 24, 24), contrast = NA_character_),
  tibble(section = "strict_DTU_units", metric = "n", value = as.numeric(expected_strict_units), contrast = contrast_levels),
  tibble(section = "strict_DTU_genes", metric = "n", value = as.numeric(expected_strict_genes), contrast = contrast_levels),
  shortread_summary %>% transmute(section = "shortread_DEG_all", metric = "n", value = n_DEG_all_Ensembl, contrast),
  shortread_summary %>% transmute(section = "shortread_DEG_shared", metric = "n", value = n_DEG_shared_unambiguous, contrast),
  overlap_results %>% transmute(section = "DTU_only", metric = "n_genes", value = n_DTU_only_genes, contrast),
  overlap_results %>% transmute(section = "DTU_plus_shortread_DE", metric = "n_genes", value = n_DTU_plus_shortread_DE_genes, contrast),
  overlap_results %>% transmute(section = "DTU_not_shortread_DE", metric = "percent", value = percent_DTU_not_shortread_DE, contrast)
)
output_paths <- c(
  output_paths,
  key_statistics = write_tsv2(key_statistics, file.path(dir_tables, "matched_shortread_key_statistics.tsv"))
)

source_index <- bind_rows(
  tibble(
    figure = "Figure 3",
    panel = c("A", "B", "C", "D"),
    source_data_path = unname(output_paths[c("panel_A", "panel_B", "panel_C", "panel_D")]),
    n_rows = c(nrow(panel_a_source), nrow(panel_b_source), nrow(panel_c_source), nrow(panel_d_source)),
    description = c(
      "All matched short-read DEG genes, frozen long-read DTE isoforms and strict long-read DTU isoforms",
      "Disjoint DTU/short-read-DE gene sets on the shared Ensembl universe",
      "Fraction of strict DTU genes not meeting the matched short-read DEG threshold",
      "Strict DTU isoform delta usage joined to matched short-read gene log2FC"
    )
  ),
  tibble(
    figure = "Supplementary internal long-read gene screen",
    panel = c("A", "B", "C", "D"),
    source_data_path = unname(input_paths[c("old_overlap", "old_fraction", "old_gene_effect", "old_isoform_effect")]),
    n_rows = c(nrow(old_overlap), nrow(old_fraction), nrow(old_gene_effect), nrow(old_isoform_effect)),
    description = c(
      "Existing strict DTU versus long-read-derived gene screen overlap",
      "Existing internal DTU-only fraction",
      "Existing internal gene log2FC among DTU genes",
      "Existing strict DTU delta usage versus internally derived gene log2FC"
    )
  )
)
output_paths <- c(
  output_paths,
  source_index = write_tsv2(source_index, file.path(dir_tables, "Figure3_matched_shortread_source_data_index.tsv"))
)

ov_ta <- overlap_results %>% filter(contrast == "T_vs_A")
ov_tn <- overlap_results %>% filter(contrast == "T_vs_N")
sr_summary_text <- shortread_summary %>%
  mutate(text = sprintf(
    "- %s：全 Ensembl 检出 %s 个 DEG；共享无歧义 Ensembl universe 中 %s 个。",
    unname(contrast_labels[contrast]), comma(n_DEG_all_Ensembl), comma(n_DEG_shared_unambiguous)
  )) %>%
  pull(text)

report_lines <- c(
  "# Matched short-read gene-level reference analysis 报告",
  "",
  "## 分析定位",
  "",
  "本分析使用与 Oxford Nanopore 长读长测序相同患者、同一批 RNA 提取物分装所得的 BGI 短读长数据，作为 matched complementary gene-level analysis。它不是独立外部验证，也不能验证完整全长 isoform 结构。",
  "",
  "## 冻结输入与停止门槛",
  "",
  "- 主分析集固定为 primary72：24 位患者、72 份组织，N/A/T 各 24 份；CRC-20 的 3 份样本已排除。",
  "- 24 位患者均为完整 N/A/T 三联配对；长短读长 sample_id 集合和顺序完全一致。",
  "- 长读长 strict DTU 唯一来源为已审计的 satuRn + stageR primary call set。严格阈值为 transcript-level adjusted FDR (`padj`) < 0.10 且 |delta usage| >= 0.05。",
  "- strict DTU unit 数量核验为 A vs N = 0、T vs A = 155、T vs N = 118；门槛通过。",
  "- 本阶段没有重跑 Kallisto、FASTQ、STAR、FLAIR、SQANTI3 或 DTU/DTE 上游分析。",
  "",
  "## 短读长 DESeq2 方法",
  "",
  "- 输入为数据契约中已存在的 Kallisto/tximport gene-level estimated counts。",
  "- 因 DESeq2 的 matrix 接口要求整数计数，对 primary72 estimated counts 统一执行一次 round-to-nearest；未改变样本、基因或阈值。取整影响详见 `tables/matched_shortread_count_rounding_audit.tsv`。",
  "- 配对模型：`design = ~ patient_id + tissue`；tissue levels 固定为 N、A、T。",
  "- contrasts：A vs N、T vs A、T vs N。",
  "- 短读长 DEG 定义固定为 adjusted P < 0.05 且 |log2FC| > 1。未按结果调整阈值。",
  "- DESeq2 `results()` 使用独立过滤生成 adjusted P；未增加结果导向的手工低表达筛选。",
  "",
  "## Ensembl 映射",
  "",
  paste0("- short-read 共 ", comma(nrow(sr_id_map)), " 个 stable Ensembl gene IDs；去版本后重复行数为 ", comma(sum(sr_id_map$n_shortread_rows > 1L)), "。"),
  paste0("- long-read tx2gene 中共有 ", comma(nrow(lr_gene_map)), " 个 stable gene IDs；与短读长共享且无歧义的 Ensembl gene universe 为 ", comma(length(shared_gene_universe)), " 个。"),
  "- overlap 只使用去版本后的 Ensembl gene ID；gene symbol 仅用于显示，从未用于匹配。",
  "- 一对多、去版本后重复、非 Ensembl 或注释歧义项均从 overlap universe 排除并记录于 `tables/matched_shortread_gene_mapping_audit.tsv`。",
  "- 99 个 T vs A strict DTU genes 和 74 个 T vs N strict DTU genes 均成功进入共享无歧义 Ensembl universe。",
  "",
  "## 短读长 DEG 结果",
  "",
  sr_summary_text,
  "",
  "## strict DTU 与 matched short-read DEG 的关系",
  "",
  sprintf(
    "- T vs A：strict DTU genes = %d；DTU only = %d；DTU + short-read DE = %d；short-read DE only = %d；%.1f%% (%d/%d) 的 DTU genes 未达到短读长 DEG 标准。",
    ov_ta$n_strict_longread_DTU_genes, ov_ta$n_DTU_only_genes,
    ov_ta$n_DTU_plus_shortread_DE_genes, ov_ta$n_shortread_DE_only_genes,
    ov_ta$percent_DTU_not_shortread_DE, ov_ta$n_DTU_only_genes, ov_ta$n_strict_longread_DTU_genes
  ),
  sprintf(
    "- T vs N：strict DTU genes = %d；DTU only = %d；DTU + short-read DE = %d；short-read DE only = %d；%.1f%% (%d/%d) 的 DTU genes 未达到短读长 DEG 标准。",
    ov_tn$n_strict_longread_DTU_genes, ov_tn$n_DTU_only_genes,
    ov_tn$n_DTU_plus_shortread_DE_genes, ov_tn$n_shortread_DE_only_genes,
    ov_tn$percent_DTU_not_shortread_DE, ov_tn$n_DTU_only_genes, ov_tn$n_strict_longread_DTU_genes
  ),
  "- 上述 overlap 的统计单位是 gene；Figure 3D 的统计单位是 strict isoform-level DTU call。两者不得混写。",
  "",
  "## Figure 3 结构",
  "",
  "- Panel A：三个 contrast 中 matched short-read DEG genes、long-read DTE isoforms 和 strict DTU isoform calls 的数量；facet 明确标注不同统计单位。",
  "- Panel B：T vs A 和 T vs N 的 DTU-only、DTU + short-read DE、short-read-DE-only 三类互斥 gene set。",
  "- Panel C：strict DTU genes 未达到 matched short-read DEG 标准的比例。",
  "- Panel D：每个 strict DTU isoform 的 delta usage 与其所属 gene 的 matched short-read log2FC；颜色仅表示该 gene 是否达到固定 DEG 阈值。",
  "- 原 long-read-derived gene-level screen 已保留并重新渲染为 Supplementary Figure；原始结果未删除。",
  "",
  "## 可写入 manuscript 的结论",
  "",
  "- 在完全匹配的同患者、同 RNA 提取物分装数据中，短读长 gene-level DE 与长读长 strict DTU 只部分重叠。",
  "- 一部分 strict DTU genes 未达到固定的 matched short-read gene-level DEG 阈值，支持 gene expression 与 isoform usage 是互补分析层面。",
  "- 该比较控制了患者来源和 RNA 输入差异，但仍受平台、文库构建、定量目标及统计模型差异影响。",
  "",
  "## 不能写的结论",
  "",
  "- 不能写成独立外部验证或独立队列复现。",
  "- 不能写成短读长验证了完整全长 isoform。",
  "- 不能把未达到 DEG 阈值解释为 gene expression 完全不变。",
  "- 不能从事件数量差异直接推断长读长平台优于短读长平台；Panel A 的统计单位不同。",
  "- 不能引入 N-A-T score、ROC、MSI、ESTIMATE、单细胞轨迹、临床关联或候选评分。",
  "",
  "## 文件路径",
  "",
  paste0("- Main PDF: `", main_figure_paths$path[main_figure_paths$format == "pdf"], "`"),
  paste0("- Main SVG: `", main_figure_paths$path[main_figure_paths$format == "svg"], "`"),
  paste0("- Main PNG: `", main_figure_paths$path[main_figure_paths$format == "png"], "`"),
  paste0("- Supplementary PDF: `", supp_figure_paths$path[supp_figure_paths$format == "pdf"], "`"),
  "- Complete DE results: `tables/matched_shortread_primary72_DE_results.tsv`",
  "- Overlap summary: `tables/DTU_shortread_DE_overlap.tsv`",
  "- DTU-to-short-read gene effects: `tables/DTU_shortread_gene_effect.tsv`"
)
report_path <- write_lines2(report_lines, file.path(dir_reports, "MATCHED_SHORTREAD_INTEGRATION_REPORT.md"))

all_outputs <- c(
  unname(output_paths), unname(main_figure_paths$path), unname(supp_figure_paths$path), report_path
)
stop_if(any(!file.exists(all_outputs)), "One or more required outputs were not created.")
stop_if(any(file.info(all_outputs)$size <= 0), "One or more required outputs are empty.")

guard_after <- file_snapshot(guard_paths)
guard_check <- full_join(guard_before, guard_after, by = "path", suffix = c("_before", "_after")) %>%
  mutate(unchanged = size_bytes_before == size_bytes_after & sha256_before == sha256_after)
stop_if(any(!guard_check$unchanged), "An unrelated Figure 1/2/4-8 main-figure artifact changed.")
write_tsv2(guard_check, file.path(dir_tables, "matched_shortread_unrelated_figure_guard.tsv"))

record_progress(
  "COMPLETE",
  paste0(
    "DESeq2 complete; strict DTU gate 0/155/118; T_vs_A missed=",
    sprintf("%.1f%%", ov_ta$percent_DTU_not_shortread_DE),
    "; T_vs_N missed=", sprintf("%.1f%%", ov_tn$percent_DTU_not_shortread_DE)
  )
)

cat("Matched short-read integration completed.\n")
cat("This is a matched complementary gene-level analysis of aliquots from the same RNA extracts, not independent external validation.\n")
cat("Report: ", report_path, "\n", sep = "")
