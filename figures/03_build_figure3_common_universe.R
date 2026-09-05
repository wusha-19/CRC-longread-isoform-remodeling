#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(20260828)

required_packages <- c(
  "ggplot2", "patchwork", "dplyr", "tidyr", "readr", "scales",
  "ragg", "systemfonts"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

v2_base <- Sys.getenv("CRC_V2_ROOT")
if (!nzchar(v2_base)) stop("Set CRC_V2_ROOT before running this script.")
study_root <- dirname(v2_base)
framework_root <- file.path(study_root, "16_paired_NAT_multi_platform_framework_20260827")
phase03_root <- file.path(study_root, "03_DTU_paired", "tables")

dir_results <- file.path(v2_base, "results")
dir_source <- file.path(v2_base, "source_data", "Figure3_v2")
dir_reports <- file.path(v2_base, "reports")
dir_logs <- file.path(v2_base, "logs")
dir_main_pdf <- file.path(v2_base, "figures", "main", "pdf")
dir_main_svg <- file.path(v2_base, "figures", "main", "svg")
dir_main_png <- file.path(v2_base, "figures", "main", "png")
dir_supp_pdf <- file.path(v2_base, "figures", "supplementary", "pdf")
dir_supp_svg <- file.path(v2_base, "figures", "supplementary", "svg")
dir_supp_png <- file.path(v2_base, "figures", "supplementary", "png")
dir_archive <- file.path(v2_base, "archive")
for (d in c(
  dir_results, dir_source, dir_reports, dir_logs,
  dir_main_pdf, dir_main_svg, dir_main_png,
  dir_supp_pdf, dir_supp_svg, dir_supp_png, dir_archive
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

progress_path <- file.path(dir_logs, "28_Figure3_common_universe_v2_progress.txt")
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

first_valid <- function(x) {
  y <- clean_chr(x)
  y <- y[!is.na(y)]
  if (length(y) == 0L) NA_character_ else y[[1]]
}

sha256_file <- function(path) {
  stop_if(!file.exists(path), paste0("Missing file for SHA256: ", path))
  out <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
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

archive_existing <- function(paths) {
  existing <- paths[file.exists(paths) & !dir.exists(paths)]
  if (length(existing) == 0L) return(NA_character_)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  archive_dir <- file.path(dir_archive, paste0("Figure3_common_universe_v2_before_", stamp))
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(existing, file.path(archive_dir, basename(existing)), overwrite = FALSE)
  stop_if(!all(ok), "Could not archive all pre-existing Figure 3 v2 outputs.")
  archive_dir
}

count_tsv_rows <- function(path) {
  max(0L, length(readLines(path, warn = FALSE)) - 1L)
}

tsv_field_audit <- function(paths) {
  bind_rows(lapply(paths, function(path) {
    lines <- readLines(path, warn = FALSE)
    counts <- lengths(regmatches(lines, gregexpr("\t", lines, fixed = TRUE))) + 1L
    tibble(
      path = path,
      n_data_rows = max(0L, length(lines) - 1L),
      header_fields = counts[[1]],
      min_fields = min(counts),
      max_fields = max(counts),
      field_count_pass = min(counts) == max(counts)
    )
  }))
}

FONT_FAMILY <- "Times New Roman"
font_match <- systemfonts::match_fonts(FONT_FAMILY)[1, , drop = FALSE]
stop_if(
  is.null(font_match$path) || !nzchar(font_match$path) || !file.exists(font_match$path),
  "Times New Roman is unavailable in the remote R graphics environment."
)

theme_manuscript <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY, colour = "#202020"),
      axis.line = element_line(linewidth = 0.28, colour = "#333333"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#333333"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.7),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 1.0, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#555555", margin = margin(b = 3)),
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
  grDevices::cairo_pdf(
    paths[["pdf"]], width = width_in, height = height_in,
    family = FONT_FAMILY, bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  grDevices::svg(
    paths[["svg"]], width = width_in, height = height_in,
    family = FONT_FAMILY, bg = "white", onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    paths[["png"]], width = width_in, height = height_in,
    units = "in", res = dpi, background = "white"
  )
  print(plot)
  grDevices::dev.off()
  tibble(format = names(paths), path = unname(paths))
}

save_top_level_bundle <- function(plot, basename_no_ext, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  paths <- c(
    pdf = file.path(v2_base, "figures", paste0(basename_no_ext, ".pdf")),
    svg = file.path(v2_base, "figures", paste0(basename_no_ext, ".svg")),
    png = file.path(v2_base, "figures", paste0(basename_no_ext, ".png"))
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

contrast_levels <- c("T_vs_A", "T_vs_N")
contrast_labels <- c(T_vs_A = "T vs A", T_vs_N = "T vs N")
all_contrasts <- c("A_vs_N", "T_vs_A", "T_vs_N")
expected_short_de <- c(A_vs_N = 4L, T_vs_A = 4615L, T_vs_N = 4388L)
expected_dte <- c(A_vs_N = 0L, T_vs_A = 4329L, T_vs_N = 4044L)
expected_strict_units <- c(A_vs_N = 0L, T_vs_A = 155L, T_vs_N = 118L)
expected_strict_genes <- c(T_vs_A = 99L, T_vs_N = 74L)
expected_overlap <- c(T_vs_A = 20L, T_vs_N = 16L)
expected_common_genes <- c(T_vs_A = 1503L, T_vs_N = 1503L)
expected_dtu_hashes <- c(
  A_vs_N = "7ac71baab936bc589bcd66ed36e1d956528dcd2c944ed7f685cdeb337325a31c",
  T_vs_A = "60074de59d6ed821a4a8be472405179beabe3aa80dc3fa7a9eaa15eff1f0e0f6",
  T_vs_N = "a3d28b4631b04b2fd0b9eef311a1af4a42d3b537c97f82fc0bb06ccb0a91d6c0"
)
expected_short_de_hash <- "2aca124a299017ff1b52aeb6eb3b22b8479ed179a37668ddbdac14c216974a68"
expected_panel_a_hash <- "1716a23d177a6be6a164ec5715b571ea02d3c78def7bb53823ba2bbfb4b4d89e"
expected_selected_event_hash <- "b2a7b90865ad07004ae6c747a2f37e5537839bc10ef15dcedc5d00a081511d82"

input_paths <- c(
  contract = file.path(framework_root, "validated_data_objects.rds"),
  matched_shortread = file.path(v2_base, "tables", "matched_shortread_primary72_DE_results.tsv"),
  panel_a_counts = file.path(v2_base, "source_data", "Figure3A_matched_shortread_longread_counts.tsv"),
  selected_event_scatter = file.path(v2_base, "source_data", "Figure3D_delta_usage_vs_shortread_gene_log2FC.tsv"),
  dtu_A_vs_N = file.path(phase03_root, "primary72_DTU_isoform_A_vs_N.tsv"),
  dtu_T_vs_A = file.path(phase03_root, "primary72_DTU_isoform_T_vs_A.tsv"),
  dtu_T_vs_N = file.path(phase03_root, "primary72_DTU_isoform_T_vs_N.tsv")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
stop_if(length(missing_inputs) > 0L, paste0("Missing required input(s): ", paste(missing_inputs, collapse = "; ")))

record_progress("INPUT_AUDIT", "Hashing frozen Figure 3 inputs")
input_manifest <- tibble(
  input_role = names(input_paths),
  path = unname(input_paths),
  size_bytes = file.info(unname(input_paths))$size,
  mtime = format(file.info(unname(input_paths))$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  sha256 = vapply(unname(input_paths), sha256_file, character(1)),
  access_mode = "read-only"
)
stop_if(sha256_file(input_paths[["matched_shortread"]]) != expected_short_de_hash,
        "Matched short-read DE source hash mismatch.")
stop_if(sha256_file(input_paths[["panel_a_counts"]]) != expected_panel_a_hash,
        "Frozen Figure 3A count source hash mismatch.")
stop_if(sha256_file(input_paths[["selected_event_scatter"]]) != expected_selected_event_hash,
        "Frozen selected-event scatter source hash mismatch.")
for (z in all_contrasts) {
  stop_if(
    sha256_file(input_paths[[paste0("dtu_", z)]]) != expected_dtu_hashes[[z]],
    paste0("Primary DTU source hash mismatch for ", z, ".")
  )
}
write_tsv2(input_manifest, file.path(dir_results, "Fig3_v2_input_manifest.tsv"))

guard_roots <- c(
  file.path(v2_base, "figures", "main"),
  file.path(v2_base, "source_data"),
  file.path(v2_base, "tables")
)
guard_candidates <- unlist(lapply(
  guard_roots,
  function(d) list.files(d, recursive = TRUE, full.names = TRUE)
))
guard_candidates <- guard_candidates[file.exists(guard_candidates) & !dir.exists(guard_candidates)]
guard_pattern <- "^(Fig|Figure)(1|2|4|5|6|7|8)([^0-9]|$)"
guard_paths <- guard_candidates[grepl(guard_pattern, basename(guard_candidates))]
guard_before <- file_snapshot(guard_paths)

all_new_output_paths <- c(
  file.path(dir_main_pdf, "Figure3_matched_shortread_longread_v2.pdf"),
  file.path(dir_main_svg, "Figure3_matched_shortread_longread_v2.svg"),
  file.path(dir_main_png, "Figure3_matched_shortread_longread_v2.png"),
  file.path(dir_supp_pdf, "Supplementary_Fig_gene_effect_vs_DTU_isoform_delta.pdf"),
  file.path(dir_supp_svg, "Supplementary_Fig_gene_effect_vs_DTU_isoform_delta.svg"),
  file.path(dir_supp_png, "Supplementary_Fig_gene_effect_vs_DTU_isoform_delta.png"),
  file.path(v2_base, "figures", "Supplementary_Fig3_threshold_sensitivity.pdf"),
  file.path(v2_base, "figures", "Supplementary_Fig3_threshold_sensitivity.svg"),
  file.path(v2_base, "figures", "Supplementary_Fig3_threshold_sensitivity.png")
)
archive_path <- archive_existing(all_new_output_paths)

record_progress("CONTRACT", "Loading validated primary72 data contract")
contract <- readRDS(input_paths[["contract"]])
required_objects <- c("lr_counts", "sr_gene_counts", "tx2gene", "meta")
stop_if(!all(required_objects %in% names(contract)), "Validated contract is missing required objects.")
tx2gene <- as.data.frame(contract[["tx2gene"]])
meta <- as.data.frame(contract[["meta"]])
assert_columns(tx2gene, c("isoform_id", "gene_id", "gene_symbol"), "tx2gene")
assert_columns(meta, c("sample_id", "patient_id", "tissue"), "meta")
stop_if(nrow(meta) != 75L, "Validated metadata must contain 75 samples before primary72 exclusion.")

meta$tissue <- as.character(meta$tissue)
crc20_flag <- grepl("^CRC[_-]?20$", meta$patient_id, ignore.case = TRUE)
stop_if(sum(crc20_flag) != 3L, "Could not uniquely identify the three CRC-20 samples.")
primary_meta <- meta[!crc20_flag, , drop = FALSE]
stop_if(nrow(primary_meta) != 72L || length(unique(primary_meta$patient_id)) != 24L,
        "Primary72 cohort is not 24 patients / 72 samples.")
stop_if(!identical(as.integer(table(factor(primary_meta$tissue, levels = c("N", "A", "T")))), c(24L, 24L, 24L)),
        "Primary72 tissue counts are not N/A/T = 24/24/24.")
pair_cells <- primary_meta %>% count(patient_id, tissue)
stop_if(any(pair_cells$n != 1L), "Primary72 contains a duplicated patient-by-tissue sample.")
pair_audit <- pair_cells %>% count(patient_id, name = "n_tissues")
stop_if(any(pair_audit$n_tissues != 3L), "Primary72 patient pairing is incomplete.")

tx_map <- tx2gene %>%
  transmute(
    isoform_id = clean_chr(isoform_id),
    gene_id = clean_chr(gene_id),
    stable_gene_id = strip_ensembl_version(gene_id),
    tx_gene_symbol = clean_chr(gene_symbol)
  )
stop_if(anyNA(tx_map$isoform_id) || anyNA(tx_map$stable_gene_id), "tx2gene contains missing identifiers.")
tx_gene_map <- tx_map %>%
  group_by(stable_gene_id) %>%
  summarise(
    n_versioned_gene_ids = n_distinct(gene_id),
    n_gene_symbols = n_distinct(tx_gene_symbol[!is.na(tx_gene_symbol)]),
    tx_gene_symbol = first_valid(tx_gene_symbol),
    .groups = "drop"
  )

record_progress("SHORT_READ", "Reading frozen matched short-read DESeq2 results")
short_de <- read_tsv(input_paths[["matched_shortread"]], show_col_types = FALSE)
assert_columns(short_de, c(
  "contrast", "stable_gene_id", "gene_symbol", "baseMean", "log2FoldChange", "padj",
  "shortread_DEG", "mapping_status", "DE_threshold", "model"
), "matched short-read DE results")
stop_if(nrow(short_de) != 236058L, "Matched short-read result table row count changed.")
short_duplicates <- short_de %>% count(contrast, stable_gene_id) %>% filter(n != 1L)
stop_if(nrow(short_duplicates) > 0L, "Matched short-read stable Ensembl IDs are duplicated within contrast.")
short_count_check <- short_de %>%
  group_by(contrast) %>%
  summarise(n_DE = sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1), .groups = "drop")
for (z in all_contrasts) {
  observed <- short_count_check$n_DE[short_count_check$contrast == z]
  stop_if(length(observed) != 1L || observed != expected_short_de[[z]],
          paste0("Frozen short-read DEG count mismatch for ", z, "."))
}

panel_a_data <- read_tsv(input_paths[["panel_a_counts"]], show_col_types = FALSE)
assert_columns(panel_a_data, c("contrast", "count", "analysis_type", "statistical_unit", "threshold", "comparison"),
               "Figure 3A frozen counts")
stop_if(nrow(panel_a_data) != 9L, "Figure 3A frozen count table must contain 9 rows.")
for (z in all_contrasts) {
  stop_if(panel_a_data$count[panel_a_data$contrast == z & panel_a_data$analysis_type == "Short-read DEG genes"] != expected_short_de[[z]],
          paste0("Panel A short-read count mismatch for ", z, "."))
  stop_if(panel_a_data$count[panel_a_data$contrast == z & panel_a_data$analysis_type == "Long-read DTE"] != expected_dte[[z]],
          paste0("Panel A DTE count mismatch for ", z, "."))
  stop_if(panel_a_data$count[panel_a_data$contrast == z & panel_a_data$analysis_type == "Long-read strict DTU"] != expected_strict_units[[z]],
          paste0("Panel A strict DTU count mismatch for ", z, "."))
}

record_progress("DTU", "Reading frozen primary72 satuRn + stageR tables")
dtu_list <- setNames(lapply(all_contrasts, function(z) {
  x <- read_tsv(input_paths[[paste0("dtu_", z)]], show_col_types = FALSE)
  assert_columns(x, c(
    "gene_id", "gene_name", "isoform_id", "contrast", "analysis_cohort",
    "mean_usage_N", "mean_usage_A", "mean_usage_T", "delta_usage", "padj"
  ), paste0("DTU ", z))
  stop_if(nrow(x) != 5477L, paste0("DTU tested isoform count changed for ", z, "."))
  stop_if(any(x$contrast != z) || any(x$analysis_cohort != "primary72"),
          paste0("DTU contrast/cohort label mismatch for ", z, "."))
  x %>% mutate(
    gene_id = clean_chr(gene_id),
    stable_gene_id = strip_ensembl_version(gene_id),
    gene_symbol_dtu = clean_chr(gene_name),
    isoform_id = clean_chr(isoform_id),
    strict_primary_call = !is.na(padj) & padj < 0.10 & abs(delta_usage) >= 0.05
  )
}), all_contrasts)

for (z in all_contrasts) {
  x <- dtu_list[[z]]
  observed_units <- sum(x$strict_primary_call)
  stop_if(observed_units != expected_strict_units[[z]], paste0("Strict DTU unit gate failed for ", z, "."))
  if (z %in% contrast_levels) {
    observed_genes <- n_distinct(x$stable_gene_id[x$strict_primary_call])
    stop_if(observed_genes != expected_strict_genes[[z]], paste0("Strict DTU gene gate failed for ", z, "."))
  }
}

record_progress("COMMON_UNIVERSE", "Building one common testable multi-isoform gene universe per tumor contrast")
common_list <- list()
mapping_audit_list <- list()

for (target_contrast in contrast_levels) {
  dtu <- dtu_list[[target_contrast]]
  expected_delta <- if (target_contrast == "T_vs_A") {
    dtu$mean_usage_T - dtu$mean_usage_A
  } else {
    dtu$mean_usage_T - dtu$mean_usage_N
  }
  stop_if(any(!is.finite(expected_delta)) || max(abs(expected_delta - dtu$delta_usage)) > 1e-10,
          paste0("Frozen delta usage is inconsistent with group mean usage for ", target_contrast, "."))
  dtu$delta_for_tvd <- expected_delta

  usage_gene <- dtu %>%
    group_by(stable_gene_id) %>%
    summarise(
      gene_id = first_valid(gene_id),
      gene_symbol_dtu = first_valid(gene_symbol_dtu),
      n_tested_isoforms = n(),
      sum_usage_N = sum(mean_usage_N),
      sum_usage_A = sum(mean_usage_A),
      sum_usage_T = sum(mean_usage_T),
      TVD_usage_change = 0.5 * sum(abs(delta_for_tvd)),
      max_abs_delta_usage = max(abs(delta_for_tvd)),
      longread_DTU_status = any(strict_primary_call),
      .groups = "drop"
    )
  stop_if(nrow(usage_gene) != 1503L || any(usage_gene$n_tested_isoforms < 2L),
          paste0("Formal multi-isoform tested universe changed for ", target_contrast, "."))
  stop_if(any(abs(usage_gene$sum_usage_N - 1) > 1e-6) ||
            any(abs(usage_gene$sum_usage_A - 1) > 1e-6) ||
            any(abs(usage_gene$sum_usage_T - 1) > 1e-6),
          paste0("Within-gene usage does not close to one for ", target_contrast, "."))
  stop_if(any(usage_gene$TVD_usage_change < -1e-12 | usage_gene$TVD_usage_change > 1 + 1e-12),
          paste0("TVD is outside [0,1] for ", target_contrast, "."))

  short_one <- short_de %>%
    filter(.data$contrast == .env$target_contrast) %>%
    transmute(
      stable_gene_id,
      short_gene_symbol = clean_chr(gene_symbol),
      shortread_log2FC = as.numeric(log2FoldChange),
      shortread_padj = as.numeric(padj),
      shortread_baseMean = as.numeric(baseMean),
      shortread_DE_status = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
      short_mapping_status = clean_chr(mapping_status)
    )

  audit <- usage_gene %>%
    left_join(short_one, by = "stable_gene_id") %>%
    left_join(tx_gene_map, by = "stable_gene_id") %>%
    mutate(
      contrast = target_contrast,
      has_shortread_result = !is.na(short_mapping_status),
      unique_tx_gene_mapping = !is.na(n_versioned_gene_ids) & n_versioned_gene_ids == 1L,
      shared_unambiguous_Ensembl = short_mapping_status == "shared_unambiguous_Ensembl",
      valid_shortread_effect = is.finite(shortread_log2FC) & is.finite(shortread_baseMean),
      valid_usage_change = is.finite(TVD_usage_change) & is.finite(max_abs_delta_usage),
      multi_isoform_tested = n_tested_isoforms >= 2L,
      included = has_shortread_result & unique_tx_gene_mapping & shared_unambiguous_Ensembl &
        valid_shortread_effect & valid_usage_change & multi_isoform_tested,
      exclusion_reason = case_when(
        !has_shortread_result ~ "missing_shortread_result",
        !unique_tx_gene_mapping ~ "ambiguous_longread_tx2gene_mapping",
        !shared_unambiguous_Ensembl ~ "not_shared_unambiguous_Ensembl",
        !valid_shortread_effect ~ "invalid_shortread_gene_effect",
        !valid_usage_change ~ "invalid_longread_usage_change",
        !multi_isoform_tested ~ "fewer_than_two_tested_isoforms",
        TRUE ~ "included"
      )
    )
  mapping_audit_list[[target_contrast]] <- audit

  common <- audit %>%
    filter(included) %>%
    transmute(
      contrast,
      gene_id,
      gene_symbol = coalesce(short_gene_symbol, tx_gene_symbol, gene_symbol_dtu),
      stable_gene_id,
      shortread_log2FC,
      shortread_padj,
      shortread_baseMean,
      shortread_DE_status,
      longread_DTU_status,
      n_tested_isoforms,
      TVD_usage_change,
      max_abs_delta_usage
    ) %>%
    arrange(stable_gene_id)

  stop_if(nrow(common) != expected_common_genes[[target_contrast]],
          paste0("Common testable gene count gate failed for ", target_contrast, "."))
  stop_if(n_distinct(common$stable_gene_id) != nrow(common),
          paste0("Common universe contains duplicate stable gene IDs for ", target_contrast, "."))
  stop_if(sum(common$longread_DTU_status) != expected_strict_genes[[target_contrast]],
          paste0("Strict DTU gene count changed after common-universe restriction for ", target_contrast, "."))
  stop_if(sum(common$longread_DTU_status & common$shortread_DE_status) != expected_overlap[[target_contrast]],
          paste0("Known DTU/DE overlap gate failed for ", target_contrast, "."))
  common_list[[target_contrast]] <- common
}

common_universe <- bind_rows(common_list)
mapping_audit <- bind_rows(mapping_audit_list) %>%
  select(
    contrast, stable_gene_id, gene_id, n_tested_isoforms,
    has_shortread_result, unique_tx_gene_mapping, shared_unambiguous_Ensembl,
    valid_shortread_effect, valid_usage_change, multi_isoform_tested,
    included, exclusion_reason
  )
stop_if(any(!mapping_audit$included),
        "At least one formally tested long-read gene failed the common-universe mapping criteria.")

common_summary <- common_universe %>%
  group_by(contrast) %>%
  summarise(
    n_common_testable_genes = n(),
    n_shortread_DE = sum(shortread_DE_status),
    background_DE_fraction = mean(shortread_DE_status),
    n_shortread_nonDE = sum(!shortread_DE_status),
    background_nonDE_fraction = mean(!shortread_DE_status),
    n_strict_DTU_genes = sum(longread_DTU_status),
    n_DTU_and_DE = sum(longread_DTU_status & shortread_DE_status),
    n_DTU_only = sum(longread_DTU_status & !shortread_DE_status),
    .groups = "drop"
  )

write_tsv2(common_universe, file.path(dir_results, "Fig3_common_testable_gene_universe.tsv"))
write_tsv2(common_summary, file.path(dir_results, "Fig3_common_testable_gene_universe_summary.tsv"))
write_tsv2(mapping_audit, file.path(dir_results, "Fig3_common_universe_mapping_audit.tsv"))

record_progress("FISHER", "Computing two-sided Fisher tests within the fixed common universe")
fisher_rows <- lapply(contrast_levels, function(target_contrast) {
  x <- common_list[[target_contrast]]
  a <- sum(x$longread_DTU_status & x$shortread_DE_status)
  b <- sum(x$longread_DTU_status & !x$shortread_DE_status)
  c0 <- sum(!x$longread_DTU_status & x$shortread_DE_status)
  d0 <- sum(!x$longread_DTU_status & !x$shortread_DE_status)
  test <- fisher.test(matrix(c(a, b, c0, d0), nrow = 2L, byrow = TRUE), alternative = "two.sided")
  direction <- if (test$p.value < 0.05 && unname(test$estimate) > 1) {
    "enrichment"
  } else if (test$p.value < 0.05 && unname(test$estimate) < 1) {
    "depletion"
  } else {
    "no clear enrichment"
  }
  tibble(
    contrast = target_contrast,
    n_common_testable_genes = nrow(x),
    a_DTU_and_DE = a,
    b_DTU_and_nonDE = b,
    c_nonDTU_and_DE = c0,
    d_nonDTU_and_nonDE = d0,
    n_shortread_DE = a + c0,
    n_strict_DTU_genes = a + b,
    n_DTU_only = b,
    n_shortread_DE_only = c0,
    all_common_DE_fraction = (a + c0) / nrow(x),
    all_common_nonDE_fraction = (b + d0) / nrow(x),
    nonDTU_DE_fraction = c0 / (c0 + d0),
    DTU_gene_DE_fraction = a / (a + b),
    DTU_gene_nonDE_fraction = b / (a + b),
    odds_ratio = unname(test$estimate),
    odds_ratio_CI_low = test$conf.int[[1]],
    odds_ratio_CI_high = test$conf.int[[2]],
    Fisher_two_sided_P = test$p.value,
    direction = direction
  )
})
fisher_summary <- bind_rows(fisher_rows)
write_tsv2(fisher_summary, file.path(dir_results, "Fig3_DTU_DE_contingency_fisher.tsv"))

bootstrap_spearman <- function(x, y, n_boot = 4000L) {
  stop_if(length(x) != length(y) || length(x) < 3L, "Invalid input for Spearman bootstrap.")
  n <- length(x)
  values <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    values[[i]] <- suppressWarnings(cor(x[idx], y[idx], method = "spearman", use = "complete.obs"))
  }
  tibble(
    bootstrap_replicates = n_boot,
    bootstrap_valid = sum(is.finite(values)),
    CI_low = unname(quantile(values, 0.025, na.rm = TRUE)),
    CI_high = unname(quantile(values, 0.975, na.rm = TRUE))
  )
}

record_progress("CONTINUOUS", "Computing all-gene TVD and max-delta Spearman analyses with 4000 gene bootstraps")
coupling_data <- common_universe %>%
  mutate(abs_shortread_log2FC = abs(shortread_log2FC))
coupling_summary_list <- list()
set.seed(20260828)
for (target_contrast in contrast_levels) {
  x <- coupling_data %>% filter(.data$contrast == .env$target_contrast)
  for (metric in c("TVD_usage_change", "max_abs_delta_usage")) {
    y <- x[[metric]]
    test <- suppressWarnings(cor.test(x$abs_shortread_log2FC, y, method = "spearman", exact = FALSE))
    ci <- bootstrap_spearman(x$abs_shortread_log2FC, y, n_boot = 4000L)
    coupling_summary_list[[paste(target_contrast, metric, sep = "__")]] <- tibble(
      contrast = target_contrast,
      metric = metric,
      n_genes = nrow(x),
      Spearman_rho = unname(test$estimate),
      two_sided_P = test$p.value,
      bootstrap_CI_low = ci$CI_low,
      bootstrap_CI_high = ci$CI_high,
      bootstrap_replicates = ci$bootstrap_replicates,
      bootstrap_valid = ci$bootstrap_valid,
      bootstrap_unit = "gene"
    )
  }
}
coupling_summary <- bind_rows(coupling_summary_list)
write_tsv2(coupling_data, file.path(dir_results, "Fig3_continuous_gene_isoform_coupling.tsv"))
write_tsv2(coupling_summary, file.path(dir_results, "Fig3_continuous_gene_isoform_coupling_summary.tsv"))

record_progress("THRESHOLDS", "Evaluating the complete prespecified DTU/DE threshold grid")
threshold_grid <- expand_grid(
  DTU_FDR_threshold = c(0.05, 0.10),
  DTU_abs_delta_threshold = c(0.03, 0.05, 0.10),
  shortread_padj_threshold = c(0.01, 0.05, 0.10),
  shortread_abs_log2FC_threshold = c(0.5, 1.0)
)
threshold_rows <- list()
row_index <- 1L
for (target_contrast in contrast_levels) {
  common_ids <- common_list[[target_contrast]]$stable_gene_id
  dtu <- dtu_list[[target_contrast]] %>% filter(stable_gene_id %in% common_ids)
  short_one <- short_de %>%
    filter(.data$contrast == .env$target_contrast, stable_gene_id %in% common_ids)
  for (i in seq_len(nrow(threshold_grid))) {
    pars <- threshold_grid[i, ]
    dtu_genes <- dtu %>%
      group_by(stable_gene_id) %>%
      summarise(
        called = any(
          !is.na(padj) & padj < pars$DTU_FDR_threshold &
            abs(delta_usage) >= pars$DTU_abs_delta_threshold
        ),
        .groups = "drop"
      ) %>%
      filter(called) %>%
      pull(stable_gene_id)
    de_genes <- short_one %>%
      filter(
        !is.na(padj),
        padj < pars$shortread_padj_threshold,
        abs(log2FoldChange) > pars$shortread_abs_log2FC_threshold
      ) %>%
      pull(stable_gene_id)
    n_dtu <- length(dtu_genes)
    n_overlap <- sum(dtu_genes %in% de_genes)
    threshold_rows[[row_index]] <- tibble(
      contrast = target_contrast,
      DTU_FDR_threshold = pars$DTU_FDR_threshold,
      DTU_abs_delta_threshold = pars$DTU_abs_delta_threshold,
      shortread_padj_threshold = pars$shortread_padj_threshold,
      shortread_abs_log2FC_threshold = pars$shortread_abs_log2FC_threshold,
      n_common_testable_genes = length(common_ids),
      n_DTU_genes = n_dtu,
      n_shortread_DE_genes = length(de_genes),
      n_overlap = n_overlap,
      n_DTU_nonDE = n_dtu - n_overlap,
      fraction_nonDE = if (n_dtu > 0L) (n_dtu - n_overlap) / n_dtu else NA_real_,
      primary_definition =
        pars$DTU_FDR_threshold == 0.10 & pars$DTU_abs_delta_threshold == 0.05 &
        pars$shortread_padj_threshold == 0.05 & pars$shortread_abs_log2FC_threshold == 1.0
    )
    row_index <- row_index + 1L
  }
}
threshold_sensitivity <- bind_rows(threshold_rows)
stop_if(nrow(threshold_sensitivity) != 72L, "Threshold sensitivity grid must contain 72 rows.")
primary_threshold_rows <- threshold_sensitivity %>% filter(primary_definition)
stop_if(nrow(primary_threshold_rows) != 2L,
        "Threshold sensitivity grid does not contain exactly two primary rows.")
for (target_contrast in contrast_levels) {
  primary_row <- primary_threshold_rows %>% filter(.data$contrast == .env$target_contrast)
  expected_summary <- fisher_summary %>% filter(.data$contrast == .env$target_contrast)
  stop_if(
    primary_row$n_DTU_genes != expected_summary$n_strict_DTU_genes ||
      primary_row$n_overlap != expected_summary$a_DTU_and_DE ||
      abs(primary_row$fraction_nonDE - expected_summary$DTU_gene_nonDE_fraction) > 1e-12,
    paste0("Primary threshold sensitivity cell mismatch for ", target_contrast, ".")
  )
}
write_tsv2(threshold_sensitivity, file.path(dir_results, "Fig3_threshold_sensitivity.tsv"))

selected_event_data <- read_tsv(input_paths[["selected_event_scatter"]], show_col_types = FALSE)
assert_columns(selected_event_data, c(
  "contrast", "gene_id", "stable_gene_id", "gene_symbol", "isoform_id",
  "delta_usage", "transcript_adjusted_FDR", "shortread_log2FC", "shortread_padj",
  "shortread_DEG", "comparison"
), "Selected-event supplementary source")
stop_if(nrow(selected_event_data) != 273L ||
          sum(selected_event_data$contrast == "T_vs_A") != 155L ||
          sum(selected_event_data$contrast == "T_vs_N") != 118L,
        "Selected-event supplementary source no longer contains 155/118 strict DTU calls.")

source_a <- panel_a_data %>% arrange(factor(analysis_type, levels = c(
  "Short-read DEG genes", "Long-read DTE", "Long-read strict DTU"
)), factor(contrast, levels = all_contrasts))
source_b <- fisher_summary
source_c <- fisher_summary %>%
  transmute(
    contrast,
    n_DTU_nonDE = b_DTU_and_nonDE,
    n_strict_DTU_genes,
    DTU_nonDE_fraction = DTU_gene_nonDE_fraction,
    n_all_common_nonDE = b_DTU_and_nonDE + d_nonDTU_and_nonDE,
    n_common_testable_genes,
    all_common_nonDE_fraction
  )
source_d <- coupling_data

source_paths <- c(
  A = write_tsv2(source_a, file.path(dir_source, "SourceData_Figure3A_v2.tsv")),
  B = write_tsv2(source_b, file.path(dir_source, "SourceData_Figure3B_v2.tsv")),
  C = write_tsv2(source_c, file.path(dir_source, "SourceData_Figure3C_v2.tsv")),
  D = write_tsv2(source_d, file.path(dir_source, "SourceData_Figure3D_v2.tsv")),
  selected_event = write_tsv2(selected_event_data, file.path(dir_source, "SourceData_Supplementary_selected_DTU_events.tsv")),
  threshold = write_tsv2(threshold_sensitivity, file.path(dir_source, "SourceData_Supplementary_threshold_sensitivity.tsv"))
)

record_progress("PLOT_A", "Rendering the frozen cross-method detection landscape")
analysis_levels <- c("Short-read DEG genes", "Long-read DTE", "Long-read strict DTU")
analysis_labels <- c(
  "Short-read DEGs\n(gene)",
  "Long-read DTE\n(isoform)",
  "Long-read strict DTU\n(DTU unit)"
)
names(analysis_labels) <- analysis_levels
analysis_colors <- c(
  "Short-read DEG genes" = "#356D87",
  "Long-read DTE" = "#D97921",
  "Long-read strict DTU" = "#7B4EA3"
)
panel_a_plot_data <- source_a %>%
  mutate(
    analysis_type = factor(analysis_type, levels = analysis_levels),
    comparison = factor(comparison, levels = c("A vs N", "T vs A", "T vs N"))
  )
panel_a <- ggplot(panel_a_plot_data, aes(comparison, count, colour = analysis_type)) +
  geom_segment(aes(xend = comparison, y = 0, yend = count), linewidth = 0.75, alpha = 0.78) +
  geom_point(size = 3.0) +
  geom_text(
    aes(label = comma(count)), colour = "#202020", fontface = "bold",
    vjust = -0.7, size = 2.55
  ) +
  facet_wrap(~ analysis_type, scales = "free_y", nrow = 1,
             labeller = as_labeller(analysis_labels)) +
  scale_colour_manual(values = analysis_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = "A  Cross-method detection landscape",
    subtitle = "Counts summarize different statistical units and are not direct platform-performance measures",
    x = NULL, y = "Count within stated unit"
  ) +
  theme_manuscript(7.2) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.spacing.x = grid::unit(5, "mm")
  )

circle_df <- function(cx, cy, radius, group_name) {
  theta <- seq(0, 2 * pi, length.out = 241)
  tibble(
    x = cx + radius * cos(theta),
    y = cy + radius * sin(theta),
    group = group_name
  )
}

format_p_relation <- function(x, digits = 3, eps = 1e-4) {
  vapply(x, function(value) {
    if (is.na(value)) return("= NA")
    if (value < eps) return(sprintf("< %.4f", eps))
    paste("=", format.pval(value, digits = digits, eps = 0))
  }, character(1))
}

make_overlap_plot <- function(row) {
  dt <- circle_df(0, 0, 0.78, "DTU")
  de <- circle_df(1.06, 0, 0.78, "DE")
  circles <- bind_rows(dt, de)
  common_de_pct <- 100 * row$all_common_DE_fraction
  dtu_de_pct <- 100 * row$DTU_gene_DE_fraction
  non_dtu_de_pct <- 100 * row$nonDTU_DE_fraction
  metrics <- c(
    sprintf("Common testable genes = %s", comma(row$n_common_testable_genes)),
    sprintf("All-common DEG rate = %.1f%%", common_de_pct),
    sprintf("DTU-gene DEG rate = %.1f%%", dtu_de_pct),
    sprintf("Non-DTU DEG rate = %.1f%%", non_dtu_de_pct),
    sprintf("OR = %.2f (95%% CI %.2f-%.2f)", row$odds_ratio, row$odds_ratio_CI_low, row$odds_ratio_CI_high),
    sprintf("Fisher P %s", format_p_relation(row$Fisher_two_sided_P, digits = 2))
  )
  ggplot(circles, aes(x, y, group = group, fill = group, colour = group)) +
    geom_polygon(alpha = 0.16, linewidth = 0.6) +
    geom_path(linewidth = 0.55) +
    annotate("text", x = 0.30, y = 1.00,
             label = sprintf("Strict DTU\n(n=%d)", row$n_strict_DTU_genes),
             colour = "#C75B12", fontface = "bold", hjust = 1, size = 2.05) +
    annotate("text", x = 0.76, y = 1.00,
             label = sprintf("Short-read DEG\n(n=%d)", row$n_shortread_DE),
             colour = "#245D82", fontface = "bold", hjust = 0, size = 2.05) +
    annotate("text", x = -0.42, y = 0,
             label = sprintf("DTU only\n%d", row$n_DTU_only), fontface = "bold", size = 1.90) +
    annotate("text", x = 0.53, y = 0,
             label = sprintf("Both\n%d", row$a_DTU_and_DE), fontface = "bold", size = 1.90) +
    annotate("text", x = 1.48, y = 0,
             label = sprintf("DE only\n%d", row$n_shortread_DE_only), fontface = "bold", size = 1.90) +
    annotate("text", x = 2.02, y = 0.88, label = metrics[[1]], hjust = 0, size = 2.05) +
    annotate("text", x = 2.02, y = 0.54, label = metrics[[2]], hjust = 0, size = 2.05) +
    annotate("text", x = 2.02, y = 0.20, label = metrics[[3]], hjust = 0, size = 2.05) +
    annotate("text", x = 2.02, y = -0.14, label = metrics[[4]], hjust = 0, size = 2.05) +
    annotate("text", x = 2.02, y = -0.48, label = metrics[[5]], hjust = 0, size = 2.05) +
    annotate("text", x = 2.02, y = -0.82, label = metrics[[6]], hjust = 0, size = 2.05) +
    annotate("text", x = 0.53, y = -1.05, label = "Schematic; areas are not proportional",
             colour = "#666666", size = 1.85) +
    scale_fill_manual(values = c(DTU = "#E66B16", DE = "#3A7095"), guide = "none") +
    scale_colour_manual(values = c(DTU = "#E66B16", DE = "#3A7095"), guide = "none") +
    coord_fixed(xlim = c(-0.90, 5.0), ylim = c(-1.22, 1.24), clip = "off") +
    labs(title = contrast_labels[[row$contrast]]) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", size = 8.0, hjust = 0.18, margin = margin(b = 2)),
      plot.margin = margin(2, 2, 2, 2)
    )
}

record_progress("PLOT_B", "Rendering non-proportional common-universe overlap schematics")
panel_b_plots <- lapply(contrast_levels, function(z) {
  make_overlap_plot(fisher_summary %>% filter(.data$contrast == .env$z))
})
panel_b <- wrap_plots(panel_b_plots, nrow = 1) +
  plot_annotation(
    title = "B  DTU-DE overlap in the common testable gene universe",
    theme = theme(
      plot.title = element_text(family = FONT_FAMILY, face = "bold", size = 8.2, margin = margin(b = 2))
    )
  )

record_progress("PLOT_C", "Rendering DTU non-DE fractions with all-common background references")
panel_c_data <- source_c %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    main_percent = 100 * DTU_nonDE_fraction,
    background_percent = 100 * all_common_nonDE_fraction,
    main_label = sprintf("%.1f%%\n%d/%d", main_percent, n_DTU_nonDE, n_strict_DTU_genes),
    background_label = sprintf("All common\n%.1f%%", background_percent)
  )
panel_c <- ggplot(panel_c_data, aes(comparison, main_percent, colour = comparison)) +
  geom_segment(aes(xend = comparison, y = 0, yend = main_percent), linewidth = 1.0, alpha = 0.82) +
  geom_point(size = 4.0) +
  geom_text(aes(label = main_label), colour = "#202020", fontface = "bold", vjust = 1.65, size = 2.45) +
  geom_segment(
    aes(x = as.numeric(comparison) - 0.28, xend = as.numeric(comparison) + 0.28,
        y = background_percent, yend = background_percent),
    inherit.aes = FALSE, colour = "#777777", linetype = "dashed", linewidth = 0.7
  ) +
  geom_text(
    aes(x = as.numeric(comparison), y = background_percent, label = background_label),
    inherit.aes = FALSE, colour = "#555555", vjust = -0.45, size = 2.15
  ) +
  scale_colour_manual(values = c("T vs A" = "#E66B16", "T vs N" = "#7B4EA3"), guide = "none") +
  scale_y_continuous(limits = c(0, 104), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title = "C  Strict DTU genes not meeting short-read DEG criteria",
    subtitle = "Grey reference: all-common-gene non-DE fraction",
    x = NULL, y = "Non-DE fraction"
  ) +
  theme_manuscript(7.0)

record_progress("PLOT_D", "Rendering all-gene abundance/TVD coupling")
panel_d_data <- coupling_data %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    highlight = ifelse(longread_DTU_status, "Strict DTU gene", "Other common gene")
  )
panel_d_ann <- coupling_summary %>%
  filter(metric == "TVD_usage_change") %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    label = sprintf(
      "Spearman rho = %.3f\n95%% CI %.3f-%.3f\nP %s; n = %s",
      Spearman_rho, bootstrap_CI_low, bootstrap_CI_high,
      format_p_relation(two_sided_P, digits = 2), comma(n_genes)
    ),
    x = 5.72,
    y = 0.302
  )
panel_d <- ggplot(panel_d_data, aes(abs_shortread_log2FC, TVD_usage_change)) +
  geom_point(
    data = panel_d_data %>% filter(!longread_DTU_status),
    colour = "#AEB4B8", alpha = 0.34, size = 0.72
  ) +
  geom_point(
    data = panel_d_data %>% filter(longread_DTU_status),
    aes(fill = comparison), shape = 21, colour = "white", stroke = 0.18,
    alpha = 0.90, size = 1.8
  ) +
  geom_text(
    data = panel_d_ann,
    aes(x = x, y = y, label = label), inherit.aes = FALSE,
    hjust = 1, vjust = 1, size = 2.25, lineheight = 1.05
  ) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = c("T vs A" = "#E66B16", "T vs N" = "#7B4EA3"), guide = "none") +
  scale_x_continuous(limits = c(0, 5.95), breaks = 0:5) +
  scale_y_continuous(limits = c(0, 0.315), breaks = seq(0, 0.3, by = 0.1)) +
  labs(
    title = "D  Continuous coupling of gene abundance and isoform composition",
    subtitle = "Statistics use all common genes; strict DTU genes are highlighted only for visualization",
    x = "|Matched short-read gene log2FC|",
    y = "Within-gene isoform composition change (TVD)"
  ) +
  theme_manuscript(7.2) +
  theme(panel.grid.major.x = element_line(linewidth = 0.15, colour = "#EFEFEF"))

panel_b_nested <- wrap_elements(full = panel_b)
middle_row <- (panel_b_nested | panel_c) + plot_layout(widths = c(1.92, 0.92))
main_figure <- panel_a / middle_row / panel_d +
  plot_layout(heights = c(1.0, 1.42, 1.35)) +
  plot_annotation(
    title = "Figure 3 | Matched short-read gene expression contextualizes long-read isoform remodeling",
    subtitle = "Paired primary72 cohort; common-universe tests use 1,503 multi-isoform genes per tumor contrast",
    caption = paste0(
      "Primary thresholds: short-read DEG adjusted P < 0.05 and |log2FC| > 1; ",
      "long-read strict DTU transcript-adjusted FDR < 0.10 and |delta usage| >= 0.05."
    ),
    theme = theme(
      plot.title = element_text(family = FONT_FAMILY, face = "bold", size = 11.6, margin = margin(b = 2)),
      plot.subtitle = element_text(family = FONT_FAMILY, size = 8.1, colour = "#555555", margin = margin(b = 6)),
      plot.caption = element_text(family = FONT_FAMILY, size = 6.2, colour = "#555555", hjust = 0, margin = margin(t = 5)),
      plot.margin = margin(5, 5, 4, 5)
    )
  )

record_progress("RENDER_MAIN", "Exporting the 180-mm Figure 3 v2 bundle")
main_paths <- save_bundle(
  main_figure, "Figure3_matched_shortread_longread_v2",
  width_mm = 180, height_mm = 184, supplementary = FALSE, dpi = 600
)

record_progress("RENDER_SELECTED", "Exporting the selected-event descriptive scatter as a supplementary figure")
selected_plot_data <- selected_event_data %>%
  mutate(
    comparison = factor(comparison, levels = c("T vs A", "T vs N")),
    shortread_status = ifelse(shortread_DEG, "Short-read DE", "Not short-read DE")
  )
selected_plot <- ggplot(selected_plot_data, aes(shortread_log2FC, delta_usage, colour = shortread_status)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", linewidth = 0.3, colour = "#777777") +
  geom_point(size = 1.45, alpha = 0.82) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_colour_manual(values = c("Not short-read DE" = "#AEB4B8", "Short-read DE" = "#C6534C")) +
  labs(
    title = "Isoform-level usage shifts among strict DTU events in relation to matched gene-level expression changes",
    subtitle = "Selected-event descriptive analysis; this is not an all-gene correlation test",
    x = "Matched short-read gene log2FC",
    y = "Long-read isoform delta usage",
    colour = NULL,
    caption = "Each point is one strict isoform-level DTU call (T vs A n=155; T vs N n=118)."
  ) +
  theme_manuscript(7.4) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linewidth = 0.15, colour = "#EFEFEF")
  )
selected_paths <- save_bundle(
  selected_plot, "Supplementary_Fig_gene_effect_vs_DTU_isoform_delta",
  width_mm = 180, height_mm = 86, supplementary = TRUE, dpi = 600
)

record_progress("RENDER_SENSITIVITY", "Exporting the threshold-sensitivity heatmap")
heatmap_data <- threshold_sensitivity %>%
  mutate(
    comparison = factor(unname(contrast_labels[contrast]), levels = c("T vs A", "T vs N")),
    dtu_label = sprintf("FDR<%.2f | delta>=%.2f", DTU_FDR_threshold, DTU_abs_delta_threshold),
    de_label = sprintf("P<%.2f\n|LFC|>%.1f", shortread_padj_threshold, shortread_abs_log2FC_threshold),
    cell_label = sprintf("%.1f%%", 100 * fraction_nonDE)
  )
dtu_order <- c(
  "FDR<0.05 | delta>=0.03", "FDR<0.05 | delta>=0.05", "FDR<0.05 | delta>=0.10",
  "FDR<0.10 | delta>=0.03", "FDR<0.10 | delta>=0.05", "FDR<0.10 | delta>=0.10"
)
de_order <- c(
  "P<0.01\n|LFC|>0.5", "P<0.01\n|LFC|>1.0",
  "P<0.05\n|LFC|>0.5", "P<0.05\n|LFC|>1.0",
  "P<0.10\n|LFC|>0.5", "P<0.10\n|LFC|>1.0"
)
heatmap_data$dtu_label <- factor(heatmap_data$dtu_label, levels = rev(dtu_order))
heatmap_data$de_label <- factor(heatmap_data$de_label, levels = de_order)
threshold_plot <- ggplot(heatmap_data, aes(de_label, dtu_label, fill = fraction_nonDE)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_tile(
    data = heatmap_data %>% filter(primary_definition),
    fill = NA, colour = "#111111", linewidth = 0.9
  ) +
  geom_text(aes(label = cell_label), size = 2.25, colour = "#202020") +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_gradientn(
    colours = c("#F3F4F4", "#F6C58E", "#D85F18"),
    limits = c(0.35, 0.85), oob = squish,
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Supplementary Figure | Threshold sensitivity of the DTU-gene non-DE fraction",
    subtitle = "Outlined cells are the frozen primary definition; primary thresholds were not changed",
    x = "Matched short-read DEG criterion",
    y = "Long-read DTU criterion",
    fill = "DTU genes\nnot short-read DE"
  ) +
  theme_manuscript(7.0) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )
threshold_paths <- save_top_level_bundle(
  threshold_plot, "Supplementary_Fig3_threshold_sensitivity",
  width_mm = 180, height_mm = 112, dpi = 600
)

record_progress("REPORTS", "Writing source-data indexes and Chinese manuscript text")
source_index <- tibble(
  figure = c("Figure 3", "Figure 3", "Figure 3", "Figure 3",
             "Supplementary selected-event scatter", "Supplementary threshold sensitivity"),
  panel = c("A", "B", "C", "D", "all", "all"),
  source_data_path = unname(source_paths),
  n_rows = vapply(unname(source_paths), count_tsv_rows, integer(1)),
  description = c(
    "Frozen cross-method detection counts and statistical units",
    "Common-universe 2x2 contingency tables, Fisher odds ratios and rates",
    "Strict-DTU non-DE fractions and all-common-gene background references",
    "All common eligible genes with short-read effect, TVD and max-delta metrics",
    "Strict isoform-level DTU calls used only for selected-event descriptive visualization",
    "Complete prespecified DTU/short-read threshold grid"
  )
)
source_index_path <- write_tsv2(source_index, file.path(dir_results, "Fig3_v2_source_data_index.tsv"))

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt_p <- function(x) format_p_relation(x, digits = 3, eps = 1e-4)
ta <- fisher_summary %>% filter(contrast == "T_vs_A")
tn <- fisher_summary %>% filter(contrast == "T_vs_N")
ta_tvd <- coupling_summary %>% filter(contrast == "T_vs_A", metric == "TVD_usage_change")
tn_tvd <- coupling_summary %>% filter(contrast == "T_vs_N", metric == "TVD_usage_change")
ta_max <- coupling_summary %>% filter(contrast == "T_vs_A", metric == "max_abs_delta_usage")
tn_max <- coupling_summary %>% filter(contrast == "T_vs_N", metric == "max_abs_delta_usage")

threshold_range <- threshold_sensitivity %>%
  group_by(contrast) %>%
  summarise(
    min_fraction = min(fraction_nonDE, na.rm = TRUE),
    max_fraction = max(fraction_nonDE, na.rm = TRUE),
    min_n_DTU = min(n_DTU_genes),
    max_n_DTU = max(n_DTU_genes),
    .groups = "drop"
  )
threshold_primary_dtu_range <- threshold_sensitivity %>%
  filter(DTU_FDR_threshold == 0.10, DTU_abs_delta_threshold == 0.05) %>%
  group_by(contrast) %>%
  summarise(
    min_fraction = min(fraction_nonDE),
    max_fraction = max(fraction_nonDE),
    .groups = "drop"
  )

legend_zh <- c(
  "# Figure 3 v2 中文图注",
  "",
  "## Figure 3 | Matched short-read gene expression contextualizes long-read isoform remodeling",
  "",
  "分析使用 primary72 队列，即24例具有完整 N/A/T 配对组织的患者；长、短读长文库来自同一 RNA 提取物分装。",
  "",
  "**A，** A vs N、T vs A 和 T vs N 中 matched short-read DEG genes、long-read DTE isoforms 和 satuRn + stageR strict DTU units 的冻结数量。三类结果的统计单位不同，因此不能用绝对数量比较测序平台性能。",
  "",
  sprintf("**B，** 在每个肿瘤 contrast 的 %s 个 common testable multi-isoform genes 中，展示 strict DTU 与 matched short-read DEG 的集合关系。T vs A 中 DTU only/overlap/DE only 为 %d/%d/%d；T vs N 为 %d/%d/%d。圆形为 non-proportional schematic。双侧 Fisher exact test 以 DTU/non-DTU 和 DEG/non-DEG 构成2x2表。T vs A: OR=%.2f, 95%% CI %.2f-%.2f, P %s；T vs N: OR=%.2f, 95%% CI %.2f-%.2f, P %s。",
          comma(ta$n_common_testable_genes), ta$n_DTU_only, ta$a_DTU_and_DE, ta$n_shortread_DE_only,
          tn$n_DTU_only, tn$a_DTU_and_DE, tn$n_shortread_DE_only,
          ta$odds_ratio, ta$odds_ratio_CI_low, ta$odds_ratio_CI_high, fmt_p(ta$Fisher_two_sided_P),
          tn$odds_ratio, tn$odds_ratio_CI_low, tn$odds_ratio_CI_high, fmt_p(tn$Fisher_two_sided_P)),
  "",
  sprintf("**C，** strict DTU genes 中未达到 matched short-read DEG 标准的比例：T vs A 为 %s (%d/%d)，T vs N 为 %s (%d/%d)。灰色参照线为同一 common universe 中全部 genes 的 non-DE 比例，分别为 %s 和 %s。",
          fmt_pct(ta$DTU_gene_nonDE_fraction), ta$b_DTU_and_nonDE, ta$n_strict_DTU_genes,
          fmt_pct(tn$DTU_gene_nonDE_fraction), tn$b_DTU_and_nonDE, tn$n_strict_DTU_genes,
          fmt_pct(ta$all_common_nonDE_fraction), fmt_pct(tn$all_common_nonDE_fraction)),
  "",
  sprintf("**D，** 对全部 common eligible genes 计算短读长基因 |log2FC| 与长读长 within-gene isoform-composition TVD 的 Spearman 相关。T vs A: rho=%.3f, 95%% bootstrap CI %.3f-%.3f, P %s, n=%s；T vs N: rho=%.3f, 95%% CI %.3f-%.3f, P %s, n=%s。CI 以 gene 为单位重抽样4000次。Strict DTU genes 仅作视觉突出，未单独用于主相关检验。",
          ta_tvd$Spearman_rho, ta_tvd$bootstrap_CI_low, ta_tvd$bootstrap_CI_high, fmt_p(ta_tvd$two_sided_P), comma(ta_tvd$n_genes),
          tn_tvd$Spearman_rho, tn_tvd$bootstrap_CI_low, tn_tvd$bootstrap_CI_high, fmt_p(tn_tvd$two_sided_P), comma(tn_tvd$n_genes)),
  "",
  "N，normal tissue；A，adjacent tissue；T，tumor tissue。Source data are provided with the figure."
)
write_lines2(legend_zh, file.path(dir_reports, "Figure3_v2_legend_zh.md"))

results_zh <- c(
  "# Figure 3 v2 Results 中文稿",
  "",
  "## Matched short-read gene expression contextualizes long-read isoform remodeling",
  "",
  sprintf("将长读长 DTU tested universe、长读长 tx2gene 映射和 matched short-read DESeq2 结果统一到稳定 Ensembl gene ID 后，T vs A 和 T vs N 均获得 %s 个 common testable multi-isoform genes。两个 universe 分别包含 %d 和 %d 个 short-read DEGs，背景 DEG 率为 %s 和 %s。",
          comma(ta$n_common_testable_genes), ta$n_shortread_DE, tn$n_shortread_DE,
          fmt_pct(ta$all_common_DE_fraction), fmt_pct(tn$all_common_DE_fraction)),
  "",
  sprintf("在 strict DTU genes 中，T vs A 有 %d/%d (%s) 同时为 short-read DEG，T vs N 为 %d/%d (%s)。与 non-DTU genes 比较，strict DTU genes 显著富集于 short-read DEGs（T vs A: OR=%.2f, 95%% CI %.2f-%.2f, P %s；T vs N: OR=%.2f, 95%% CI %.2f-%.2f, P %s）。",
          ta$a_DTU_and_DE, ta$n_strict_DTU_genes, fmt_pct(ta$DTU_gene_DE_fraction),
          tn$a_DTU_and_DE, tn$n_strict_DTU_genes, fmt_pct(tn$DTU_gene_DE_fraction),
          ta$odds_ratio, ta$odds_ratio_CI_low, ta$odds_ratio_CI_high, fmt_p(ta$Fisher_two_sided_P),
          tn$odds_ratio, tn$odds_ratio_CI_low, tn$odds_ratio_CI_high, fmt_p(tn$Fisher_two_sided_P)),
  "",
  sprintf("尽管存在上述富集，多数 strict DTU genes 仍未达到预设 short-read DEG 标准：T vs A 为 %d/%d (%s)，T vs N 为 %d/%d (%s)。但 common universe 整体的 non-DE 比例更高，分别为 %s 和 %s。因此，%s/%s 不应被孤立解释为 DTU genes 在 DEGs 中 depleted。",
          ta$b_DTU_and_nonDE, ta$n_strict_DTU_genes, fmt_pct(ta$DTU_gene_nonDE_fraction),
          tn$b_DTU_and_nonDE, tn$n_strict_DTU_genes, fmt_pct(tn$DTU_gene_nonDE_fraction),
          fmt_pct(ta$all_common_nonDE_fraction), fmt_pct(tn$all_common_nonDE_fraction),
          fmt_pct(ta$DTU_gene_nonDE_fraction), fmt_pct(tn$DTU_gene_nonDE_fraction)),
  "",
  sprintf("在不基于显著性预筛选的全 common-gene 连续效应分析中，|gene log2FC| 与 TVD 呈弱正相关（T vs A: rho=%.3f, 95%% CI %.3f-%.3f, P %s；T vs N: rho=%.3f, 95%% CI %.3f-%.3f, P %s）。以 max absolute delta usage 作为敏感性指标时方向一致（rho=%.3f 和 %.3f）。这表明基因丰度变化幅度与 isoform-composition remodeling 存在统计关联，但效应较弱，两者不能互相替代。",
          ta_tvd$Spearman_rho, ta_tvd$bootstrap_CI_low, ta_tvd$bootstrap_CI_high, fmt_p(ta_tvd$two_sided_P),
          tn_tvd$Spearman_rho, tn_tvd$bootstrap_CI_low, tn_tvd$bootstrap_CI_high, fmt_p(tn_tvd$two_sided_P),
          ta_max$Spearman_rho, tn_max$Spearman_rho)
)
write_lines2(results_zh, file.path(dir_reports, "Figure3_v2_results_zh.md"))

methods_zh <- c(
  "# Figure 3 v2 Methods 中文稿",
  "",
  "## Common testable gene universe and matched gene-isoform analysis",
  "",
  "分析限定于 primary72 队列，即排除 CRC-20 后24例患者的72份完整配对 N/A/T 组织。长、短读长数据来自同一 RNA 提取物分装。本阶段仅复用已冻结的 matched short-read DESeq2 结果和 primary satuRn + stageR DTU 结果，未重跑上游定量或统计模型。",
  "",
  "T vs A 和 T vs N 的 common testable multi-isoform gene universe 分别定义。纳入基因必须：存在于 matched short-read 结果；通过去版本稳定 Ensembl gene ID 与长读长 tx2gene 唯一映射；属于正式 DTU tested universe；至少含2个正式 usage-analysis isoforms；具有有效 short-read log2FC 和 long-read usage change。Gene symbol 仅用于显示，不用于匹配。",
  "",
  "Matched short-read DEG 固定定义为 DESeq2 adjusted P < 0.05 且 |log2FC| > 1。Strict DTU 固定定义为 transcript-level adjusted FDR < 0.10 且 |delta usage| >= 0.05。在 common universe 内按 DTU/non-DTU 和 DEG/non-DEG 构建2x2列联表，使用双侧 Fisher exact test，报告 odds ratio、95% CI 和 P 值。",
  "",
  "Within-gene isoform-composition change 使用 total variation distance (TVD) 度量：TVD = 0.5 x sum_i |usage_group2,i - usage_group1,i|。T vs A 中 group2/group1 为 Tumor/Adjacent，T vs N 为 Tumor/Normal。求和遍历该基因所有正式 eligible isoforms，而不仅是 strict DTU isoforms。敏感性指标为基因内 max_i |delta usage_i|。",
  "",
  "分别在 T vs A 和 T vs N 的全部 common genes 上计算 |matched short-read gene log2FC| 与 TVD 的 Spearman 相关，使用双侧 P 值。95% CI 通过以 gene 为重抽样单位的4000次 bootstrap 获得，随机种子为20260828。Max absolute delta usage 使用相同方法进行敏感性分析。",
  "",
  "阈值敏感性分析穷举组合 DTU FDR 0.05/0.10、|delta usage| 0.03/0.05/0.10、short-read padj 0.01/0.05/0.10 和 |log2FC| 0.5/1.0。每个组合均使用同一 common universe，并记录 DTU genes、overlap 和 DTU-gene non-DE fraction。Primary definition 始终保持为 DTU FDR < 0.10 + |delta usage| >= 0.05 与 DEG padj < 0.05 + |log2FC| > 1。"
)
write_lines2(methods_zh, file.path(dir_reports, "Figure3_v2_methods_zh.md"))

ta_thr <- threshold_range %>% filter(contrast == "T_vs_A")
tn_thr <- threshold_range %>% filter(contrast == "T_vs_N")
ta_primary_thr <- threshold_primary_dtu_range %>% filter(contrast == "T_vs_A")
tn_primary_thr <- threshold_primary_dtu_range %>% filter(contrast == "T_vs_N")

redesign_report <- c(
  "# FIG3 STATISTICAL REDESIGN REPORT",
  "",
  "## 执行边界",
  "",
  "- 本轮只重算 Figure 3 所需派生统计并生成新 `_v2` 图件。",
  "- 未重跑 FASTQ、Kallisto、FLAIR、SQANTI3、DTE 或 primary DTU。",
  "- Figure 1/2/4-8 由 SHA256 guard 保护。",
  "- Matched short-read 是同一队列、同一 RNA 提取物分装的 complementary analysis，不是独立外部验证。",
  "",
  "## 1. Common testable gene universe",
  "",
  sprintf("- T vs A: %s genes。", comma(ta$n_common_testable_genes)),
  sprintf("- T vs N: %s genes。", comma(tn$n_common_testable_genes)),
  "- 两个 universe 均来自正式 DTU tested multi-isoform genes，所有基因均完成稳定 Ensembl ID 唯一映射。",
  "",
  "## 2. Common-universe short-read DEG baseline rate",
  "",
  sprintf("- T vs A: %d/%d = %s。", ta$n_shortread_DE, ta$n_common_testable_genes, fmt_pct(ta$all_common_DE_fraction)),
  sprintf("- T vs N: %d/%d = %s。", tn$n_shortread_DE, tn$n_common_testable_genes, fmt_pct(tn$all_common_DE_fraction)),
  "",
  "## 3. Strict DTU genes 中的 DEG rate",
  "",
  sprintf("- T vs A: %d/%d = %s。", ta$a_DTU_and_DE, ta$n_strict_DTU_genes, fmt_pct(ta$DTU_gene_DE_fraction)),
  sprintf("- T vs N: %d/%d = %s。", tn$a_DTU_and_DE, tn$n_strict_DTU_genes, fmt_pct(tn$DTU_gene_DE_fraction)),
  "",
  "## 4. Fisher exact test",
  "",
  sprintf("- T vs A: OR=%.4f, 95%% CI %.4f-%.4f, two-sided P=%s。", ta$odds_ratio, ta$odds_ratio_CI_low, ta$odds_ratio_CI_high, fmt_p(ta$Fisher_two_sided_P)),
  sprintf("- T vs N: OR=%.4f, 95%% CI %.4f-%.4f, two-sided P=%s。", tn$odds_ratio, tn$odds_ratio_CI_low, tn$odds_ratio_CI_high, fmt_p(tn$Fisher_two_sided_P)),
  "",
  "## 5. Enrichment/depletion 方向",
  "",
  sprintf("- T vs A: %s。", ta$direction),
  sprintf("- T vs N: %s。", tn$direction),
  "- Strict DTU genes 在 short-read DEGs 中是 enrichment，不是 depletion。",
  "",
  "## 6. |gene log2FC| 与 TVD 的连续关联",
  "",
  sprintf("- T vs A: rho=%.4f, 95%% bootstrap CI %.4f-%.4f, P=%s, n=%d。", ta_tvd$Spearman_rho, ta_tvd$bootstrap_CI_low, ta_tvd$bootstrap_CI_high, fmt_p(ta_tvd$two_sided_P), ta_tvd$n_genes),
  sprintf("- T vs N: rho=%.4f, 95%% bootstrap CI %.4f-%.4f, P=%s, n=%d。", tn_tvd$Spearman_rho, tn_tvd$bootstrap_CI_low, tn_tvd$bootstrap_CI_high, fmt_p(tn_tvd$two_sided_P), tn_tvd$n_genes),
  "- 两个 contrast 均为统计显著但效应较弱的正相关。",
  "",
  "## 7. Max absolute delta usage sensitivity",
  "",
  sprintf("- T vs A: rho=%.4f, 95%% CI %.4f-%.4f, P=%s。", ta_max$Spearman_rho, ta_max$bootstrap_CI_low, ta_max$bootstrap_CI_high, fmt_p(ta_max$two_sided_P)),
  sprintf("- T vs N: rho=%.4f, 95%% CI %.4f-%.4f, P=%s。", tn_max$Spearman_rho, tn_max$bootstrap_CI_low, tn_max$bootstrap_CI_high, fmt_p(tn_max$two_sided_P)),
  "- 两个敏感性结果与 TVD 主分析方向一致。",
  "",
  "## 8. 79.8%/78.4% 的阈值稳定性",
  "",
  sprintf("- 所有72个预设组合中，T vs A fraction_nonDE 范围为 %s-%s，T vs N 为 %s-%s。",
          fmt_pct(ta_thr$min_fraction), fmt_pct(ta_thr$max_fraction), fmt_pct(tn_thr$min_fraction), fmt_pct(tn_thr$max_fraction)),
  sprintf("- 固定 primary DTU 定义、仅变化 short-read DEG 阈值时，T vs A 范围为 %s-%s，T vs N 为 %s-%s。",
          fmt_pct(ta_primary_thr$min_fraction), fmt_pct(ta_primary_thr$max_fraction), fmt_pct(tn_primary_thr$min_fraction), fmt_pct(tn_primary_thr$max_fraction)),
  "- 因此79.8%/78.4%仅对冻结 primary thresholds 成立，在全部合理阈值网格下并非数值稳定常数。",
  "",
  "## 9. 推荐主文措辞",
  "",
  "> Within the common testable multi-isoform gene universe, most strict DTU genes did not meet the stringent matched short-read DEG criterion; however, DTU genes were enriched among short-read DEGs relative to non-DTU genes, and gene-level abundance-change magnitude showed only a weak positive association with isoform-composition remodeling.",
  "",
  "中文：在统一可检验的多转录本基因集中，多数 strict DTU genes 未达到严格的 matched short-read DEG 标准；但与 non-DTU genes 相比，DTU genes 在 short-read DEGs 中显著富集，且基因丰度变化幅度与 isoform-composition remodeling 仅呈弱正相关。",
  "",
  "## 10. 禁止的过度表述",
  "",
  "- 不得写 gene-level expression 完全看不到 isoform remodeling。",
  "- 不得写 DTU genes 在 DEGs 中 depleted；Fisher 结果方向为 enrichment。",
  "- 不得将 79.8%/78.4% 写成与阈值无关的稳定生物学常数。",
  "- 不得将 matched short-read 写成 independent external validation。",
  "- 不得写 short reads validate complete full-length isoforms。",
  "- 不得宣称测序平台优越性、因果机制或功能验证。",
  "",
  "## 新输出",
  "",
  paste0("- Main Figure 3 v2: `", main_paths$path[main_paths$format == "pdf"], "`"),
  paste0("- Selected-event supplementary: `", selected_paths$path[selected_paths$format == "pdf"], "`"),
  paste0("- Threshold sensitivity: `", threshold_paths$path[threshold_paths$format == "pdf"], "`"),
  paste0("- Previous-v2 archive: `", ifelse(is.na(archive_path), "none_required", archive_path), "`")
)
write_lines2(redesign_report, file.path(dir_reports, "FIG3_STATISTICAL_REDESIGN_REPORT.md"))

record_progress("QA", "Running source-data, field-count, boundary and unrelated-figure guards")
source_index_check <- source_index %>%
  mutate(
    exists_nonempty = file.exists(source_data_path) & file.info(source_data_path)$size > 0,
    observed_rows = vapply(source_data_path, count_tsv_rows, integer(1)),
    row_count_match = observed_rows == n_rows
  )
stop_if(any(!source_index_check$exists_nonempty | !source_index_check$row_count_match),
        "Figure 3 v2 source-data index completeness QA failed.")

tsv_paths <- c(
  file.path(dir_results, "Fig3_v2_input_manifest.tsv"),
  file.path(dir_results, "Fig3_common_testable_gene_universe.tsv"),
  file.path(dir_results, "Fig3_common_testable_gene_universe_summary.tsv"),
  file.path(dir_results, "Fig3_common_universe_mapping_audit.tsv"),
  file.path(dir_results, "Fig3_DTU_DE_contingency_fisher.tsv"),
  file.path(dir_results, "Fig3_continuous_gene_isoform_coupling.tsv"),
  file.path(dir_results, "Fig3_continuous_gene_isoform_coupling_summary.tsv"),
  file.path(dir_results, "Fig3_threshold_sensitivity.tsv"),
  source_index_path,
  unname(source_paths)
)
field_qa <- tsv_field_audit(tsv_paths)
stop_if(any(!field_qa$field_count_pass), "One or more Figure 3 v2 TSV files have inconsistent field counts.")
field_qa_path <- write_tsv2(field_qa, file.path(dir_results, "Fig3_v2_TSV_field_QA.tsv"))

figure_paths <- c(main_paths$path, selected_paths$path, threshold_paths$path)
stop_if(any(!file.exists(figure_paths) | file.info(figure_paths)$size <= 0),
        "One or more Figure 3 v2 figure files are missing or empty.")

svg_paths <- figure_paths[grepl("[.]svg$", figure_paths)]
forbidden_terms <- c(
  "candidate prioritization", "assay feasibility", "RT-ready", "validation readiness",
  "clinical association", "prognostic model", "biomarker", "independent external validation",
  "validated full-length isoform", "confirmed mechanism"
)
boundary_hits <- bind_rows(lapply(svg_paths, function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  tibble(
    path = path,
    term = forbidden_terms,
    hit = vapply(forbidden_terms, function(term) grepl(tolower(term), tolower(text), fixed = TRUE), logical(1))
  )
}))
stop_if(any(boundary_hits$hit), "A Figure 3 v2 SVG contains a prohibited manuscript-boundary term.")
write_tsv2(boundary_hits, file.path(dir_results, "Fig3_v2_boundary_grep.tsv"))

guard_after <- file_snapshot(guard_paths)
guard_check <- guard_before %>%
  rename(size_bytes_before = size_bytes, mtime_before = mtime, sha256_before = sha256) %>%
  left_join(
    guard_after %>% rename(size_bytes_after = size_bytes, mtime_after = mtime, sha256_after = sha256),
    by = "path"
  ) %>%
  mutate(
    unchanged = !is.na(sha256_after) & size_bytes_before == size_bytes_after & sha256_before == sha256_after
  )
stop_if(any(!guard_check$unchanged), "Figure 1/2/4-8 scope guard failed.")
guard_path <- write_tsv2(guard_check, file.path(dir_results, "Fig3_v2_unrelated_figure_guard.tsv"))

output_manifest <- tibble(
  output_path = c(
    figure_paths, tsv_paths, field_qa_path, guard_path,
    file.path(dir_reports, c(
      "FIG3_STATISTICAL_REDESIGN_REPORT.md", "Figure3_v2_legend_zh.md",
      "Figure3_v2_results_zh.md", "Figure3_v2_methods_zh.md"
    ))
  )
) %>%
  distinct(output_path) %>%
  mutate(
    exists = file.exists(output_path),
    size_bytes = file.info(output_path)$size,
    sha256 = vapply(output_path, sha256_file, character(1))
  )
stop_if(any(!output_manifest$exists | output_manifest$size_bytes <= 0), "Output manifest contains a missing/empty file.")
write_tsv2(output_manifest, file.path(dir_results, "Fig3_v2_output_manifest.tsv"))

qa_report <- c(
  "# Figure 3 v2 programmatic QA",
  "",
  "- Input SHA256 gates: PASS.",
  "- Primary72 cohort gate (24 patients / 72 N/A/T tissues): PASS.",
  "- Strict DTU unit/gene gates (0/155/118; 99/74): PASS.",
  "- Common-universe gate (1,503 genes per tumor contrast): PASS.",
  "- Stable Ensembl mapping uses no gene-symbol matching; ambiguity exclusions are audited.",
  "- Within-gene usage sums close to one and TVD lies in [0,1]: PASS.",
  "- Fisher tables and known 20/16 overlap gates: PASS.",
  "- Spearman bootstrap: 4,000 gene-level resamples per contrast and metric.",
  "- Threshold grid: 72/72 combinations complete; primary cells match frozen results.",
  "- TSV field counts and source-data row counts: PASS.",
  "- Figure PDF/SVG/PNG existence and non-empty checks: PASS.",
  "- Figure SVG boundary grep: PASS.",
  sprintf("- Figure 1/2/4-8 guard: %d/%d unchanged.", sum(guard_check$unchanged), nrow(guard_check)),
  "- Final local visual QA remains required after PNG download."
)
write_lines2(qa_report, file.path(dir_reports, "Figure3_v2_QA_report_zh.md"))

record_progress("COMPLETE", "Figure 3 common-universe statistical redesign completed")
cat("Figure 3 v2 statistical redesign completed.\n")
cat("Main figure:", main_paths$path[main_paths$format == "pdf"], "\n")
cat("Report:", file.path(dir_reports, "FIG3_STATISTICAL_REDESIGN_REPORT.md"), "\n")
cat("No upstream quantification, DTE or primary DTU model was rerun.\n")
