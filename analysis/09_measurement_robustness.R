#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(20260829)

required_packages <- c(
  "DESeq2", "dplyr", "tidyr", "readr", "tibble", "ggplot2",
  "patchwork", "scales", "ragg", "systemfonts"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

ROOT <- Sys.getenv(
  "V2_BASE",
  Sys.getenv("CRC_V2_ROOT")
)
STUDY_ROOT <- dirname(ROOT)
FRAMEWORK_ROOT <- file.path(STUDY_ROOT, "16_paired_NAT_multi_platform_framework_20260827")
PHASE03_ROOT <- file.path(STUDY_ROOT, "03_DTU_paired", "tables")

DIR_RESULTS <- file.path(ROOT, "results")
DIR_FIGURES <- file.path(ROOT, "figures")
DIR_REPORTS <- file.path(ROOT, "reports")
DIR_LOGS <- file.path(ROOT, "logs")
DIR_CHECKPOINT <- file.path(ROOT, "checkpoints", "Fig3_measurement_robustness")
DIR_ARCHIVE <- file.path(ROOT, "archive")
for (d in c(DIR_RESULTS, DIR_FIGURES, DIR_REPORTS, DIR_LOGS, DIR_CHECKPOINT, DIR_ARCHIVE)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PROGRESS_FILE <- file.path(DIR_LOGS, "29_Fig3_measurement_robustness_progress.txt")
record_progress <- function(stage, detail) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), stage, detail, sep = "\t")
  cat(line, "\n", file = PROGRESS_FILE, append = TRUE)
  message("[", stage, "] ", detail)
}
stop_if <- function(condition, message_text) {
  if (isTRUE(condition)) {
    record_progress("FAIL", message_text)
    stop(message_text, call. = FALSE)
  }
}
assert_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  stop_if(length(missing) > 0L, paste0(label, " missing columns: ", paste(missing, collapse = ", ")))
}
clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x) | x == "." | toupper(x) == "NA"] <- NA_character_
  x
}
strip_ensembl_version <- function(x) sub("[.][0-9]+$", "", clean_chr(x))
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
count_tsv_rows <- function(path) max(0L, length(readLines(path, warn = FALSE)) - 1L)
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
  target <- file.path(DIR_ARCHIVE, paste0("Fig3_robustness_candidate_before_", stamp))
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(existing, file.path(target, basename(existing)), overwrite = FALSE)
  stop_if(!all(ok), "Could not archive all pre-existing robustness outputs.")
  target
}
atomic_save_rds <- function(object, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(object, tmp, compress = FALSE)
  stop_if(!file.rename(tmp, path), paste0("Could not atomically save checkpoint: ", path))
  invisible(path)
}
safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = method))
}
safe_cor_test_p <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) return(NA_real_)
  suppressWarnings(cor.test(x[keep], y[keep], method = method, exact = FALSE)$p.value)
}
bootstrap_cor <- function(x, y, B = 4000L, seed = 20260829L, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  stop_if(n < 4L, "Insufficient observations for bootstrap correlation.")
  set.seed(seed)
  vals <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    safe_cor(x[idx], y[idx], method = method)
  })
  vals <- vals[is.finite(vals)]
  stop_if(length(vals) < floor(0.95 * B), "Too few valid bootstrap correlations.")
  c(
    estimate = safe_cor(x, y, method = method),
    p = safe_cor_test_p(x, y, method = method),
    ci_low = unname(quantile(vals, 0.025)),
    ci_high = unname(quantile(vals, 0.975)),
    bootstrap_valid = length(vals),
    n = n
  )
}
normal_score <- function(x) {
  n <- length(x)
  qnorm((rank(x, ties.method = "average") - 0.5) / n)
}
quartile_by_rank <- function(x, stable_id) {
  stop_if(length(x) != length(stable_id), "Quartile inputs differ in length.")
  ord <- order(x, stable_id, na.last = NA)
  out <- rep(NA_integer_, length(x))
  out[ord] <- ceiling(seq_along(ord) * 4 / length(ord))
  factor(paste0("Q", out), levels = paste0("Q", 1:4))
}
multiplicity_group <- function(x) {
  factor(
    ifelse(x == 2, "2", ifelse(x == 3, "3", ifelse(x <= 5, "4-5", ">=6"))),
    levels = c("2", "3", "4-5", ">=6")
  )
}
summarise_distribution <- function(x) {
  x <- x[is.finite(x)]
  tibble(
    n_replicates = length(x),
    median = median(x),
    Q1 = unname(quantile(x, 0.25)),
    Q3 = unname(quantile(x, 0.75)),
    percentile_2.5 = unname(quantile(x, 0.025)),
    percentile_97.5 = unname(quantile(x, 0.975)),
    mean = mean(x),
    SD = sd(x)
  )
}
tsv_field_audit <- function(paths) {
  bind_rows(lapply(paths, function(path) {
    lines <- readLines(path, warn = FALSE)
    stop_if(length(lines) == 0L, paste0("Empty TSV: ", path))
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

EXPECTED_HASHES <- c(
  common_universe = "b5e7b1a109bb9783b54d5a61fae581e6a8284a752abf59817e5a2cac8af044fd",
  frozen_coupling = "3733f7fc591a39d8f160aaec286a7f4ef9c2e09e40c4efc048835348473f9825",
  frozen_coupling_summary = "9b70adf27dfde449c5507c80859032d3656ea55c37ad0bf670464c20f62545f5",
  matched_shortread = "2aca124a299017ff1b52aeb6eb3b22b8479ed179a37668ddbdac14c216974a68",
  contract = "4663ad2f344c5ecdd9ebc6c9faca879ad9a8ed074783df06e49494d474098f51",
  dtu_T_vs_A = "60074de59d6ed821a4a8be472405179beabe3aa80dc3fa7a9eaa15eff1f0e0f6",
  dtu_T_vs_N = "a3d28b4631b04b2fd0b9eef311a1af4a42d3b537c97f82fc0bb06ccb0a91d6c0"
)
INPUT_PATHS <- c(
  common_universe = file.path(DIR_RESULTS, "Fig3_common_testable_gene_universe.tsv"),
  frozen_coupling = file.path(DIR_RESULTS, "Fig3_continuous_gene_isoform_coupling.tsv"),
  frozen_coupling_summary = file.path(DIR_RESULTS, "Fig3_continuous_gene_isoform_coupling_summary.tsv"),
  matched_shortread = file.path(ROOT, "tables", "matched_shortread_primary72_DE_results.tsv"),
  contract = file.path(FRAMEWORK_ROOT, "validated_data_objects.rds"),
  dtu_T_vs_A = file.path(PHASE03_ROOT, "primary72_DTU_isoform_T_vs_A.tsv"),
  dtu_T_vs_N = file.path(PHASE03_ROOT, "primary72_DTU_isoform_T_vs_N.tsv")
)

record_progress("INPUT_AUDIT", "Checking frozen Figure 3 sources and hashes")
stop_if(any(!file.exists(INPUT_PATHS)), paste0("Missing required input(s): ", paste(INPUT_PATHS[!file.exists(INPUT_PATHS)], collapse = "; ")))
observed_hashes <- vapply(INPUT_PATHS, sha256_file, character(1))
stop_if(!identical(unname(observed_hashes[names(EXPECTED_HASHES)]), unname(EXPECTED_HASHES)), "One or more frozen input hashes changed.")
input_manifest <- tibble(
  input_role = names(INPUT_PATHS),
  path = unname(INPUT_PATHS),
  size_bytes = file.info(unname(INPUT_PATHS))$size,
  mtime = format(file.info(unname(INPUT_PATHS))$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  sha256 = unname(observed_hashes),
  access_mode = "read-only",
  gate_status = "PASS"
)
write_tsv2(input_manifest, file.path(DIR_RESULTS, "Fig3_robustness_input_manifest.tsv"))

guard_roots <- c(file.path(ROOT, "figures", "main"), file.path(ROOT, "source_data"), file.path(ROOT, "tables"))
guard_candidates <- unlist(lapply(guard_roots, function(d) list.files(d, recursive = TRUE, full.names = TRUE)))
guard_candidates <- guard_candidates[file.exists(guard_candidates) & !dir.exists(guard_candidates)]
guard_pattern <- "^(Fig|Figure)(1|2|4|5|6|7|8)([^0-9]|$)|^Figure3_matched_shortread_longread_v2[.]"
guard_paths <- guard_candidates[grepl(guard_pattern, basename(guard_candidates))]
guard_before <- file_snapshot(guard_paths)

record_progress("LOAD", "Loading frozen universe, eligible isoforms and primary72 count matrices")
common <- read_tsv(INPUT_PATHS[["common_universe"]], show_col_types = FALSE)
assert_columns(common, c(
  "contrast", "gene_id", "gene_symbol", "stable_gene_id", "shortread_log2FC",
  "shortread_padj", "shortread_baseMean", "shortread_DE_status", "longread_DTU_status",
  "n_tested_isoforms", "TVD_usage_change", "max_abs_delta_usage"
), "Frozen common universe")
contrast_levels <- c("T_vs_A", "T_vs_N")
common <- common %>% filter(contrast %in% contrast_levels)
common_counts <- common %>% count(contrast)
stop_if(!identical(common_counts$contrast, contrast_levels) || !identical(common_counts$n, c(1503L, 1503L)), "Common universe gate is not 1,503 genes per tumor contrast.")
common_ta <- common %>% filter(contrast == "T_vs_A") %>% arrange(stable_gene_id)
common_tn <- common %>% filter(contrast == "T_vs_N") %>% arrange(stable_gene_id)
stop_if(!identical(common_ta$stable_gene_id, common_tn$stable_gene_id), "T vs A and T vs N common-gene sets differ.")
COMMON_IDS <- common_ta$stable_gene_id

rho_gate_ta <- safe_cor(abs(common_ta$shortread_log2FC), common_ta$TVD_usage_change)
rho_gate_tn <- safe_cor(abs(common_tn$shortread_log2FC), common_tn$TVD_usage_change)
stop_if(abs(rho_gate_ta - 0.13303076963618166) > 1e-12 || abs(rho_gate_tn - 0.10173880855452026) > 1e-12,
        "Frozen observed Spearman rho could not be reproduced.")

contract <- readRDS(INPUT_PATHS[["contract"]])
stop_if(!all(c("lr_counts", "sr_gene_counts", "tx2gene", "meta") %in% names(contract)), "Validated contract is incomplete.")
lr_counts <- contract$lr_counts
sr_counts <- contract$sr_gene_counts
tx2gene <- as.data.frame(contract$tx2gene)
meta <- as.data.frame(contract$meta)
assert_columns(meta, c("sample_id", "patient_id", "tissue"), "Metadata")
assert_columns(tx2gene, c("isoform_id", "gene_id", "gene_symbol"), "tx2gene")
stop_if(!identical(colnames(lr_counts), meta$sample_id) || !identical(colnames(sr_counts), meta$sample_id), "Count matrix sample order differs from metadata.")

primary_keep <- meta$patient_id != "CRC_20"
primary_meta <- meta[primary_keep, , drop = FALSE]
primary_meta$tissue <- factor(as.character(primary_meta$tissue), levels = c("N", "A", "T"))
primary_meta$patient_id <- factor(primary_meta$patient_id)
PRIMARY_SAMPLES <- primary_meta$sample_id
rownames(primary_meta) <- PRIMARY_SAMPLES
stop_if(nrow(primary_meta) != 72L || nlevels(primary_meta$patient_id) != 24L, "primary72 is not 24 patients / 72 samples.")
stop_if(!identical(as.integer(table(primary_meta$tissue)), c(24L, 24L, 24L)), "primary72 N/A/T counts are not 24/24/24.")
pairing <- table(primary_meta$patient_id, primary_meta$tissue)
stop_if(!all(pairing == 1L), "primary72 patient pairing is incomplete.")
stop_if(!identical(PRIMARY_SAMPLES, colnames(lr_counts)[primary_keep]) || !identical(PRIMARY_SAMPLES, colnames(sr_counts)[primary_keep]),
        "Long-read and short-read primary72 sample order mismatch.")
sample_audit <- primary_meta %>%
  transmute(
    sample_id, patient_id = as.character(patient_id), tissue = as.character(tissue),
    primary72_included = TRUE, longread_present = sample_id %in% colnames(lr_counts),
    shortread_present = sample_id %in% colnames(sr_counts), paired_triplet = TRUE
  )
write_tsv2(sample_audit, file.path(DIR_RESULTS, "Fig3_robustness_sample_audit.tsv"))

dtu_ta <- read_tsv(INPUT_PATHS[["dtu_T_vs_A"]], show_col_types = FALSE)
dtu_tn <- read_tsv(INPUT_PATHS[["dtu_T_vs_N"]], show_col_types = FALSE)
assert_columns(dtu_ta, c("gene_id", "gene_name", "isoform_id", "mean_usage_N", "mean_usage_A", "mean_usage_T"), "T vs A DTU table")
stop_if(nrow(dtu_ta) != 5477L || nrow(dtu_tn) != 5477L, "Eligible isoform row count changed from 5,477.")
stop_if(!setequal(dtu_ta$isoform_id, dtu_tn$isoform_id), "Eligible isoform memberships differ across tumor contrasts.")
dtu_ta <- dtu_ta %>%
  mutate(stable_gene_id = strip_ensembl_version(gene_id), isoform_id = clean_chr(isoform_id)) %>%
  arrange(stable_gene_id, isoform_id)
stop_if(anyDuplicated(dtu_ta$isoform_id) > 0L, "Eligible isoform IDs are duplicated.")
stop_if(!setequal(unique(dtu_ta$stable_gene_id), COMMON_IDS), "Eligible DTU genes differ from the frozen common universe.")
eligible_index <- match(dtu_ta$isoform_id, rownames(lr_counts))
stop_if(anyNA(eligible_index), "At least one eligible isoform is absent from lr_counts.")
eligible_counts <- lr_counts[eligible_index, PRIMARY_SAMPLES, drop = FALSE]
stop_if(any(abs(eligible_counts - round(eligible_counts)) > 1e-8), "Long-read eligible counts are not integer.")
storage.mode(eligible_counts) <- "integer"
ELIGIBLE_GENE <- dtu_ta$stable_gene_id
ELIGIBLE_GENE_LEVELS <- unique(ELIGIBLE_GENE)
stop_if(!setequal(ELIGIBLE_GENE_LEVELS, COMMON_IDS), "Eligible gene order audit failed.")

tx_clean <- tx2gene %>%
  transmute(
    isoform_id = clean_chr(isoform_id),
    stable_gene_id = strip_ensembl_version(gene_id),
    gene_symbol_tx = clean_chr(gene_symbol)
  )
stop_if(anyDuplicated(tx_clean$isoform_id) > 0L, "tx2gene isoform mapping is duplicated.")
lr_map_index <- match(rownames(lr_counts), tx_clean$isoform_id)
stop_if(anyNA(lr_map_index), "lr_counts rows are not fully covered by tx2gene.")
lr_stable_gene <- tx_clean$stable_gene_id[lr_map_index]
all_common_lr_keep <- lr_stable_gene %in% COMMON_IDS
all_common_lr_counts <- lr_counts[all_common_lr_keep, PRIMARY_SAMPLES, drop = FALSE]
all_common_lr_gene <- lr_stable_gene[all_common_lr_keep]
all_lr_gene_totals <- rowsum(all_common_lr_counts, group = all_common_lr_gene, reorder = FALSE)
all_lr_gene_totals <- all_lr_gene_totals[match(COMMON_IDS, rownames(all_lr_gene_totals)), , drop = FALSE]
stop_if(anyNA(rownames(all_lr_gene_totals)), "Could not order all long-read gene totals.")
lr_library_totals <- colSums(lr_counts[, PRIMARY_SAMPLES, drop = FALSE])
stop_if(any(lr_library_totals <= 0), "A primary72 long-read library has zero total count.")
all_lr_gene_cpm <- sweep(all_lr_gene_totals, 2L, lr_library_totals, "/") * 1e6

stable_sr_ids <- strip_ensembl_version(rownames(sr_counts))
sr_duplicated <- duplicated(stable_sr_ids) | duplicated(stable_sr_ids, fromLast = TRUE)
sr_index <- match(COMMON_IDS, stable_sr_ids)
stop_if(anyNA(sr_index) || any(sr_duplicated[sr_index]), "Common genes do not uniquely map to short-read count rows.")
common_sr_counts <- sr_counts[sr_index, PRIMARY_SAMPLES, drop = FALSE]
rownames(common_sr_counts) <- COMMON_IDS
stop_if(!any(abs(common_sr_counts - round(common_sr_counts)) > 1e-8), "Short-read estimated counts unexpectedly contain no fractional values.")
stop_if(any(colSums(common_sr_counts) <= 0), "A common-universe short-read pseudo-library has zero depth.")

compute_usage_tvd <- function(count_matrix) {
  gene_totals <- rowsum(count_matrix, group = ELIGIBLE_GENE, reorder = FALSE)
  denominator <- gene_totals[match(ELIGIBLE_GENE, rownames(gene_totals)), , drop = FALSE]
  usage <- count_matrix / pmax(denominator, 1)
  usage[denominator == 0] <- NA_real_
  group_mean <- lapply(c("N", "A", "T"), function(tissue_name) {
    rowMeans(usage[, primary_meta$tissue == tissue_name, drop = FALSE], na.rm = TRUE)
  })
  names(group_mean) <- c("N", "A", "T")
  stop_if(any(!is.finite(unlist(group_mean))), "A pseudo-half produced a non-finite tissue mean usage.")
  make_tvd <- function(delta) {
    x <- rowsum(matrix(abs(delta), ncol = 1L), group = ELIGIBLE_GENE, reorder = FALSE)[, 1] * 0.5
    x[match(COMMON_IDS, names(x))]
  }
  list(
    T_vs_A = make_tvd(group_mean$T - group_mean$A),
    T_vs_N = make_tvd(group_mean$T - group_mean$N)
  )
}

record_progress("COVARIATES", "Reconstructing frozen TVD and defining measurement covariates")
reconstructed_tvd <- compute_usage_tvd(eligible_counts)
stop_if(max(abs(reconstructed_tvd$T_vs_A - common_ta$TVD_usage_change)) > 1e-12 ||
          max(abs(reconstructed_tvd$T_vs_N - common_tn$TVD_usage_change)) > 1e-12,
        "TVD could not be reconstructed from frozen eligible isoform counts.")

eligible_positive_samples <- rowSums(eligible_counts > 0)
EFFECTIVE_SUPPORT_SAMPLES <- 6L
effective_isoform_counts <- rowsum(
  matrix(as.integer(eligible_positive_samples >= EFFECTIVE_SUPPORT_SAMPLES), ncol = 1L),
  group = ELIGIBLE_GENE, reorder = FALSE
)[, 1]
effective_isoform_counts <- effective_isoform_counts[match(COMMON_IDS, names(effective_isoform_counts))]
zero_fraction_gene <- rowsum(
  matrix(rowMeans(eligible_counts == 0), ncol = 1L),
  group = ELIGIBLE_GENE, reorder = FALSE
)[, 1]
n_eligible_by_gene <- as.integer(table(factor(ELIGIBLE_GENE, levels = COMMON_IDS)))
zero_fraction_gene <- zero_fraction_gene[match(COMMON_IDS, names(zero_fraction_gene))] / n_eligible_by_gene

top_usage <- function(group_a, group_b) {
  top_a <- dtu_ta %>% group_by(stable_gene_id) %>% summarise(top = max(.data[[group_a]]), .groups = "drop")
  top_b <- dtu_ta %>% group_by(stable_gene_id) %>% summarise(top = max(.data[[group_b]]), .groups = "drop")
  (top_a$top[match(COMMON_IDS, top_a$stable_gene_id)] + top_b$top[match(COMMON_IDS, top_b$stable_gene_id)]) / 2
}

gene_base <- tibble(
  stable_gene_id = COMMON_IDS,
  lr_mean_gene_count = rowMeans(all_lr_gene_totals),
  lr_mean_gene_CPM = rowMeans(all_lr_gene_cpm),
  isoform_multiplicity = n_eligible_by_gene,
  zero_fraction = zero_fraction_gene,
  n_effective_isoforms = as.integer(effective_isoform_counts),
  effective_isoform_support_samples = EFFECTIVE_SUPPORT_SAMPLES
) %>%
  mutate(
    lr_abundance_quartile = quartile_by_rank(lr_mean_gene_count, stable_gene_id),
    multiplicity_group = multiplicity_group(isoform_multiplicity)
  )

covariates <- bind_rows(lapply(contrast_levels, function(z) {
  x <- if (z == "T_vs_A") common_ta else common_tn
  top <- if (z == "T_vs_A") top_usage("mean_usage_T", "mean_usage_A") else top_usage("mean_usage_T", "mean_usage_N")
  x %>%
    select(contrast, gene_id, gene_symbol, stable_gene_id, shortread_log2FC, shortread_padj,
           shortread_baseMean, shortread_DE_status, longread_DTU_status, n_tested_isoforms,
           TVD_usage_change, max_abs_delta_usage) %>%
    left_join(gene_base, by = "stable_gene_id") %>%
    mutate(
      mean_top_isoform_usage = top[match(stable_gene_id, COMMON_IDS)],
      sr_abundance_quartile = quartile_by_rank(shortread_baseMean, stable_gene_id),
      abs_shortread_log2FC = abs(shortread_log2FC),
      n_effective_isoforms_definition = paste0("eligible isoform count > 0 in at least ", EFFECTIVE_SUPPORT_SAMPLES, " of 72 primary samples")
    )
}))
stop_if(nrow(covariates) != 3006L || anyNA(covariates$lr_mean_gene_count) || anyNA(covariates$zero_fraction),
        "Measurement covariate table is incomplete.")
write_tsv2(covariates, file.path(DIR_RESULTS, "Fig3_measurement_covariates.tsv"))

lr_depth_stratum <- gene_base$lr_abundance_quartile
lr_mult_stratum <- gene_base$multiplicity_group

LR_REPLICATES <- 500L
SR_REPLICATES <- 200L
N_CORES <- max(1L, as.integer(Sys.getenv("ROBUSTNESS_CORES", "6")))

run_checkpointed <- function(kind, replicates, worker, checkpoint_dir, batch_size = 12L) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  path_for <- function(i) file.path(checkpoint_dir, sprintf("%s_rep_%04d.rds", kind, i))
  todo <- replicates[!file.exists(vapply(replicates, path_for, character(1)))]
  if (length(todo) > 0L) {
    batches <- split(todo, ceiling(seq_along(todo) / batch_size))
    for (batch in batches) {
      record_progress(toupper(kind), paste0("Running replicates ", min(batch), "-", max(batch), " using ", min(N_CORES, length(batch)), " workers"))
      status <- parallel::mclapply(batch, function(i) {
        tryCatch({
          result <- worker(i)
          atomic_save_rds(result, path_for(i))
          list(ok = TRUE, replicate = i, error = NA_character_)
        }, error = function(e) {
          list(ok = FALSE, replicate = i, error = conditionMessage(e))
        })
      }, mc.cores = min(N_CORES, length(batch)), mc.preschedule = FALSE)
      failed <- status[!vapply(status, `[[`, logical(1), "ok")]
      if (length(failed) > 0L) {
        msg <- paste(vapply(failed, function(x) paste0(x$replicate, ": ", x$error), character(1)), collapse = "; ")
        stop_if(TRUE, paste0(kind, " checkpoint batch failed: ", msg))
      }
    }
  } else {
    record_progress(toupper(kind), "All requested checkpoints already exist; resuming from disk")
  }
  paths <- vapply(replicates, path_for, character(1))
  stop_if(any(!file.exists(paths)), paste0("Missing ", kind, " checkpoint(s) after run."))
  lapply(paths, readRDS)
}

lr_worker <- function(rep_id) {
  seed <- 20260829L + rep_id
  set.seed(seed)
  half1 <- matrix(
    rbinom(length(eligible_counts), size = as.vector(eligible_counts), prob = 0.5),
    nrow = nrow(eligible_counts), ncol = ncol(eligible_counts),
    dimnames = dimnames(eligible_counts)
  )
  half2 <- eligible_counts - half1
  tvd1 <- compute_usage_tvd(half1)
  tvd2 <- compute_usage_tvd(half2)
  all_rows <- bind_rows(lapply(contrast_levels, function(z) {
    tibble(
      replicate = rep_id, replicate_seed = seed, contrast = z,
      n_genes = length(COMMON_IDS), reliability_rho = safe_cor(tvd1[[z]], tvd2[[z]])
    )
  }))
  strat_rows <- bind_rows(lapply(contrast_levels, function(z) {
    bind_rows(
      lapply(levels(lr_depth_stratum), function(s) {
        idx <- which(lr_depth_stratum == s)
        tibble(replicate = rep_id, replicate_seed = seed, contrast = z,
               stratum_type = "long-read abundance quartile", stratum = s,
               n_genes = length(idx), reliability_rho = safe_cor(tvd1[[z]][idx], tvd2[[z]][idx]))
      }),
      lapply(levels(lr_mult_stratum), function(s) {
        idx <- which(lr_mult_stratum == s)
        tibble(replicate = rep_id, replicate_seed = seed, contrast = z,
               stratum_type = "eligible isoform multiplicity", stratum = s,
               n_genes = length(idx), reliability_rho = safe_cor(tvd1[[z]][idx], tvd2[[z]][idx]))
      })
    )
  }))
  list(all = all_rows, stratified = strat_rows)
}

record_progress("LR_RELIABILITY", "Starting 500 complementary binomial pseudo-replicate splits")
lr_checkpoint_objects <- run_checkpointed(
  "lr", seq_len(LR_REPLICATES), lr_worker,
  file.path(DIR_CHECKPOINT, "longread"), batch_size = 18L
)
lr_reliability_all <- bind_rows(lapply(lr_checkpoint_objects, `[[`, "all"))
lr_reliability_strat_all <- bind_rows(lapply(lr_checkpoint_objects, `[[`, "stratified"))
stop_if(nrow(lr_reliability_all) != LR_REPLICATES * 2L || any(!is.finite(lr_reliability_all$reliability_rho)),
        "Long-read reliability replicate table is incomplete.")

sr_library_depths <- tibble(
  sample_id = PRIMARY_SAMPLES,
  common_universe_estimated_total = colSums(common_sr_counts),
  rounded_total = as.integer(round(colSums(common_sr_counts))),
  half1_depth = as.integer(floor(round(colSums(common_sr_counts)) / 2)),
  half2_depth = as.integer(round(colSums(common_sr_counts))) - as.integer(floor(round(colSums(common_sr_counts)) / 2))
)
stop_if(any(abs(sr_library_depths$half1_depth - sr_library_depths$half2_depth) > 1L), "Short-read pseudo-library depths differ by more than one read.")
write_tsv2(sr_library_depths, file.path(DIR_RESULTS, "Fig3_SR_pseudolibrary_depth_audit.tsv"))

fit_sr_half <- function(count_matrix) {
  rownames(count_matrix) <- COMMON_IDS
  colnames(count_matrix) <- PRIMARY_SAMPLES
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = primary_meta, design = ~ patient_id + tissue)
  dds <- suppressMessages(DESeq(dds, quiet = TRUE, parallel = FALSE))
  cbind(
    T_vs_A = results(dds, contrast = c("tissue", "T", "A"), alpha = 0.05)$log2FoldChange,
    T_vs_N = results(dds, contrast = c("tissue", "T", "N"), alpha = 0.05)$log2FoldChange
  )
}

sr_worker <- function(rep_id) {
  seed <- 20260829L + 100000L + rep_id
  set.seed(seed)
  half1 <- matrix(0L, nrow = length(COMMON_IDS), ncol = length(PRIMARY_SAMPLES))
  half2 <- matrix(0L, nrow = length(COMMON_IDS), ncol = length(PRIMARY_SAMPLES))
  for (j in seq_along(PRIMARY_SAMPLES)) {
    probability <- common_sr_counts[, j]
    probability <- probability / sum(probability)
    half1[, j] <- as.vector(rmultinom(1L, sr_library_depths$half1_depth[[j]], probability))
    half2[, j] <- as.vector(rmultinom(1L, sr_library_depths$half2_depth[[j]], probability))
  }
  fit1 <- fit_sr_half(half1)
  fit2 <- fit_sr_half(half2)
  bind_rows(lapply(contrast_levels, function(z) {
    keep <- is.finite(fit1[, z]) & is.finite(fit2[, z])
    tibble(
      replicate = rep_id, replicate_seed = seed, contrast = z,
      n_common_genes = length(COMMON_IDS), n_valid_genes = sum(keep),
      reliability_rho = safe_cor(abs(fit1[, z]), abs(fit2[, z]))
    )
  }))
}

record_progress("SR_RELIABILITY", "Starting 200 paired-DESeq2 multinomial pseudo-library splits")
sr_checkpoint_objects <- run_checkpointed(
  "sr", seq_len(SR_REPLICATES), sr_worker,
  file.path(DIR_CHECKPOINT, "shortread"), batch_size = 12L
)
sr_reliability_all <- bind_rows(sr_checkpoint_objects)
stop_if(nrow(sr_reliability_all) != SR_REPLICATES * 2L || any(!is.finite(sr_reliability_all$reliability_rho)),
        "Short-read reliability replicate table is incomplete.")
stop_if(any(sr_reliability_all$n_valid_genes != 1503L), "A short-read pseudo-replicate did not yield finite effects for all 1,503 genes.")

lr_reliability_summary <- lr_reliability_all %>%
  group_by(contrast) %>%
  group_modify(~ summarise_distribution(.x$reliability_rho)) %>%
  ungroup() %>%
  mutate(
    measurement = "long-read TVD",
    replicate_definition = "complementary Binomial(count, 0.5) pseudo-replicates; not technical replicates"
  ) %>%
  select(measurement, contrast, everything())
lr_reliability_stratified <- lr_reliability_strat_all %>%
  group_by(contrast, stratum_type, stratum, n_genes) %>%
  group_modify(~ summarise_distribution(.x$reliability_rho)) %>%
  ungroup()
sr_reliability_summary <- sr_reliability_all %>%
  group_by(contrast) %>%
  group_modify(~ summarise_distribution(.x$reliability_rho)) %>%
  ungroup() %>%
  mutate(
    measurement = "matched short-read absolute log2FC",
    replicate_definition = "independent equal-depth multinomial pseudo-libraries followed by paired DESeq2; not technical replicates"
  ) %>%
  select(measurement, contrast, everything())

write_tsv2(lr_reliability_all, file.path(DIR_RESULTS, "Fig3_LR_TVD_reliability_all_replicates.tsv"))
write_tsv2(lr_reliability_summary, file.path(DIR_RESULTS, "Fig3_LR_TVD_reliability_summary.tsv"))
write_tsv2(lr_reliability_stratified, file.path(DIR_RESULTS, "Fig3_LR_TVD_reliability_stratified.tsv"))
write_tsv2(sr_reliability_all, file.path(DIR_RESULTS, "Fig3_SR_log2FC_reliability_all_replicates.tsv"))
write_tsv2(sr_reliability_summary, file.path(DIR_RESULTS, "Fig3_SR_log2FC_reliability_summary.tsv"))

record_progress("ADJUSTED_COUPLING", "Computing observed and approximate reliability-adjusted coupling")
reliability_adjusted <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  dat <- covariates %>% filter(contrast == z) %>% arrange(stable_gene_id)
  observed_boot <- bootstrap_cor(
    dat$abs_shortread_log2FC, dat$TVD_usage_change,
    B = 4000L, seed = 20260829L + 200000L + i
  )
  rel_x_values <- sr_reliability_all$reliability_rho[sr_reliability_all$contrast == z]
  rel_y_values <- lr_reliability_all$reliability_rho[lr_reliability_all$contrast == z]
  rel_x <- median(rel_x_values)
  rel_y <- median(rel_y_values)
  stop_if(rel_x <= 0 || rel_y <= 0, paste0("Non-positive median reliability for ", z, "."))
  adjusted <- max(-1, min(1, observed_boot[["estimate"]] / sqrt(rel_x * rel_y)))
  set.seed(20260829L + 210000L + i)
  B <- 4000L
  adjusted_boot <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
    rho_b <- safe_cor(dat$abs_shortread_log2FC[idx], dat$TVD_usage_change[idx])
    rx <- sample(rel_x_values[rel_x_values > 0], 1L)
    ry <- sample(rel_y_values[rel_y_values > 0], 1L)
    adjusted_boot[[b]] <- max(-1, min(1, rho_b / sqrt(rx * ry)))
  }
  stop_if(sum(is.finite(adjusted_boot)) < 3800L, paste0("Insufficient adjusted bootstrap estimates for ", z, "."))
  tibble(
    contrast = z,
    n_genes = nrow(dat),
    rho_observed = observed_boot[["estimate"]],
    observed_P = observed_boot[["p"]],
    observed_CI_low = observed_boot[["ci_low"]],
    observed_CI_high = observed_boot[["ci_high"]],
    rel_X_median_SR_abs_log2FC = rel_x,
    rel_Y_median_LR_TVD = rel_y,
    rho_reliability_adjusted = adjusted,
    adjusted_CI_low = unname(quantile(adjusted_boot, 0.025, na.rm = TRUE)),
    adjusted_CI_high = unname(quantile(adjusted_boot, 0.975, na.rm = TRUE)),
    bootstrap_replicates = B,
    interpretation_boundary = "approximate reliability-adjusted sensitivity estimate; not a true correlation estimate"
  )
}))
write_tsv2(reliability_adjusted, file.path(DIR_RESULTS, "Fig3_reliability_adjusted_coupling.tsv"))

run_stratified_coupling <- function(data, stratum_col, stratum_type, outcome_col, seed_offset) {
  strata <- levels(data[[stratum_col]])
  if (is.null(strata)) strata <- sort(unique(as.character(data[[stratum_col]])))
  bind_rows(lapply(seq_along(strata), function(i) {
    s <- strata[[i]]
    one <- data[as.character(data[[stratum_col]]) == s, , drop = FALSE]
    boot <- bootstrap_cor(
      one$abs_shortread_log2FC, one[[outcome_col]],
      B = 4000L, seed = 20260829L + seed_offset + i
    )
    tibble(
      contrast = unique(data$contrast), outcome = outcome_col,
      stratum_type = stratum_type, stratum = s, n_genes = boot[["n"]],
      rho = boot[["estimate"]], P = boot[["p"]],
      CI_low = boot[["ci_low"]], CI_high = boot[["ci_high"]],
      bootstrap_replicates = 4000L
    )
  }))
}

record_progress("STRATIFIED", "Computing abundance- and multiplicity-stratified coupling")
coupling_lr_abundance <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  run_stratified_coupling(
    covariates %>% filter(contrast == z), "lr_abundance_quartile",
    "long-read mean gene-count quartile", "TVD_usage_change", 220000L + i * 100L
  )
}))
coupling_sr_abundance <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  run_stratified_coupling(
    covariates %>% filter(contrast == z), "sr_abundance_quartile",
    "short-read baseMean quartile", "TVD_usage_change", 230000L + i * 100L
  )
}))
coupling_multiplicity <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  run_stratified_coupling(
    covariates %>% filter(contrast == z), "multiplicity_group",
    "eligible isoform multiplicity", "TVD_usage_change", 240000L + i * 100L
  )
}))
write_tsv2(coupling_lr_abundance, file.path(DIR_RESULTS, "Fig3_coupling_by_LR_abundance.tsv"))
write_tsv2(coupling_sr_abundance, file.path(DIR_RESULTS, "Fig3_coupling_by_SR_abundance.tsv"))
write_tsv2(coupling_multiplicity, file.path(DIR_RESULTS, "Fig3_coupling_by_isoform_multiplicity.tsv"))

partial_once <- function(data, outcome_col, covariate_set) {
  z <- data.frame(
    X = normal_score(data$abs_shortread_log2FC),
    Y = normal_score(data[[outcome_col]])
  )
  for (nm in covariate_set) z[[nm]] <- normal_score(data[[nm]])
  design <- model.matrix(reformulate(covariate_set), z)
  residual_x <- lm.fit(design, z$X)$residuals
  residual_y <- lm.fit(design, z$Y)$residuals
  c(rho = cor(residual_x, residual_y), P = suppressWarnings(cor.test(residual_x, residual_y)$p.value))
}
bootstrap_partial <- function(data, outcome_col, covariate_set, B, seed) {
  observed <- partial_once(data, outcome_col, covariate_set)
  set.seed(seed)
  vals <- replicate(B, {
    idx <- sample.int(nrow(data), nrow(data), replace = TRUE)
    partial_once(data[idx, , drop = FALSE], outcome_col, covariate_set)[["rho"]]
  })
  vals <- vals[is.finite(vals)]
  stop_if(length(vals) < floor(0.95 * B), "Too few valid partial-correlation bootstrap values.")
  c(
    rho = observed[["rho"]], P = observed[["P"]],
    CI_low = unname(quantile(vals, 0.025)), CI_high = unname(quantile(vals, 0.975)),
    bootstrap_valid = length(vals), n = nrow(data)
  )
}

partial_covariate_sets <- list(
  primary = c("log_lr_count", "log_sr_baseMean", "isoform_multiplicity", "zero_fraction"),
  sensitivity = c("log_lr_count", "isoform_multiplicity")
)
prepare_adjusted_data <- function(data) {
  data %>% mutate(
    log_lr_count = log10(lr_mean_gene_count + 1),
    log_sr_baseMean = log10(shortread_baseMean + 1)
  )
}

record_progress("PARTIAL", "Computing rank-residual partial correlations with 4,000 gene bootstraps")
partial_correlation <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  dat <- prepare_adjusted_data(covariates %>% filter(contrast == z))
  bind_rows(lapply(seq_along(partial_covariate_sets), function(j) {
    model_name <- names(partial_covariate_sets)[[j]]
    covs <- partial_covariate_sets[[j]]
    out <- bootstrap_partial(dat, "TVD_usage_change", covs, 4000L, 20260829L + 250000L + i * 100L + j)
    tibble(
      contrast = z, outcome = "TVD_usage_change", model = model_name,
      covariates = paste(covs, collapse = ";"), n_genes = out[["n"]],
      partial_rho = out[["rho"]], P = out[["P"]],
      CI_low = out[["CI_low"]], CI_high = out[["CI_high"]],
      bootstrap_replicates = 4000L, rank_method = "normal-score ranks with Pearson correlation of covariate-adjusted rank residuals"
    )
  }))
}))
write_tsv2(partial_correlation, file.path(DIR_RESULTS, "Fig3_partial_correlation.tsv"))

cv_explained_variance <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  denominator <- sum((observed - mean(observed))^2)
  if (length(observed) < 4L || denominator <= 0) return(NA_real_)
  1 - sum((observed - predicted)^2) / denominator
}

prepare_incremental_data <- function(data, outcome_col, model_set) {
  n <- nrow(data)
  out <- data.frame(
    stable_gene_id = data$stable_gene_id,
    abs_shortread_log2FC = data$abs_shortread_log2FC,
    log_lr_count = log10(data$lr_mean_gene_count + 1),
    log_sr_baseMean = log10(data$shortread_baseMean + 1),
    isoform_multiplicity = data$isoform_multiplicity,
    zero_fraction = data$zero_fraction
  )
  if (model_set == "rank_linear") {
    out$Y <- normal_score(data[[outcome_col]])
  } else if (model_set == "quasibinomial") {
    out$Y <- (data[[outcome_col]] * (n - 1) + 0.5) / n
  } else {
    stop("Unknown model set: ", model_set)
  }
  out
}

fit_incremental_models <- function(model_data, model_set) {
  baseline_formula <- Y ~ log_lr_count + log_sr_baseMean + isoform_multiplicity + zero_fraction
  extended_formula <- update(baseline_formula, . ~ . + abs_shortread_log2FC)
  if (model_set == "rank_linear") {
    baseline <- lm(baseline_formula, data = model_data)
    extended <- lm(extended_formula, data = model_data)
    baseline_metric <- summary(baseline)$adj.r.squared
    extended_metric <- summary(extended)$adj.r.squared
    metric_name <- "adjusted_R2"
    added_effect_P <- coef(summary(extended))["abs_shortread_log2FC", "Pr(>|t|)"]
  } else {
    baseline <- suppressWarnings(glm(baseline_formula, data = model_data, family = quasibinomial(link = "logit")))
    extended <- suppressWarnings(glm(extended_formula, data = model_data, family = quasibinomial(link = "logit")))
    baseline_metric <- 1 - baseline$deviance / baseline$null.deviance
    extended_metric <- 1 - extended$deviance / extended$null.deviance
    metric_name <- "explained_deviance"
    added_effect_P <- coef(summary(extended))["abs_shortread_log2FC", "Pr(>|t|)"]
  }
  list(
    baseline = baseline, extended = extended, metric_name = metric_name,
    baseline_metric = baseline_metric, extended_metric = extended_metric,
    delta_metric = extended_metric - baseline_metric, added_effect_P = added_effect_P
  )
}

cross_validate_incremental <- function(data, outcome_col, model_set, repeats = 100L, seed_offset = 0L) {
  model_data <- prepare_incremental_data(data, outcome_col, model_set)
  bind_rows(lapply(seq_len(repeats), function(rep_id) {
    seed <- 20260829L + seed_offset + rep_id
    set.seed(seed)
    folds <- sample(rep(seq_len(10L), length.out = nrow(model_data)))
    pred_baseline <- rep(NA_real_, nrow(model_data))
    pred_extended <- rep(NA_real_, nrow(model_data))
    for (fold in seq_len(10L)) {
      train <- model_data[folds != fold, , drop = FALSE]
      test <- model_data[folds == fold, , drop = FALSE]
      fits <- fit_incremental_models(train, model_set)
      pred_baseline[folds == fold] <- predict(fits$baseline, newdata = test, type = "response")
      pred_extended[folds == fold] <- predict(fits$extended, newdata = test, type = "response")
    }
    baseline_ev <- cv_explained_variance(model_data$Y, pred_baseline)
    extended_ev <- cv_explained_variance(model_data$Y, pred_extended)
    tibble(
      cv_repeat = rep_id, repeat_seed = seed, n_genes = nrow(model_data),
      outcome = outcome_col, model_set = model_set, folds = 10L,
      baseline_CV_explained_variance = baseline_ev,
      extended_CV_explained_variance = extended_ev,
      delta_CV_explained_variance = extended_ev - baseline_ev
    )
  }))
}

run_incremental_analysis <- function(data, outcome_col, contrast, universe_label, seed_offset, model_sets = c("rank_linear", "quasibinomial")) {
  all_cv <- list()
  summary_rows <- list()
  for (j in seq_along(model_sets)) {
    model_set <- model_sets[[j]]
    model_data <- prepare_incremental_data(data, outcome_col, model_set)
    apparent <- fit_incremental_models(model_data, model_set)
    cv <- cross_validate_incremental(data, outcome_col, model_set, 100L, seed_offset + j * 1000L) %>%
      mutate(contrast = contrast, universe = universe_label, .before = 1)
    all_cv[[model_set]] <- cv
    summary_rows[[model_set]] <- tibble(
      contrast = contrast, universe = universe_label, outcome = outcome_col, model_set = model_set,
      n_genes = nrow(data), apparent_metric = apparent$metric_name,
      baseline_apparent = apparent$baseline_metric,
      extended_apparent = apparent$extended_metric,
      delta_apparent = apparent$delta_metric,
      added_abs_log2FC_P = apparent$added_effect_P,
      baseline_CV_median = median(cv$baseline_CV_explained_variance),
      baseline_CV_CI_low = unname(quantile(cv$baseline_CV_explained_variance, 0.025)),
      baseline_CV_CI_high = unname(quantile(cv$baseline_CV_explained_variance, 0.975)),
      extended_CV_median = median(cv$extended_CV_explained_variance),
      extended_CV_CI_low = unname(quantile(cv$extended_CV_explained_variance, 0.025)),
      extended_CV_CI_high = unname(quantile(cv$extended_CV_explained_variance, 0.975)),
      delta_CV_median = median(cv$delta_CV_explained_variance),
      delta_CV_CI_low = unname(quantile(cv$delta_CV_explained_variance, 0.025)),
      delta_CV_CI_high = unname(quantile(cv$delta_CV_explained_variance, 0.975)),
      CV_repeats = 100L, CV_folds = 10L
    )
  }
  list(summary = bind_rows(summary_rows), cv = bind_rows(all_cv))
}

record_progress("INCREMENTAL", "Running rank-linear and quasibinomial incremental models with repeated 10-fold CV")
incremental_objects <- lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  dat <- covariates %>% filter(contrast == z)
  run_incremental_analysis(dat, "TVD_usage_change", z, "all_common_genes", 300000L + i * 10000L)
})
incremental_summary <- bind_rows(lapply(incremental_objects, `[[`, "summary"))
incremental_cv <- bind_rows(lapply(incremental_objects, `[[`, "cv"))

record_progress("HIGH_CONFIDENCE", "Applying pre-specified quality subset without outcome-based tuning")
high_confidence_membership <- covariates %>%
  group_by(contrast) %>%
  mutate(
    lr_mean_gene_count_median = median(lr_mean_gene_count),
    zero_fraction_median = median(zero_fraction),
    high_confidence = lr_mean_gene_count >= lr_mean_gene_count_median &
      isoform_multiplicity <= 5 & zero_fraction <= zero_fraction_median,
    high_confidence_definition = "LR mean gene count >= common-universe median; eligible isoform multiplicity <=5; zero fraction <= common-universe median"
  ) %>%
  ungroup()
write_tsv2(high_confidence_membership, file.path(DIR_RESULTS, "Fig3_high_confidence_subset.tsv"))

high_confidence_summary_rows <- list()
high_confidence_cv_rows <- list()
for (i in seq_along(contrast_levels)) {
  z <- contrast_levels[[i]]
  dat <- high_confidence_membership %>% filter(contrast == z, high_confidence) %>% prepare_adjusted_data()
  stop_if(nrow(dat) < 100L, paste0("High-confidence subset is unexpectedly small for ", z, "."))
  raw <- bootstrap_cor(dat$abs_shortread_log2FC, dat$TVD_usage_change, 4000L, 20260829L + 320000L + i)
  partial <- bootstrap_partial(dat, "TVD_usage_change", partial_covariate_sets$primary, 4000L, 20260829L + 321000L + i)
  inc <- run_incremental_analysis(
    dat, "TVD_usage_change", z, "high_confidence", 330000L + i * 10000L,
    model_sets = "rank_linear"
  )
  high_confidence_cv_rows[[z]] <- inc$cv
  high_confidence_summary_rows[[z]] <- tibble(
    contrast = z, n_high_confidence_genes = nrow(dat),
    observed_rho = raw[["estimate"]], observed_P = raw[["p"]],
    observed_CI_low = raw[["ci_low"]], observed_CI_high = raw[["ci_high"]],
    partial_rho = partial[["rho"]], partial_P = partial[["P"]],
    partial_CI_low = partial[["CI_low"]], partial_CI_high = partial[["CI_high"]],
    incremental_baseline_CV_median = inc$summary$baseline_CV_median,
    incremental_extended_CV_median = inc$summary$extended_CV_median,
    incremental_delta_CV_median = inc$summary$delta_CV_median
  )
}
high_confidence_summary <- bind_rows(high_confidence_summary_rows)
high_confidence_cv <- bind_rows(high_confidence_cv_rows)
write_tsv2(high_confidence_summary, file.path(DIR_RESULTS, "Fig3_high_confidence_summary.tsv"))

record_progress("MAX_DELTA", "Repeating non-reliability analyses with max absolute delta usage")
maxdelta_raw <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  dat <- covariates %>% filter(contrast == z)
  out <- bootstrap_cor(dat$abs_shortread_log2FC, dat$max_abs_delta_usage, 4000L, 20260829L + 350000L + i)
  tibble(
    analysis = "raw correlation", contrast = z, stratum_type = "all", stratum = "all",
    n_genes = out[["n"]], estimate = out[["estimate"]], P = out[["p"]],
    CI_low = out[["ci_low"]], CI_high = out[["ci_high"]], detail = "Spearman rho"
  )
}))
maxdelta_lr <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  run_stratified_coupling(covariates %>% filter(contrast == z), "lr_abundance_quartile",
                          "long-read mean gene-count quartile", "max_abs_delta_usage", 351000L + i * 100L)
})) %>% transmute(analysis = "LR abundance stratification", contrast, stratum_type, stratum,
                  n_genes, estimate = rho, P, CI_low, CI_high, detail = "Spearman rho")
maxdelta_sr <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  run_stratified_coupling(covariates %>% filter(contrast == z), "sr_abundance_quartile",
                          "short-read baseMean quartile", "max_abs_delta_usage", 352000L + i * 100L)
})) %>% transmute(analysis = "SR abundance stratification", contrast, stratum_type, stratum,
                  n_genes, estimate = rho, P, CI_low, CI_high, detail = "Spearman rho")
maxdelta_mult <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  run_stratified_coupling(covariates %>% filter(contrast == z), "multiplicity_group",
                          "eligible isoform multiplicity", "max_abs_delta_usage", 353000L + i * 100L)
})) %>% transmute(analysis = "multiplicity stratification", contrast, stratum_type, stratum,
                  n_genes, estimate = rho, P, CI_low, CI_high, detail = "Spearman rho")
maxdelta_partial <- bind_rows(lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  dat <- prepare_adjusted_data(covariates %>% filter(contrast == z))
  bind_rows(lapply(seq_along(partial_covariate_sets), function(j) {
    model_name <- names(partial_covariate_sets)[[j]]
    out <- bootstrap_partial(dat, "max_abs_delta_usage", partial_covariate_sets[[j]], 4000L, 20260829L + 354000L + i * 100L + j)
    tibble(
      analysis = "partial correlation", contrast = z, stratum_type = "model", stratum = model_name,
      n_genes = out[["n"]], estimate = out[["rho"]], P = out[["P"]],
      CI_low = out[["CI_low"]], CI_high = out[["CI_high"]], detail = "rank-residual partial rho"
    )
  }))
}))
maxdelta_incremental_objects <- lapply(seq_along(contrast_levels), function(i) {
  z <- contrast_levels[[i]]
  run_incremental_analysis(
    covariates %>% filter(contrast == z), "max_abs_delta_usage", z,
    "all_common_genes", 360000L + i * 10000L
  )
})
maxdelta_incremental_summary <- bind_rows(lapply(maxdelta_incremental_objects, `[[`, "summary"))
maxdelta_incremental_cv <- bind_rows(lapply(maxdelta_incremental_objects, `[[`, "cv"))
maxdelta_incremental_rows <- maxdelta_incremental_summary %>%
  transmute(
    analysis = "incremental CV", contrast, stratum_type = "model_set", stratum = model_set,
    n_genes, estimate = delta_CV_median, P = added_abs_log2FC_P,
    CI_low = delta_CV_CI_low, CI_high = delta_CV_CI_high,
    detail = "delta CV explained variance after adding absolute short-read log2FC"
  )
maxdelta_sensitivity <- bind_rows(
  maxdelta_raw, maxdelta_lr, maxdelta_sr, maxdelta_mult, maxdelta_partial, maxdelta_incremental_rows
)
write_tsv2(maxdelta_sensitivity, file.path(DIR_RESULTS, "Fig3_maxdelta_sensitivity.tsv"))

incremental_cv_all <- bind_rows(incremental_cv, high_confidence_cv)
write_tsv2(incremental_summary, file.path(DIR_RESULTS, "Fig3_incremental_explanatory_value.tsv"))
write_tsv2(incremental_cv_all, file.path(DIR_RESULTS, "Fig3_incremental_CV_all_replicates.tsv"))

record_progress("LEVEL", "Applying the pre-specified LEVEL A/B/C decision rules")
primary_partial <- partial_correlation %>% filter(model == "primary")
primary_incremental <- incremental_summary %>%
  filter(universe == "all_common_genes", model_set == "rank_linear")
high_quality_strata <- bind_rows(
  coupling_lr_abundance %>% filter(stratum == "Q4") %>% mutate(quality_layer = "LR abundance Q4"),
  coupling_sr_abundance %>% filter(stratum == "Q4") %>% mutate(quality_layer = "SR abundance Q4"),
  coupling_multiplicity %>% filter(stratum %in% c("2", "3", "4-5"), n_genes >= 100L) %>%
    mutate(quality_layer = paste0("multiplicity ", stratum))
)

decision_values <- tibble(
  contrast = contrast_levels,
  LR_reliability = lr_reliability_summary$median[match(contrast_levels, lr_reliability_summary$contrast)],
  SR_reliability = sr_reliability_summary$median[match(contrast_levels, sr_reliability_summary$contrast)],
  observed_rho = reliability_adjusted$rho_observed[match(contrast_levels, reliability_adjusted$contrast)],
  adjusted_rho = reliability_adjusted$rho_reliability_adjusted[match(contrast_levels, reliability_adjusted$contrast)],
  partial_rho = primary_partial$partial_rho[match(contrast_levels, primary_partial$contrast)],
  high_confidence_rho = high_confidence_summary$observed_rho[match(contrast_levels, high_confidence_summary$contrast)],
  delta_CV = primary_incremental$delta_CV_median[match(contrast_levels, primary_incremental$contrast)]
)
max_high_quality_rho <- max(abs(high_quality_strata$rho), na.rm = TRUE)

level_a_conditions <- c(
  LR_reliability_ge_0.75 = all(decision_values$LR_reliability >= 0.75),
  SR_reliability_ge_0.80 = all(decision_values$SR_reliability >= 0.80),
  observed_rho_le_0.20 = all(abs(decision_values$observed_rho) <= 0.20),
  adjusted_rho_le_0.25 = all(abs(decision_values$adjusted_rho) <= 0.25),
  partial_rho_le_0.25 = all(abs(decision_values$partial_rho) <= 0.25),
  high_confidence_rho_le_0.30 = all(abs(decision_values$high_confidence_rho) <= 0.30),
  delta_CV_lt_0.05 = all(decision_values$delta_CV < 0.05),
  no_high_quality_stratum_rho_ge_0.50 = max_high_quality_rho < 0.50
)
level_c_conditions <- c(
  LR_reliability_lt_0.50 = any(decision_values$LR_reliability < 0.50),
  high_quality_stratum_rho_ge_0.50 = max_high_quality_rho >= 0.50,
  adjusted_rho_gt_0.40 = any(abs(decision_values$adjusted_rho) > 0.40),
  delta_CV_ge_0.10 = any(decision_values$delta_CV >= 0.10)
)
FINAL_LEVEL <- if (all(level_a_conditions)) "LEVEL A" else if (any(level_c_conditions)) "LEVEL C" else "LEVEL B"
level_interpretation <- switch(
  FINAL_LEVEL,
  "LEVEL A" = "strong support for limited cross-layer coupling",
  "LEVEL B" = "partially overlapping, non-equivalent molecular layers",
  "LEVEL C" = "measurement noise or confounding materially contributes; weak coupling cannot be the primary claim"
)
level_decision <- bind_rows(
  tibble(level = "LEVEL A", criterion = names(level_a_conditions), pass = unname(level_a_conditions)),
  tibble(level = "LEVEL C trigger", criterion = names(level_c_conditions), pass = unname(level_c_conditions))
) %>% mutate(final_level = FINAL_LEVEL, interpretation = level_interpretation)
write_tsv2(decision_values, file.path(DIR_RESULTS, "Fig3_robustness_decision_values.tsv"))
write_tsv2(level_decision, file.path(DIR_RESULTS, "Fig3_robustness_level_decision.tsv"))

maxdelta_raw_check <- maxdelta_raw %>% arrange(match(contrast, contrast_levels))
maxdelta_partial_primary <- maxdelta_partial %>% filter(stratum == "primary") %>% arrange(match(contrast, contrast_levels))
MAXDELTA_DIRECTION_CONSISTENT <- all(sign(maxdelta_raw_check$estimate) == sign(decision_values$observed_rho)) &&
  all(sign(maxdelta_partial_primary$estimate) == sign(decision_values$partial_rho))

FONT_FAMILY <- "Times New Roman"
font_match <- systemfonts::match_fonts(FONT_FAMILY)[1, , drop = FALSE]
stop_if(is.null(font_match$path) || !nzchar(font_match$path) || !file.exists(font_match$path),
        "Times New Roman is unavailable in the remote R graphics environment.")
palette <- c(
  T_vs_A = "#E66B16", T_vs_N = "#7B4EA3",
  LR = "#3A7095", SR = "#B45A53", raw = "#777777",
  adjusted = "#3A7095", partial = "#6A9F58", baseline = "#B9BEC2", extended = "#D77A2A"
)
theme_robustness <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY, colour = "#202020"),
      axis.line = element_line(linewidth = 0.28, colour = "#333333"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#333333"),
      plot.title = element_text(face = "bold", size = base_size + 0.8, margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#555555", margin = margin(b = 4)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_blank(), legend.key.height = grid::unit(3, "mm"),
      panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.18),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
save_bundle <- function(plot, base_path, width_mm, height_mm, dpi = 600L) {
  paths <- c(pdf = paste0(base_path, ".pdf"), svg = paste0(base_path, ".svg"), png = paste0(base_path, ".png"))
  grDevices::cairo_pdf(paths[["pdf"]], width = width_mm / 25.4, height = height_mm / 25.4, family = FONT_FAMILY)
  print(plot)
  dev.off()
  grDevices::svg(paths[["svg"]], width = width_mm / 25.4, height = height_mm / 25.4,
                 family = FONT_FAMILY, onefile = TRUE, bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(paths[["png"]], width = width_mm, height = height_mm, units = "mm", res = dpi, background = "white")
  print(plot)
  dev.off()
  stop_if(any(!file.exists(paths) | file.info(paths)$size <= 0), paste0("Figure export failed for ", base_path))
  tibble(format = names(paths), path = unname(paths), size_bytes = file.info(paths)$size)
}

candidate_paths_expected <- c(
  file.path(DIR_FIGURES, paste0("Figure3_robustness_candidate.", c("pdf", "svg", "png"))),
  file.path(DIR_FIGURES, paste0("Supplementary_Fig3_reliability.", c("pdf", "svg", "png"))),
  file.path(DIR_FIGURES, paste0("Supplementary_Fig3_coupling_by_abundance.", c("pdf", "svg", "png"))),
  file.path(DIR_FIGURES, paste0("Supplementary_Fig3_coupling_by_multiplicity.", c("pdf", "svg", "png"))),
  file.path(DIR_FIGURES, paste0("Supplementary_Fig3_incremental_explanatory_value.", c("pdf", "svg", "png"))),
  file.path(DIR_FIGURES, "Fig3_incremental_explanatory_value.pdf")
)
archive_path <- archive_existing(candidate_paths_expected)

record_progress("PLOT", "Rendering internal-review candidate and supplementary robustness figures")
reliability_plot_data <- bind_rows(
  lr_reliability_summary %>% transmute(
    contrast, measurement = "LR TVD", median, low = percentile_2.5, high = percentile_97.5, platform = "LR"
  ),
  sr_reliability_summary %>% transmute(
    contrast, measurement = "SR |log2FC|", median, low = percentile_2.5, high = percentile_97.5, platform = "SR"
  )
) %>%
  mutate(
    comparison = factor(ifelse(contrast == "T_vs_A", "T vs A", "T vs N"), levels = c("T vs A", "T vs N")),
    measurement = factor(measurement, levels = c("LR TVD", "SR |log2FC|"))
  )
panel_a <- ggplot(reliability_plot_data, aes(measurement, median, colour = platform)) +
  geom_hline(yintercept = c(0.75, 0.80), colour = "#B0B0B0", linetype = "dashed", linewidth = 0.35) +
  geom_linerange(aes(ymin = low, ymax = high), linewidth = 0.65) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", median)), colour = "#202020", vjust = -1.1, size = 2.2) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_colour_manual(values = palette[c("LR", "SR")], guide = "none") +
  scale_y_continuous(limits = c(0, 1.04), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "A  Pseudo-replicate reliability",
    subtitle = "Median and 2.5-97.5 percentiles; pseudo-replicates are count-resampling constructs",
    x = NULL, y = "Reliability (Spearman rho)"
  ) + theme_robustness()

partial_for_plot <- primary_partial %>%
  transmute(
    contrast, method = "Partial", estimate = partial_rho,
    low = CI_low, high = CI_high
  )
coupling_plot_data <- bind_rows(
  reliability_adjusted %>% transmute(contrast, method = "Observed", estimate = rho_observed, low = observed_CI_low, high = observed_CI_high),
  reliability_adjusted %>% transmute(contrast, method = "Reliability-adjusted", estimate = rho_reliability_adjusted, low = adjusted_CI_low, high = adjusted_CI_high),
  partial_for_plot
) %>%
  mutate(
    comparison = factor(ifelse(contrast == "T_vs_A", "T vs A", "T vs N"), levels = c("T vs A", "T vs N")),
    method = factor(method, levels = c("Observed", "Reliability-adjusted", "Partial"))
  )
panel_b <- ggplot(coupling_plot_data, aes(method, estimate, colour = method)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
  geom_hline(yintercept = 0.25, colour = "#B0B0B0", linetype = "dashed", linewidth = 0.35) +
  geom_linerange(aes(ymin = low, ymax = high), linewidth = 0.65) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", estimate)), colour = "#202020", vjust = -1.1, size = 2.2) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_colour_manual(values = c(Observed = palette[["raw"]], `Reliability-adjusted` = palette[["adjusted"]], Partial = palette[["partial"]]), guide = "none") +
  scale_y_continuous(limits = c(-0.05, max(0.46, max(coupling_plot_data$high) + 0.05)), breaks = seq(0, 0.5, 0.1)) +
  labs(
    title = "B  Cross-layer coupling after robustness adjustment",
    subtitle = "Reliability adjustment is an approximate sensitivity estimate",
    x = NULL, y = "Correlation estimate"
  ) + theme_robustness() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

incremental_plot_data <- incremental_summary %>%
  filter(universe == "all_common_genes", model_set == "rank_linear") %>%
  select(contrast, baseline_CV_median, baseline_CV_CI_low, baseline_CV_CI_high,
         extended_CV_median, extended_CV_CI_low, extended_CV_CI_high, delta_CV_median) %>%
  pivot_longer(
    cols = starts_with(c("baseline_CV", "extended_CV")),
    names_pattern = "(baseline|extended)_CV_(median|CI_low|CI_high)",
    names_to = c("model", ".value")
  ) %>%
  mutate(
    comparison = factor(ifelse(contrast == "T_vs_A", "T vs A", "T vs N"), levels = c("T vs A", "T vs N")),
    model = factor(model, levels = c("baseline", "extended"), labels = c("Baseline", "+ |gene log2FC|"))
  )
panel_c <- ggplot(incremental_plot_data, aes(model, median, colour = model)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high), linewidth = 0.65) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", median)), colour = "#202020", vjust = -1.1, size = 2.2) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_colour_manual(values = c(Baseline = palette[["baseline"]], `+ |gene log2FC|` = palette[["extended"]]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
  labs(
    title = "C  Incremental cross-validated explanatory value",
    subtitle = "Repeated 10-fold CV; rank-linear model; error bars show 2.5-97.5 percentiles",
    x = NULL, y = "CV explained variance"
  ) + theme_robustness() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

strata_plot_data <- bind_rows(
  coupling_lr_abundance %>% mutate(display_group = paste0("LR ", stratum), family = "LR abundance"),
  coupling_multiplicity %>% mutate(display_group = paste0("Isoforms ", stratum), family = "Multiplicity")
) %>%
  mutate(
    comparison = factor(ifelse(contrast == "T_vs_A", "T vs A", "T vs N"), levels = c("T vs A", "T vs N")),
    display_group = factor(display_group, levels = c(paste0("LR Q", 1:4), "Isoforms 2", "Isoforms 3", "Isoforms 4-5", "Isoforms >=6"))
  )
panel_d <- ggplot(strata_plot_data, aes(display_group, rho, colour = comparison)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
  geom_hline(yintercept = 0.50, colour = "#B0B0B0", linetype = "dashed", linewidth = 0.35) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high), linewidth = 0.55, position = position_dodge(width = 0.4)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.4)) +
  facet_grid(. ~ family, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c(`T vs A` = palette[["T_vs_A"]], `T vs N` = palette[["T_vs_N"]])) +
  labs(
    title = "D  Coupling across abundance and isoform multiplicity strata",
    subtitle = "Points are Spearman rho; intervals are 4,000-gene-bootstrap 95% CIs",
    x = NULL, y = "Stratified rho", colour = NULL
  ) + theme_robustness() +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "top")

candidate_figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(heights = c(1, 1.08)) +
  plot_annotation(
    title = "Figure 3 robustness candidate | Measurement reliability and adjusted cross-layer coupling",
    subtitle = paste0("Internal review only; formal Figure 3 is unchanged. Pre-specified decision: ", FINAL_LEVEL, " (", level_interpretation, ")"),
    caption = "Primary72 matched cohort; 1,503 common multi-isoform genes per tumor contrast. All intervals and resampling units are defined in Source Data.",
    theme = theme(
      text = element_text(family = FONT_FAMILY),
      plot.title = element_text(face = "bold", size = 11.2),
      plot.subtitle = element_text(size = 7.8, colour = "#555555"),
      plot.caption = element_text(size = 6.1, colour = "#555555", hjust = 0)
    )
  )
candidate_paths <- save_bundle(candidate_figure, file.path(DIR_FIGURES, "Figure3_robustness_candidate"), 183, 155)

reliability_distribution <- bind_rows(
  lr_reliability_all %>% transmute(contrast, measurement = "LR TVD", replicate, reliability_rho),
  sr_reliability_all %>% transmute(contrast, measurement = "SR |log2FC|", replicate, reliability_rho)
) %>%
  mutate(comparison = ifelse(contrast == "T_vs_A", "T vs A", "T vs N"))
p_rel_1 <- ggplot(reliability_distribution, aes(measurement, reliability_rho, fill = measurement)) +
  geom_violin(width = 0.85, alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = c(`LR TVD` = palette[["LR"]], `SR |log2FC|` = palette[["SR"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(title = "A  Overall pseudo-replicate reliability", x = NULL, y = "Spearman rho") + theme_robustness()
p_rel_2 <- ggplot(lr_reliability_stratified %>% filter(stratum_type == "long-read abundance quartile"),
                  aes(stratum, median, colour = contrast)) +
  geom_linerange(aes(ymin = percentile_2.5, ymax = percentile_97.5), linewidth = 0.6, position = position_dodge(width = 0.35)) +
  geom_point(size = 2.3, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = palette[c("T_vs_A", "T_vs_N")], labels = c("T vs A", "T vs N")) +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(title = "B  Long-read TVD reliability by abundance", x = "LR mean gene-count quartile", y = "Median reliability", colour = NULL) +
  theme_robustness() + theme(legend.position = "top")
supp_reliability <- p_rel_1 | p_rel_2
supp_reliability_paths <- save_bundle(supp_reliability, file.path(DIR_FIGURES, "Supplementary_Fig3_reliability"), 183, 88)

abundance_plot_data <- bind_rows(coupling_lr_abundance, coupling_sr_abundance) %>%
  mutate(
    comparison = ifelse(contrast == "T_vs_A", "T vs A", "T vs N"),
    abundance_metric = ifelse(grepl("long-read", stratum_type), "Long-read mean gene count", "Short-read baseMean")
  )
supp_abundance <- ggplot(abundance_plot_data, aes(stratum, rho, colour = comparison, group = comparison)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high), linewidth = 0.55, position = position_dodge(width = 0.35)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.35)) +
  facet_wrap(~ abundance_metric, nrow = 1) +
  scale_colour_manual(values = c(`T vs A` = palette[["T_vs_A"]], `T vs N` = palette[["T_vs_N"]])) +
  labs(
    title = "Supplementary Figure | Cross-layer coupling by expression depth",
    subtitle = "Spearman rho with 4,000-gene-bootstrap 95% confidence intervals",
    x = "Within-universe abundance quartile", y = "Spearman rho", colour = NULL
  ) + theme_robustness(7.5) + theme(legend.position = "top")
supp_abundance_paths <- save_bundle(supp_abundance, file.path(DIR_FIGURES, "Supplementary_Fig3_coupling_by_abundance"), 183, 82)

supp_multiplicity <- ggplot(coupling_multiplicity %>% mutate(
                              comparison = ifelse(contrast == "T_vs_A", "T vs A", "T vs N"),
                              stratum = factor(stratum, levels = c("2", "3", "4-5", ">=6"))
                            ),
                            aes(stratum, rho, colour = comparison, group = comparison)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high), linewidth = 0.55, position = position_dodge(width = 0.35)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = c(`T vs A` = palette[["T_vs_A"]], `T vs N` = palette[["T_vs_N"]])) +
  labs(
    title = "Supplementary Figure | Cross-layer coupling by isoform multiplicity",
    subtitle = "Eligible usage-analysis isoforms per gene; Spearman rho with 95% bootstrap confidence intervals",
    x = "Eligible isoform multiplicity", y = "Spearman rho", colour = NULL
  ) + theme_robustness(7.5) + theme(legend.position = "top")
supp_multiplicity_paths <- save_bundle(supp_multiplicity, file.path(DIR_FIGURES, "Supplementary_Fig3_coupling_by_multiplicity"), 120, 82)

incremental_cv_plot <- incremental_cv %>%
  select(contrast, model_set, cv_repeat, baseline_CV_explained_variance, extended_CV_explained_variance) %>%
  pivot_longer(cols = c(baseline_CV_explained_variance, extended_CV_explained_variance), names_to = "model", values_to = "CV_explained_variance") %>%
  mutate(
    comparison = ifelse(contrast == "T_vs_A", "T vs A", "T vs N"),
    model = factor(model, levels = c("baseline_CV_explained_variance", "extended_CV_explained_variance"), labels = c("Baseline", "+ |gene log2FC|")),
    model_set = factor(model_set, levels = c("rank_linear", "quasibinomial"), labels = c("Rank-linear", "Quasibinomial sensitivity"))
  )
supp_incremental <- ggplot(incremental_cv_plot, aes(model, CV_explained_variance, fill = model)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.35) +
  facet_grid(model_set ~ comparison, scales = "free_y") +
  scale_fill_manual(values = c(Baseline = palette[["baseline"]], `+ |gene log2FC|` = palette[["extended"]]), guide = "none") +
  labs(
    title = "Supplementary Figure | Incremental explanatory value of matched gene-abundance change",
    subtitle = "One hundred repeated 10-fold cross-validations; each box contains repeat-level held-out explained variance",
    x = NULL, y = "CV explained variance"
  ) + theme_robustness(7.2) + theme(axis.text.x = element_text(angle = 18, hjust = 1))
supp_incremental_paths <- save_bundle(supp_incremental, file.path(DIR_FIGURES, "Supplementary_Fig3_incremental_explanatory_value"), 183, 115)
stop_if(!file.copy(
  file.path(DIR_FIGURES, "Supplementary_Fig3_incremental_explanatory_value.pdf"),
  file.path(DIR_FIGURES, "Fig3_incremental_explanatory_value.pdf"), overwrite = TRUE
), "Could not write the requested incremental explanatory-value PDF alias.")

fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 1e-4, "< 0.0001", paste0("= ", formatC(x, format = "g", digits = 3))))
}
lookup_one <- function(data, contrast, column, ...) {
  x <- data %>% filter(.data$contrast == .env$contrast, ...) %>% pull(all_of(column))
  stop_if(length(x) != 1L, paste0("Expected one lookup value for ", contrast, "/", column, "."))
  x[[1]]
}
contrast_label <- function(x) ifelse(x == "T_vs_A", "T vs A", "T vs N")

key_statistics <- bind_rows(lapply(contrast_levels, function(z) {
  tibble(
    contrast = z,
    n_common_genes = 1503L,
    LR_TVD_reliability_median = lookup_one(lr_reliability_summary, z, "median"),
    SR_abs_log2FC_reliability_median = lookup_one(sr_reliability_summary, z, "median"),
    observed_rho = lookup_one(reliability_adjusted, z, "rho_observed"),
    reliability_adjusted_rho = lookup_one(reliability_adjusted, z, "rho_reliability_adjusted"),
    primary_partial_rho = lookup_one(primary_partial, z, "partial_rho"),
    high_confidence_rho = lookup_one(high_confidence_summary, z, "observed_rho"),
    rank_linear_baseline_CV = lookup_one(primary_incremental, z, "baseline_CV_median"),
    rank_linear_extended_CV = lookup_one(primary_incremental, z, "extended_CV_median"),
    rank_linear_delta_CV = lookup_one(primary_incremental, z, "delta_CV_median"),
    maxdelta_raw_rho = maxdelta_raw$estimate[maxdelta_raw$contrast == z],
    maxdelta_primary_partial_rho = maxdelta_partial_primary$estimate[maxdelta_partial_primary$contrast == z],
    final_level = FINAL_LEVEL,
    level_interpretation = level_interpretation
  )
}))
write_tsv2(key_statistics, file.path(DIR_RESULTS, "Fig3_robustness_key_statistics.tsv"))

lr_q_lookup <- function(z, q) {
  lr_reliability_stratified %>%
    filter(contrast == z, stratum_type == "long-read abundance quartile", stratum == q) %>%
    pull(median)
}
coupling_q_lookup <- function(data, z, q) data %>% filter(contrast == z, stratum == q) %>% pull(rho)
mult_rows_text <- function(z) {
  x <- coupling_multiplicity %>% filter(contrast == z)
  paste0(x$stratum, ": ", fmt_num(x$rho), " (n=", x$n_genes, ")", collapse = "; ")
}

strongest_sentence <- switch(
  FINAL_LEVEL,
  "LEVEL A" = "在高内部计数信度和预设质量控制条件下，matched gene-abundance change 与 long-read isoform-composition remodeling 仅呈有限耦合，支持二者作为相互补充而不可互换的分子层面。",
  "LEVEL B" = "Matched gene-abundance change 与 long-read isoform-composition remodeling 部分重叠但并不等价，其耦合程度受到测量信度和基因特征影响。",
  "LEVEL C" = "Measurement-aware sensitivity analyses indicate that the apparent weak cross-layer coupling is materially influenced by measurement depth or complexity and therefore requires cautious interpretation."
)
centerpiece_decision <- if (FINAL_LEVEL == "LEVEL A") {
  "支持继续作为 conceptual centerpiece，但主张应限定为 limited/incomplete coupling，而不是生物学独立。"
} else if (FINAL_LEVEL == "LEVEL B") {
  "可保留为 conceptual component，但应降格为 partially overlapping, non-equivalent layers，不能作为强独立性主张。"
} else {
  "不支持继续把 weak coupling 作为全文 conceptual centerpiece；应转为测量限制和敏感性结果。"
}

report_lines <- c(
  "# Figure 3 measurement robustness report",
  "",
  "## 执行边界",
  "",
  "- 本轮仅使用冻结 primary72、1,503-gene common universe、已有 long-read raw counts 和 matched short-read estimated counts进行派生稳健性分析。",
  "- 未重跑 FASTQ、Kallisto、FLAIR、SQANTI3、DTE 或 primary DTU，未修改正式 Figure 3，也未修改 Figure 1/2/4-8。",
  "- Long-read split 是互补 Binomial count-resampling pseudo-replicates；short-read split 是等深度 multinomial pseudo-libraries。二者均不是技术重复。",
  "- Reliability-adjusted rho 仅是 approximate sensitivity estimate，不是‘真实相关’。",
  "",
  "## 预设分析定义",
  "",
  "- common universe：T vs A 和 T vs N 各1,503个基因，两个集合完全一致。",
  "- TVD：对每个基因全部5,477个正式 eligible isoforms 中属于该基因的转录本计算；由原始计数重构的 TVD 与冻结表最大误差小于 1e-12。",
  "- long-read abundance：全部映射到该基因的 long-read isoform counts 在72个样本中的基因总count均值。",
  "- zero_fraction：正式 eligible isoform-by-sample cells 中的零计数比例。",
  paste0("- n_effective_isoforms：在至少", EFFECTIVE_SUPPORT_SAMPLES, "/72个样本中count>0的 eligible isoform 数；该阈值在结果计算前固定。"),
  "- high-confidence subset：LR mean gene count不低于共同universe中位数、eligible isoform multiplicity <=5、zero_fraction不高于共同universe中位数。",
  "",
  "## 1. Long-read TVD pseudo-replicate reliability",
  "",
  paste0("- T vs A：median=", fmt_num(lookup_one(lr_reliability_summary, "T_vs_A", "median")),
         "，IQR ", fmt_num(lookup_one(lr_reliability_summary, "T_vs_A", "Q1")), "-", fmt_num(lookup_one(lr_reliability_summary, "T_vs_A", "Q3")),
         "，2.5-97.5 percentile ", fmt_num(lookup_one(lr_reliability_summary, "T_vs_A", "percentile_2.5")), "-", fmt_num(lookup_one(lr_reliability_summary, "T_vs_A", "percentile_97.5")), "。"),
  paste0("- T vs N：median=", fmt_num(lookup_one(lr_reliability_summary, "T_vs_N", "median")),
         "，IQR ", fmt_num(lookup_one(lr_reliability_summary, "T_vs_N", "Q1")), "-", fmt_num(lookup_one(lr_reliability_summary, "T_vs_N", "Q3")),
         "，2.5-97.5 percentile ", fmt_num(lookup_one(lr_reliability_summary, "T_vs_N", "percentile_2.5")), "-", fmt_num(lookup_one(lr_reliability_summary, "T_vs_N", "percentile_97.5")), "。"),
  paste0("- abundance dependence（Q1 -> Q4）：T vs A ", fmt_num(lr_q_lookup("T_vs_A", "Q1")), " -> ", fmt_num(lr_q_lookup("T_vs_A", "Q4")),
         "；T vs N ", fmt_num(lr_q_lookup("T_vs_N", "Q1")), " -> ", fmt_num(lr_q_lookup("T_vs_N", "Q4")), "。"),
  "",
  "## 2. Short-read |log2FC| pseudo-replicate reliability",
  "",
  paste0("- T vs A：median=", fmt_num(lookup_one(sr_reliability_summary, "T_vs_A", "median")),
         "，2.5-97.5 percentile ", fmt_num(lookup_one(sr_reliability_summary, "T_vs_A", "percentile_2.5")), "-", fmt_num(lookup_one(sr_reliability_summary, "T_vs_A", "percentile_97.5")), "。"),
  paste0("- T vs N：median=", fmt_num(lookup_one(sr_reliability_summary, "T_vs_N", "median")),
         "，2.5-97.5 percentile ", fmt_num(lookup_one(sr_reliability_summary, "T_vs_N", "percentile_2.5")), "-", fmt_num(lookup_one(sr_reliability_summary, "T_vs_N", "percentile_97.5")), "。"),
  "- Fractional tximport estimated counts 未直接使用 rbinom；每个样本在1,503-gene universe内归一化为概率，并以近似等深度 rmultinom pseudo-libraries 重新运行 paired DESeq2。",
  "",
  "## 3. Observed and reliability-adjusted coupling",
  "",
  paste0("- T vs A：observed rho=", fmt_num(lookup_one(reliability_adjusted, "T_vs_A", "rho_observed")),
         "；approximate adjusted rho=", fmt_num(lookup_one(reliability_adjusted, "T_vs_A", "rho_reliability_adjusted")),
         "（95% CI ", fmt_num(lookup_one(reliability_adjusted, "T_vs_A", "adjusted_CI_low")), "-", fmt_num(lookup_one(reliability_adjusted, "T_vs_A", "adjusted_CI_high")), "）。"),
  paste0("- T vs N：observed rho=", fmt_num(lookup_one(reliability_adjusted, "T_vs_N", "rho_observed")),
         "；approximate adjusted rho=", fmt_num(lookup_one(reliability_adjusted, "T_vs_N", "rho_reliability_adjusted")),
         "（95% CI ", fmt_num(lookup_one(reliability_adjusted, "T_vs_N", "adjusted_CI_low")), "-", fmt_num(lookup_one(reliability_adjusted, "T_vs_N", "adjusted_CI_high")), "）。"),
  "",
  "## 4. Expression-depth and multiplicity stratification",
  "",
  paste0("- LR abundance quartiles T vs A Q1-Q4 rho：", paste(fmt_num(coupling_lr_abundance$rho[coupling_lr_abundance$contrast == "T_vs_A"]), collapse = ", "), "。"),
  paste0("- LR abundance quartiles T vs N Q1-Q4 rho：", paste(fmt_num(coupling_lr_abundance$rho[coupling_lr_abundance$contrast == "T_vs_N"]), collapse = ", "), "。"),
  paste0("- SR baseMean quartiles T vs A Q1-Q4 rho：", paste(fmt_num(coupling_sr_abundance$rho[coupling_sr_abundance$contrast == "T_vs_A"]), collapse = ", "), "。"),
  paste0("- SR baseMean quartiles T vs N Q1-Q4 rho：", paste(fmt_num(coupling_sr_abundance$rho[coupling_sr_abundance$contrast == "T_vs_N"]), collapse = ", "), "。"),
  paste0("- Multiplicity T vs A：", mult_rows_text("T_vs_A"), "。"),
  paste0("- Multiplicity T vs N：", mult_rows_text("T_vs_N"), "。"),
  paste0("- 预设 major high-quality strata 的最大 |rho|=", fmt_num(max_high_quality_rho), "；LEVEL阈值为0.50。"),
  "",
  "## 5. Partial correlation",
  "",
  paste0("- Primary adjustment T vs A：partial rho=", fmt_num(lookup_one(primary_partial, "T_vs_A", "partial_rho")),
         "（95% CI ", fmt_num(lookup_one(primary_partial, "T_vs_A", "CI_low")), "-", fmt_num(lookup_one(primary_partial, "T_vs_A", "CI_high")), "）。"),
  paste0("- Primary adjustment T vs N：partial rho=", fmt_num(lookup_one(primary_partial, "T_vs_N", "partial_rho")),
         "（95% CI ", fmt_num(lookup_one(primary_partial, "T_vs_N", "CI_low")), "-", fmt_num(lookup_one(primary_partial, "T_vs_N", "CI_high")), "）。"),
  "- Primary model 对 X/Y及4个covariates作normal-score rank转换后，分别残差化并计算rank-residual correlation；sensitivity model仅调整LR count和isoform multiplicity。",
  "",
  "## 6. Incremental explanatory value",
  "",
  paste0("- T vs A rank-linear：baseline CV explained variance=", fmt_num(lookup_one(primary_incremental, "T_vs_A", "baseline_CV_median")),
         "，extended=", fmt_num(lookup_one(primary_incremental, "T_vs_A", "extended_CV_median")),
         "，delta=", fmt_num(lookup_one(primary_incremental, "T_vs_A", "delta_CV_median")), "。"),
  paste0("- T vs N rank-linear：baseline CV explained variance=", fmt_num(lookup_one(primary_incremental, "T_vs_N", "baseline_CV_median")),
         "，extended=", fmt_num(lookup_one(primary_incremental, "T_vs_N", "extended_CV_median")),
         "，delta=", fmt_num(lookup_one(primary_incremental, "T_vs_N", "delta_CV_median")), "。"),
  "- Quasibinomial sensitivity 使用边界压缩后的TVD并独立重复100次10-fold CV；它用于检查模型形式，不替代rank-linear主结果。",
  "",
  "## 7. High-confidence subset",
  "",
  paste0("- T vs A：n=", lookup_one(high_confidence_summary, "T_vs_A", "n_high_confidence_genes"),
         "，observed rho=", fmt_num(lookup_one(high_confidence_summary, "T_vs_A", "observed_rho")),
         "，partial rho=", fmt_num(lookup_one(high_confidence_summary, "T_vs_A", "partial_rho")),
         "，delta CV=", fmt_num(lookup_one(high_confidence_summary, "T_vs_A", "incremental_delta_CV_median")), "。"),
  paste0("- T vs N：n=", lookup_one(high_confidence_summary, "T_vs_N", "n_high_confidence_genes"),
         "，observed rho=", fmt_num(lookup_one(high_confidence_summary, "T_vs_N", "observed_rho")),
         "，partial rho=", fmt_num(lookup_one(high_confidence_summary, "T_vs_N", "partial_rho")),
         "，delta CV=", fmt_num(lookup_one(high_confidence_summary, "T_vs_N", "incremental_delta_CV_median")), "。"),
  "",
  "## 8. Max-absolute-delta sensitivity",
  "",
  paste0("- Raw rho T vs A/T vs N：", paste(fmt_num(maxdelta_raw_check$estimate), collapse = "/"), "。"),
  paste0("- Primary partial rho T vs A/T vs N：", paste(fmt_num(maxdelta_partial_primary$estimate), collapse = "/"), "。"),
  paste0("- 与TVD主分析方向是否一致：", ifelse(MAXDELTA_DIRECTION_CONSISTENT, "是", "否"), "。"),
  "",
  "## 9. Threshold overlap as a secondary result",
  "",
  "- 保留冻结结果：T vs A中20/99 strict DTU genes为matched short-read DEG，Fisher OR=2.66；T vs N为16/74，OR=3.21。",
  "- 固定逻辑顺序：DTU genes在matched short-read DEGs中富集；集合重叠仍不完整；连续效应量的耦合有限。",
  "- 79.8%/78.4%不能单独写为short-read misses most DTU，也不是short-read error rate。",
  "",
  "## 10. 预设LEVEL判定",
  "",
  paste0("- 最终等级：**", FINAL_LEVEL, "**。"),
  paste0("- 解释：", level_interpretation, "。"),
  paste0("- Figure 3地位：", centerpiece_decision),
  "",
  "## 11. 可写入摘要的最强句子",
  "",
  paste0("> ", strongest_sentence),
  "",
  "## 12. 禁止的过度结论",
  "",
  "- 不得将pseudo-replicates称为technical replicates。",
  "- 不得将reliability-adjusted rho称为真实相关或去噪后的真实效应。",
  "- 不得把低rho自动解释为两个生物层面相互独立。",
  "- 不得把79.8%/78.4%解释为short-read错误率或完整full-length isoform漏检率。",
  "- 不得把matched short-read数据写成独立外部验证队列。",
  "- 不得根据分层结果重新定义1,503-gene universe、DTU/DEG阈值或high-confidence cutoff。",
  "- 不得把统计耦合、增量解释度或可靠性校正写成功能机制、因果关系或测序平台优越性。",
  "",
  "## 13. 图件定位",
  "",
  "- `figures/Figure3_robustness_candidate.*`仅供内部审阅，不自动替换正式Figure 3。",
  "- 四套Supplementary robustness figures分别展示reliability、abundance strata、multiplicity strata和incremental CV。",
  paste0("- 覆盖前候选输出归档：", ifelse(is.na(archive_path), "none", archive_path), "。")
)
writeLines(report_lines, file.path(DIR_REPORTS, "FIG3_MEASUREMENT_ROBUSTNESS_REPORT.md"), useBytes = TRUE)

brief_lines <- c(
  "# ChatGPT review brief: Figure 3 measurement robustness",
  "",
  paste0("- **Pre-specified decision:** ", FINAL_LEVEL, " — ", level_interpretation, "."),
  "- **Frozen universe:** 1,503 multi-isoform genes per T vs A/T vs N; frozen observed rho reproduced exactly.",
  paste0("- **LR TVD reliability:** ", paste(paste0(contrast_label(lr_reliability_summary$contrast), "=", fmt_num(lr_reliability_summary$median)), collapse = "; "), "."),
  paste0("- **SR |log2FC| reliability:** ", paste(paste0(contrast_label(sr_reliability_summary$contrast), "=", fmt_num(sr_reliability_summary$median)), collapse = "; "), "."),
  paste0("- **Observed -> adjusted rho:** T vs A ", fmt_num(decision_values$observed_rho[[1]]), " -> ", fmt_num(decision_values$adjusted_rho[[1]]),
         "; T vs N ", fmt_num(decision_values$observed_rho[[2]]), " -> ", fmt_num(decision_values$adjusted_rho[[2]]), "."),
  paste0("- **Primary partial rho:** ", paste(paste0(contrast_label(decision_values$contrast), "=", fmt_num(decision_values$partial_rho)), collapse = "; "), "."),
  paste0("- **High-confidence rho:** ", paste(paste0(contrast_label(decision_values$contrast), "=", fmt_num(decision_values$high_confidence_rho)), collapse = "; "), "."),
  paste0("- **Incremental rank-linear delta CV:** ", paste(paste0(contrast_label(decision_values$contrast), "=", fmt_num(decision_values$delta_CV)), collapse = "; "), "."),
  paste0("- **Highest pre-defined high-quality stratum |rho|:** ", fmt_num(max_high_quality_rho), "."),
  paste0("- **Max-delta direction consistency:** ", MAXDELTA_DIRECTION_CONSISTENT, "."),
  paste0("- **Conceptual centerpiece decision:** ", centerpiece_decision),
  paste0("- **Strongest abstract-safe sentence:** ", strongest_sentence),
  "- **Boundary:** pseudo-replicates are count-resampling constructs; attenuation adjustment is approximate; matched short read is not external validation.",
  "- **Scope:** no upstream quantification/model rerun and no formal Figure 3 replacement."
)
writeLines(brief_lines, file.path(DIR_REPORTS, "CHATGPT_REVIEW_BRIEF_FIG3_ROBUSTNESS.md"), useBytes = TRUE)

panel_source_paths <- c(
  panel_A = write_tsv2(reliability_plot_data, file.path(DIR_RESULTS, "Fig3_robustness_candidate_panel_A.tsv")),
  panel_B = write_tsv2(coupling_plot_data, file.path(DIR_RESULTS, "Fig3_robustness_candidate_panel_B.tsv")),
  panel_C = write_tsv2(incremental_plot_data, file.path(DIR_RESULTS, "Fig3_robustness_candidate_panel_C.tsv")),
  panel_D = write_tsv2(strata_plot_data, file.path(DIR_RESULTS, "Fig3_robustness_candidate_panel_D.tsv"))
)

result_paths <- c(file.path(DIR_RESULTS, c(
  "Fig3_measurement_covariates.tsv",
  "Fig3_LR_TVD_reliability_all_replicates.tsv",
  "Fig3_LR_TVD_reliability_summary.tsv",
  "Fig3_LR_TVD_reliability_stratified.tsv",
  "Fig3_SR_log2FC_reliability_all_replicates.tsv",
  "Fig3_SR_log2FC_reliability_summary.tsv",
  "Fig3_reliability_adjusted_coupling.tsv",
  "Fig3_coupling_by_LR_abundance.tsv",
  "Fig3_coupling_by_SR_abundance.tsv",
  "Fig3_coupling_by_isoform_multiplicity.tsv",
  "Fig3_partial_correlation.tsv",
  "Fig3_incremental_explanatory_value.tsv",
  "Fig3_incremental_CV_all_replicates.tsv",
  "Fig3_high_confidence_subset.tsv",
  "Fig3_high_confidence_summary.tsv",
  "Fig3_maxdelta_sensitivity.tsv",
  "Fig3_robustness_key_statistics.tsv",
  "Fig3_robustness_level_decision.tsv",
  "Fig3_robustness_decision_values.tsv",
  "Fig3_SR_pseudolibrary_depth_audit.tsv",
  "Fig3_robustness_input_manifest.tsv",
  "Fig3_robustness_sample_audit.tsv"
)), unname(panel_source_paths))
figure_paths <- c(
  candidate_paths$path, supp_reliability_paths$path, supp_abundance_paths$path,
  supp_multiplicity_paths$path, supp_incremental_paths$path,
  file.path(DIR_FIGURES, "Fig3_incremental_explanatory_value.pdf")
)
stop_if(any(!file.exists(result_paths) | file.info(result_paths)$size <= 0), "A required robustness result table is missing or empty.")
stop_if(any(!file.exists(figure_paths) | file.info(figure_paths)$size <= 0), "A required robustness figure is missing or empty.")

source_index <- tibble(
  figure = c(
    rep("Figure3 robustness candidate", 4),
    "Supplementary reliability", "Supplementary abundance", "Supplementary multiplicity", "Supplementary incremental CV"
  ),
  panel = c("A", "B", "C", "D", "all", "all", "all", "all"),
  source_data_path = c(
    panel_source_paths[["panel_A"]],
    panel_source_paths[["panel_B"]],
    panel_source_paths[["panel_C"]],
    panel_source_paths[["panel_D"]],
    file.path(DIR_RESULTS, "Fig3_LR_TVD_reliability_all_replicates.tsv"),
    file.path(DIR_RESULTS, "Fig3_coupling_by_SR_abundance.tsv"),
    file.path(DIR_RESULTS, "Fig3_coupling_by_isoform_multiplicity.tsv"),
    file.path(DIR_RESULTS, "Fig3_incremental_CV_all_replicates.tsv")
  ),
  n_rows = NA_integer_,
  statistical_unit = c(
    "pseudo-split replicate summary", "gene plus reliability distribution", "CV repeat", "gene stratum",
    "long-read pseudo-split replicate", "gene stratum", "gene stratum", "CV repeat"
  ),
  description = c(
    "LR and SR pseudo-replicate reliability", "Observed, approximate adjusted and partial coupling",
    "Baseline versus extended incremental explanatory value", "LR abundance and multiplicity stratified rho",
    "All long-read reliability replicates", "LR/SR abundance quartile coupling",
    "Multiplicity-stratified coupling", "Repeated held-out explained variance"
  )
) %>%
  mutate(n_rows = vapply(source_data_path, count_tsv_rows, integer(1)))
source_index_path <- write_tsv2(source_index, file.path(DIR_RESULTS, "Fig3_robustness_source_data_index.tsv"))

record_progress("QA", "Running output, field-count, source-index, boundary and scope guards")
all_tsv_paths <- c(result_paths, source_index_path)
field_qa <- tsv_field_audit(all_tsv_paths)
stop_if(any(!field_qa$field_count_pass), "One or more robustness TSV files have inconsistent field counts.")
field_qa_path <- write_tsv2(field_qa, file.path(DIR_RESULTS, "Fig3_robustness_TSV_field_QA.tsv"))

source_index_check <- source_index %>%
  mutate(
    exists_nonempty = file.exists(source_data_path) & file.info(source_data_path)$size > 0,
    observed_rows = vapply(source_data_path, count_tsv_rows, integer(1)),
    row_count_match = observed_rows == n_rows
  )
stop_if(any(!source_index_check$exists_nonempty | !source_index_check$row_count_match), "Robustness source-data index QA failed.")

svg_paths <- figure_paths[grepl("[.]svg$", figure_paths)]
forbidden_terms <- c(
  "technical replicate", "true correlation", "independent external validation",
  "validated full-length isoform", "biological independence", "candidate prioritization",
  "clinical association", "confirmed mechanism"
)
boundary_grep <- bind_rows(lapply(svg_paths, function(path) {
  txt <- tolower(paste(readLines(path, warn = FALSE), collapse = "\n"))
  tibble(path = path, term = forbidden_terms,
         hit = vapply(forbidden_terms, function(term) grepl(tolower(term), txt, fixed = TRUE), logical(1)))
}))
stop_if(any(boundary_grep$hit), "A robustness SVG contains a prohibited claim-boundary phrase.")
boundary_path <- write_tsv2(boundary_grep, file.path(DIR_RESULTS, "Fig3_robustness_boundary_grep.tsv"))

guard_after <- file_snapshot(guard_paths)
guard_check <- guard_before %>%
  rename(size_before = size_bytes, mtime_before = mtime, sha256_before = sha256) %>%
  left_join(
    guard_after %>% rename(size_after = size_bytes, mtime_after = mtime, sha256_after = sha256),
    by = "path"
  ) %>%
  mutate(unchanged = !is.na(sha256_after) & size_before == size_after & sha256_before == sha256_after)
stop_if(any(!guard_check$unchanged), "Figure 1/2/formal Figure 3/Figure 4-8 scope guard failed.")
guard_path <- write_tsv2(guard_check, file.path(DIR_RESULTS, "Fig3_robustness_unrelated_figure_guard.tsv"))

post_input_hashes <- vapply(INPUT_PATHS, sha256_file, character(1))
stop_if(!identical(unname(post_input_hashes), unname(observed_hashes)), "A frozen input changed during robustness analysis.")

report_paths <- file.path(DIR_REPORTS, c(
  "FIG3_MEASUREMENT_ROBUSTNESS_REPORT.md",
  "CHATGPT_REVIEW_BRIEF_FIG3_ROBUSTNESS.md"
))
output_manifest <- tibble(
  output_path = unique(c(result_paths, source_index_path, field_qa_path, boundary_path, guard_path, figure_paths, report_paths))
) %>%
  mutate(
    exists = file.exists(output_path), size_bytes = file.info(output_path)$size,
    sha256 = vapply(output_path, sha256_file, character(1))
  )
stop_if(any(!output_manifest$exists | output_manifest$size_bytes <= 0), "Output manifest contains missing or empty files.")
write_tsv2(output_manifest, file.path(DIR_RESULTS, "Fig3_robustness_output_manifest.tsv"))

qa_lines <- c(
  "# Figure 3 measurement robustness programmatic QA",
  "",
  "- Frozen input SHA256 gates: PASS.",
  "- primary72 gate (24 patients / 72 matched N/A/T samples): PASS.",
  "- Common universe gate (1,503 genes per contrast): PASS.",
  "- Frozen observed rho reproduction (0.13303077 / 0.10173881): PASS.",
  "- TVD reconstruction from 5,477 eligible isoforms (maximum error <1e-12): PASS.",
  "- LR complementary binomial pseudo-splits: 500/500 per requested analysis.",
  "- SR multinomial pseudo-library splits plus paired DESeq2: 200/200; 1,503/1,503 finite genes in every split.",
  "- Reliability-adjusted bootstrap: 4,000 gene resamples combined with reliability distributions.",
  "- Abundance/multiplicity and partial-correlation bootstrap: 4,000 gene resamples per requested analysis.",
  "- Incremental models: rank-linear plus quasibinomial sensitivity; 100 repeated 10-fold CV.",
  "- High-confidence subset cutoffs applied before outcome inspection: PASS.",
  "- TSV field counts and source-data row counts: PASS.",
  "- PDF/SVG/PNG existence and non-empty checks: PASS.",
  "- SVG claim-boundary grep: PASS.",
  paste0("- Figure 1/2/formal Figure 3/Figure 4-8 guard: ", sum(guard_check$unchanged), "/", nrow(guard_check), " unchanged."),
  "- Final local visual QA remains required after PNG download."
)
writeLines(qa_lines, file.path(DIR_REPORTS, "FIG3_MEASUREMENT_ROBUSTNESS_QA.md"), useBytes = TRUE)

record_progress("COMPLETE", paste0("Figure 3 measurement robustness completed: ", FINAL_LEVEL))
cat("Figure 3 measurement robustness completed.\n")
cat("Decision:", FINAL_LEVEL, "-", level_interpretation, "\n")
cat("Report:", file.path(DIR_REPORTS, "FIG3_MEASUREMENT_ROBUSTNESS_REPORT.md"), "\n")
cat("Candidate figure:", file.path(DIR_FIGURES, "Figure3_robustness_candidate.pdf"), "\n")
cat("Formal Figure 3 and all upstream analyses were not modified.\n")
