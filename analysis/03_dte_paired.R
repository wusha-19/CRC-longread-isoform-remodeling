#!/usr/bin/env Rscript
set.seed(20260624)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

root <- Sys.getenv("CRC_ROOT")
if (!nzchar(root)) stop("Set CRC_ROOT before running this script.")
phase01 <- file.path(root, "05_isoform_main_study/01_master_catalog")
out <- file.path(root, "05_isoform_main_study/02_DTE_paired")
dirs <- c("code", "config", "tables", "figures", "reports", "logs", "provenance", "matrices")
for (d in dirs) dir.create(file.path(out, d), recursive=TRUE, showWarnings=FALSE)

counts_primary_path <- file.path(phase01, "matrices/counts_primary72_isoforms.tsv")
meta_primary_path <- file.path(phase01, "matrices/metadata_primary72.tsv")
counts_sens_path <- file.path(phase01, "matrices/counts_sensitivity75_isoforms.tsv")
meta_sens_path <- file.path(phase01, "matrices/metadata_sensitivity75.tsv")
catalog_path <- file.path(phase01, "tables/CRC_isoform_master_catalog.tsv")

contrast_defs <- list(
  A_vs_N = c("tissue_group", "A", "N"),
  T_vs_A = c("tissue_group", "T", "A"),
  T_vs_N = c("tissue_group", "T", "N")
)

read_counts <- function(path) {
  x <- read.delim(path, check.names=FALSE)
  rownames(x) <- x[[1]]
  x <- x[, -1, drop=FALSE]
  x <- as.matrix(x)
  storage.mode(x) <- "integer"
  x
}

robust_z <- function(x) {
  med <- median(x, na.rm=TRUE)
  madv <- median(abs(x - med), na.rm=TRUE)
  if (is.na(madv) || madv == 0) {
    s <- sd(x, na.rm=TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    return((x - mean(x, na.rm=TRUE)) / s)
  }
  0.6745 * (x - med) / madv
}

check_inputs <- function(counts, meta, cohort) {
  stopifnot("sample_id" %in% names(meta))
  if (!identical(colnames(counts), meta$sample_id)) {
    stop(sprintf("%s sample order mismatch between count matrix and metadata", cohort))
  }
  if (any(meta$tissue_group == "D")) stop(sprintf("%s contains D samples", cohort))
  if (cohort == "primary72" && any(meta$CRC20_status == "CRC20")) stop("primary72 contains CRC_20")
  if (cohort == "sensitivity75" && !any(meta$CRC20_status == "CRC20")) stop("sensitivity75 lacks CRC_20")
  tab <- table(meta$patient_id, meta$tissue_group)
  if (!all(colnames(tab) %in% c("N","A","T"))) stop(sprintf("%s has unexpected tissue group", cohort))
  complete <- apply(tab[, c("N","A","T"), drop=FALSE], 1, function(z) all(z == 1))
  if (!all(complete)) stop(sprintf("%s has incomplete N/A/T patient set", cohort))
}

save_pdf <- function(path, width, height, expr) {
  pdf(path, width=width, height=height)
  force(expr)
  dev.off()
}

plot_volcano <- function(df, title, path) {
  dat <- data.frame(log2FoldChange=df$log2FoldChange, padj=df$padj, baseMean=df$baseMean)
  dat$neglog10padj <- -log10(pmax(dat$padj, .Machine$double.xmin))
  dat$sig_lfc1 <- !is.na(dat$padj) & dat$padj < 0.05 & abs(dat$log2FoldChange) >= 1 & dat$baseMean >= 10
  p <- ggplot(dat, aes(log2FoldChange, neglog10padj)) +
    geom_point(aes(color=sig_lfc1), alpha=0.45, size=0.35) +
    scale_color_manual(values=c("FALSE"="grey70", "TRUE"="#B23A48")) +
    theme_bw(base_size=10) + labs(title=title, x="log2 fold-change", y="-log10 adjusted P", color="sig")
  ggsave(path, p, width=6.5, height=5)
  write.table(dat, sub("\\.pdf$", ".data.tsv", path), sep="\t", quote=FALSE, row.names=FALSE)
}

plot_ma <- function(df, title, path) {
  dat <- data.frame(baseMean=df$baseMean, log2FoldChange=df$log2FoldChange, padj=df$padj)
  dat$sig_lfc1 <- !is.na(dat$padj) & dat$padj < 0.05 & abs(dat$log2FoldChange) >= 1 & dat$baseMean >= 10
  p <- ggplot(dat, aes(log10(baseMean + 1), log2FoldChange)) +
    geom_point(aes(color=sig_lfc1), alpha=0.45, size=0.35) +
    scale_color_manual(values=c("FALSE"="grey70", "TRUE"="#B23A48")) +
    theme_bw(base_size=10) + labs(title=title, x="log10(baseMean + 1)", y="log2 fold-change", color="sig")
  ggsave(path, p, width=6.5, height=5)
  write.table(dat, sub("\\.pdf$", ".data.tsv", path), sep="\t", quote=FALSE, row.names=FALSE)
}

patient_consistency <- function(res_df, logmat, meta, contrast, sig_ids) {
  pieces <- strsplit(contrast, "_vs_")[[1]]
  g1 <- pieces[1]; g0 <- pieces[2]
  rows <- list()
  for (id in sig_ids) {
    overall <- res_df$log2FoldChange[match(id, res_df$isoform_id)]
    if (is.na(overall) || overall == 0) next
    same <- 0; evaln <- 0
    for (pid in unique(meta$patient_id)) {
      s1 <- meta$sample_id[meta$patient_id == pid & meta$tissue_group == g1]
      s0 <- meta$sample_id[meta$patient_id == pid & meta$tissue_group == g0]
      if (length(s1) == 1 && length(s0) == 1 && id %in% rownames(logmat)) {
        diff <- logmat[id, s1] - logmat[id, s0]
        if (!is.na(diff) && diff != 0) {
          evaln <- evaln + 1
          same <- same + as.integer(sign(diff) == sign(overall))
        }
      }
    }
    rows[[length(rows)+1]] <- data.frame(
      isoform_id=id, contrast=contrast, overall_log2FC=overall,
      n_patients_evaluable=evaln, n_patients_same_direction=same,
      patient_direction_consistency=ifelse(evaln > 0, same/evaln, NA_real_)
    )
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

run_cohort <- function(cohort, counts_path, meta_path, catalog) {
  counts <- read_counts(counts_path)
  meta <- read.delim(meta_path, check.names=FALSE)
  check_inputs(counts, meta, cohort)
  meta$tissue_group <- factor(meta$tissue_group, levels=c("N","A","T"))
  meta$patient_id <- factor(meta$patient_id)
  rownames(meta) <- meta$sample_id

  dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=~ patient_id + tissue_group)
  dds <- estimateSizeFactors(dds, type="poscounts")
  sf <- sizeFactors(dds)
  dds <- DESeq(dds, quiet=FALSE)
  vsd <- vst(dds, blind=FALSE)
  vst_mat <- assay(vsd)

  lib <- colSums(counts)
  detected <- colSums(counts > 0)
  qc <- data.frame(
    sample_id=colnames(counts),
    analysis_cohort=cohort,
    library_size=as.numeric(lib),
    detected_isoforms=as.numeric(detected),
    size_factor=as.numeric(sf),
    patient_id=meta$patient_id,
    tissue_group=meta$tissue_group,
    CRC20_status=meta$CRC20_status
  )
  qc$library_size_robust_z <- robust_z(qc$library_size)
  qc$detected_isoforms_robust_z <- robust_z(qc$detected_isoforms)

  vars <- apply(vst_mat, 1, var)
  top <- names(sort(vars, decreasing=TRUE))[seq_len(min(5000, length(vars)))]
  pca <- prcomp(t(vst_mat[top,]), center=TRUE, scale.=TRUE)
  pca_df <- data.frame(sample_id=rownames(pca$x), PC1=pca$x[,1], PC2=pca$x[,2],
                       PC3=pca$x[,3], PC4=pca$x[,4], analysis_cohort=cohort)
  pca_df$PC1_var <- summary(pca)$importance[2,1]
  pca_df$PC2_var <- summary(pca)$importance[2,2]
  pca_df <- merge(pca_df, meta, by.x="sample_id", by.y="row.names", sort=FALSE)
  centers <- aggregate(cbind(PC1,PC2,PC3,PC4) ~ tissue_group, data=pca_df, mean)
  dist <- rep(NA_real_, nrow(pca_df))
  for (i in seq_len(nrow(pca_df))) {
    c0 <- centers[centers$tissue_group == pca_df$tissue_group[i], c("PC1","PC2","PC3","PC4")]
    dist[i] <- sqrt(sum((as.numeric(pca_df[i,c("PC1","PC2","PC3","PC4")]) - as.numeric(c0))^2))
  }
  pca_df$distance_to_tissue_center_PC1_4 <- dist
  pca_df$distance_to_tissue_center_robust_z <- ave(pca_df$distance_to_tissue_center_PC1_4, pca_df$tissue_group, FUN=robust_z)
  qc <- merge(qc, pca_df[,c("sample_id","distance_to_tissue_center_PC1_4","distance_to_tissue_center_robust_z")], by="sample_id", sort=FALSE)
  qc$new_qc_flag <- abs(qc$library_size_robust_z) >= 3 | abs(qc$detected_isoforms_robust_z) >= 3 | abs(qc$distance_to_tissue_center_robust_z) >= 3

  write.table(qc, file.path(out, "tables", paste0(cohort, "_DTE_preQC_sample_metrics.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(pca_df, file.path(out, "tables", paste0(cohort, "_DTE_PCA_VST_coordinates.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  corr <- cor(vst_mat[top,], method="pearson")
  write.table(corr, file.path(out, "tables", paste0(cohort, "_DTE_sample_correlation.tsv")), sep="\t", quote=FALSE)

  p <- ggplot(pca_df, aes(PC1, PC2, color=tissue_group, label=sample_id)) +
    geom_point(size=2) + theme_bw(base_size=10) +
    labs(title=paste0(cohort, " VST PCA"), x=sprintf("PC1 (%.1f%%)", 100*pca_df$PC1_var[1]), y=sprintf("PC2 (%.1f%%)", 100*pca_df$PC2_var[1]))
  ggsave(file.path(out, "figures", paste0(cohort, "_DTE_PCA_VST.pdf")), p, width=6.5, height=5.5)
  p2 <- ggplot(qc, aes(tissue_group, library_size, color=tissue_group)) + geom_jitter(width=0.08, height=0, size=2) + theme_bw(base_size=10) + labs(title=paste0(cohort, " library size"))
  ggsave(file.path(out, "figures", paste0(cohort, "_DTE_library_size.pdf")), p2, width=6, height=4.5)
  corr_df <- as.data.frame(as.table(corr)); names(corr_df) <- c("sample1","sample2","correlation")
  p3 <- ggplot(corr_df, aes(sample1, sample2, fill=correlation)) + geom_tile() + theme_bw(base_size=6) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0.8) + labs(title=paste0(cohort, " sample correlation"), x="", y="")
  ggsave(file.path(out, "figures", paste0(cohort, "_DTE_sample_correlation.pdf")), p3, width=8.5, height=8)
  write.table(corr_df, file.path(out, "figures", paste0(cohort, "_DTE_sample_correlation.data.tsv")), sep="\t", quote=FALSE, row.names=FALSE)

  all_cons <- list()
  result_list <- list()
  for (nm in names(contrast_defs)) {
    res <- as.data.frame(results(dds, contrast=contrast_defs[[nm]], alpha=0.05))
    res$isoform_id <- rownames(res)
    res$contrast <- nm
    res$analysis_cohort <- cohort
    res <- merge(res, catalog[,c("isoform_id","gene_id","gene_name","sqanti3_category")], by="isoform_id", all.x=TRUE, sort=FALSE)
    res <- res[,c("isoform_id","gene_id","gene_name","sqanti3_category","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","contrast","analysis_cohort")]
    result_list[[nm]] <- res
    write.table(res, file.path(out, "tables", paste0(cohort, "_DTE_", nm, "_all.tsv")), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    sig1 <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1 & res$baseMean >= 10,]
    sig058 <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 0.58 & res$baseMean >= 10,]
    write.table(sig1, file.path(out, "tables", paste0(cohort, "_DTE_", nm, "_sig_lfc1.tsv")), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    write.table(sig058, file.path(out, "tables", paste0(cohort, "_DTE_", nm, "_sig_lfc0.58.tsv")), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    plot_volcano(res, paste0(cohort, " ", nm), file.path(out, "figures", paste0(cohort, "_volcano_", nm, ".pdf")))
    plot_ma(res, paste0(cohort, " ", nm), file.path(out, "figures", paste0(cohort, "_MA_", nm, ".pdf")))
    cons <- patient_consistency(res, log2(counts(dds, normalized=TRUE) + 1), meta, nm, sig1$isoform_id)
    if (nrow(cons) > 0) all_cons[[nm]] <- cons
  }
  cons_df <- if (length(all_cons)) do.call(rbind, all_cons) else data.frame()
  write.table(cons_df, file.path(out, "tables", paste0(cohort, "_DTE_patient_direction_consistency.tsv")), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  list(dds=dds, qc=qc, results=result_list, consistency=cons_df)
}

catalog <- read.delim(catalog_path, check.names=FALSE)
primary <- run_cohort("primary72", counts_primary_path, meta_primary_path, catalog)
sensitivity <- run_cohort("sensitivity75", counts_sens_path, meta_sens_path, catalog)

summary_rows <- list()
cat_rows <- list()
for (cohort in c("primary72","sensitivity75")) {
  reslist <- if (cohort == "primary72") primary$results else sensitivity$results
  for (nm in names(reslist)) {
    res <- reslist[[nm]]
    sig1 <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1 & res$baseMean >= 10
    sig058 <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 0.58 & res$baseMean >= 10
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      analysis_cohort=cohort, contrast=nm, n_tested_isoforms=nrow(res),
      n_padj_lt_0.05=sum(!is.na(res$padj) & res$padj < 0.05),
      n_sig_lfc1=sum(sig1), n_up_lfc1=sum(sig1 & res$log2FoldChange > 0), n_down_lfc1=sum(sig1 & res$log2FoldChange < 0),
      n_sig_lfc0.58=sum(sig058), n_up_lfc0.58=sum(sig058 & res$log2FoldChange > 0), n_down_lfc0.58=sum(sig058 & res$log2FoldChange < 0),
      median_baseMean=median(res$baseMean, na.rm=TRUE), notes=""
    )
    for (cat in sort(unique(res$sqanti3_category))) {
      sub <- res[res$sqanti3_category == cat,]
      s1 <- !is.na(sub$padj) & sub$padj < 0.05 & abs(sub$log2FoldChange) >= 1 & sub$baseMean >= 10
      s058 <- !is.na(sub$padj) & sub$padj < 0.05 & abs(sub$log2FoldChange) >= 0.58 & sub$baseMean >= 10
      cat_rows[[length(cat_rows)+1]] <- data.frame(
        analysis_cohort=cohort, contrast=nm, sqanti3_category=cat, n_tested=nrow(sub),
        n_sig_lfc1=sum(s1), n_up_lfc1=sum(s1 & sub$log2FoldChange > 0), n_down_lfc1=sum(s1 & sub$log2FoldChange < 0),
        n_sig_lfc0.58=sum(s058)
      )
    }
  }
}
summary_df <- do.call(rbind, summary_rows)
cat_summary <- do.call(rbind, cat_rows)
write.table(summary_df, file.path(out, "tables", "DTE_summary_by_contrast.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(cat_summary, file.path(out, "tables", "DTE_summary_by_sqanti3_category.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

conc_rows <- list()
for (nm in names(contrast_defs)) {
  p <- primary$results[[nm]]
  s <- sensitivity$results[[nm]]
  m <- merge(p[,c("isoform_id","log2FoldChange","stat","padj","baseMean")], s[,c("isoform_id","log2FoldChange","stat","padj","baseMean")], by="isoform_id", suffixes=c("_primary","_sensitivity"))
  psig <- p$isoform_id[!is.na(p$padj) & p$padj < 0.05 & abs(p$log2FoldChange) >= 1 & p$baseMean >= 10]
  ssig <- s$isoform_id[!is.na(s$padj) & s$padj < 0.05 & abs(s$log2FoldChange) >= 1 & s$baseMean >= 10]
  overlap <- intersect(psig, ssig)
  mo <- m[m$isoform_id %in% overlap,]
  dir_conc <- ifelse(nrow(mo) > 0, mean(sign(mo$log2FoldChange_primary) == sign(mo$log2FoldChange_sensitivity), na.rm=TRUE), NA_real_)
  conc_rows[[length(conc_rows)+1]] <- data.frame(
    contrast=nm, n_primary_sig_lfc1=length(psig), n_sensitivity_sig_lfc1=length(ssig), n_overlap_sig_lfc1=length(overlap),
    direction_concordance_overlap=dir_conc,
    cor_log2FC_all_shared=cor(m$log2FoldChange_primary, m$log2FoldChange_sensitivity, use="complete.obs"),
    cor_stat_all_shared=cor(m$stat_primary, m$stat_sensitivity, use="complete.obs"),
    notes=""
  )
  write.table(m, file.path(out, "figures", paste0("DTE_log2FC_primary_vs_sensitivity_", nm, ".data.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  gp <- ggplot(m, aes(log2FoldChange_primary, log2FoldChange_sensitivity)) + geom_point(alpha=0.35, size=0.35) + theme_bw(base_size=10) + geom_hline(yintercept=0, linetype=2, color="grey50") + geom_vline(xintercept=0, linetype=2, color="grey50") + labs(title=paste0("Primary vs sensitivity log2FC: ", nm), x="primary72 log2FC", y="sensitivity75 log2FC")
  ggsave(file.path(out, "figures", paste0("DTE_log2FC_primary_vs_sensitivity_", nm, ".pdf")), gp, width=5.5, height=5.5)
}
concordance <- do.call(rbind, conc_rows)
write.table(concordance, file.path(out, "tables", "DTE_primary72_vs_sensitivity75_concordance.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

qc_all <- rbind(primary$qc, sensitivity$qc)
write.table(qc_all, file.path(out, "tables", "DTE_preQC_all_sample_metrics.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

report <- paste0(
"# Phase 02 paired DTE report\n\n",
"## Inputs\n\n",
"- Primary counts: `", counts_primary_path, "`\n",
"- Primary metadata: `", meta_primary_path, "`\n",
"- Sensitivity counts: `", counts_sens_path, "`\n",
"- Sensitivity metadata: `", meta_sens_path, "`\n",
"- Master catalog: `", catalog_path, "`\n\n",
"## Cohort definitions\n\n",
"Primary72 uses 24 complete N/A/T patients and excludes CRC_20 and D samples. Sensitivity75 uses 25 complete N/A/T patients, includes CRC_20, and excludes D samples.\n\n",
"## DTE model\n\n",
"DESeq2 isoform-level paired model: `design = ~ patient_id + tissue_group`. Contrasts: A_vs_N = A - N; T_vs_A = T - A; T_vs_N = T - N.\n\n",
"## Summary by contrast\n\n",
paste(capture.output(print(summary_df)), collapse="\n"), "\n\n",
"## Primary vs sensitivity concordance\n\n",
paste(capture.output(print(concordance)), collapse="\n"), "\n\n",
"## QC\n\n",
"Pre-DTE QC tables and plots were generated for library size, detected isoforms, size factors, VST PCA, and sample correlation. New QC flags are recorded in `tables/DTE_preQC_all_sample_metrics.tsv`; flags are not used to remove samples.\n\n",
"## Patient direction consistency\n\n",
"Patient direction consistency for main-threshold significant isoforms is reported in `tables/primary72_DTE_patient_direction_consistency.tsv` and `tables/sensitivity75_DTE_patient_direction_consistency.tsv`.\n\n",
"## SQANTI3 category distribution\n\n",
"See `tables/DTE_summary_by_sqanti3_category.tsv`.\n\n",
"## Risks and limitations\n\n",
"This is DTE, not DTU or isoform-switch analysis. Gene-level names are `NA` where unavailable from the Phase 01 catalog. CRC_20 is not reinterpreted in this phase.\n\n",
"## Next step\n\n",
"Review DTE results and QC before entering DTU/SUPPA2 or candidate prioritization phases.\n"
)
writeLines(report, file.path(out, "reports", "phase02_DTE_report.md"))

worklog <- paste0(
"# Phase 02 DTE worklog\n\n",
"## Date\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
"## Task goal\n\nDTE pre-QC and paired isoform-level DTE with DESeq2.\n\n",
"## Input files\n\nSee report and provenance JSON.\n\n",
"## Output files\n\nFull DTE tables, significant threshold tables, QC tables, plots, concordance tables, and reports under `02_DTE_paired`.\n\n",
"## Scripts and commands\n\n`Rscript code/phase02_dte_paired.R`\n\n",
"## Software environment\n\nRscript: `Rscript`; DESeq2 used.\n\n",
"## Key parameters\n\nDesign `~ patient_id + tissue_group`; thresholds padj < 0.05, |log2FC| >= 1 or 0.58, baseMean >= 10.\n\n",
"## Run status\n\nSUCCESS\n\n",
"## Errors and fixes\n\nNo sample was deleted. D and CRC_20 handling inherited from Phase 00/01.\n\n",
"## Next step\n\nStop here; do not enter DTU until this report is reviewed.\n"
)
writeLines(worklog, file.path(out, "reports", "phase02_DTE_worklog.md"))

capture.output(sessionInfo(), file=file.path(out, "provenance", "sessionInfo_phase02.txt"))
prov <- list(
  timestamp=as.character(Sys.time()),
  primary_counts=counts_primary_path,
  primary_metadata=meta_primary_path,
  sensitivity_counts=counts_sens_path,
  sensitivity_metadata=meta_sens_path,
  catalog=catalog_path,
  model="~ patient_id + tissue_group",
  contrasts=names(contrast_defs)
)
json_lines <- c(
  "{",
  paste0('  "timestamp": "', prov$timestamp, '",'),
  paste0('  "primary_counts": "', prov$primary_counts, '",'),
  paste0('  "primary_metadata": "', prov$primary_metadata, '",'),
  paste0('  "sensitivity_counts": "', prov$sensitivity_counts, '",'),
  paste0('  "sensitivity_metadata": "', prov$sensitivity_metadata, '",'),
  paste0('  "catalog": "', prov$catalog, '",'),
  paste0('  "model": "', prov$model, '",'),
  paste0('  "contrasts": ["', paste(prov$contrasts, collapse='", "'), '"]'),
  "}"
)
writeLines(json_lines, file.path(out, "provenance", "phase02_provenance.json"))
message("PHASE02_DTE_COMPLETE")
