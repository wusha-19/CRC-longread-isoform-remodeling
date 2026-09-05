#!/usr/bin/env Rscript
set.seed(20260625)

suppressPackageStartupMessages({
  library(satuRn)
  library(stageR)
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
})

root <- Sys.getenv("CRC_ROOT")
if (!nzchar(root)) stop("Set CRC_ROOT before running this script.")
phase01 <- file.path(root, "05_isoform_main_study/01_master_catalog")
phase02 <- file.path(root, "05_isoform_main_study/02_DTE_paired")
phase02b <- file.path(root, "05_isoform_main_study/02B_DTE_QC_sensitivity")
out <- file.path(root, "05_isoform_main_study/03_DTU_paired")
for (d in c("code", "config", "tables", "figures", "reports", "logs", "matrices",
            "provenance", "figures/top_DTU_usage_plots")) {
  dir.create(file.path(out, d), recursive=TRUE, showWarnings=FALSE)
}

counts_primary_path <- file.path(phase01, "matrices/counts_primary72_isoforms.tsv")
meta_primary_path <- file.path(phase01, "matrices/metadata_primary72.tsv")
counts_sens_path <- file.path(phase01, "matrices/counts_sensitivity75_isoforms.tsv")
meta_sens_path <- file.path(phase01, "matrices/metadata_sensitivity75.tsv")
catalog_path <- file.path(phase01, "tables/CRC_isoform_master_catalog.tsv")

contrast_info <- data.frame(
  contrast=c("A_vs_N", "T_vs_A", "T_vs_N"),
  group1=c("A", "T", "T"),
  group0=c("N", "A", "N"),
  stringsAsFactors=FALSE
)

read_counts <- function(path) {
  x <- read.delim(path, check.names=FALSE)
  rownames(x) <- x[[1]]
  x <- as.matrix(x[, -1, drop=FALSE])
  storage.mode(x) <- "integer"
  x
}

first_nonempty <- function(x) {
  x <- unique(as.character(x[!is.na(x) & x != ""]))
  if (length(x) == 0) return(NA_character_)
  x[1]
}

validate_inputs <- function(cohort, counts, meta) {
  if (!identical(colnames(counts), meta$sample_id)) stop(sprintf("%s sample order mismatch", cohort))
  if (any(meta$tissue_group == "D")) stop(sprintf("%s contains D samples", cohort))
  if (cohort == "primary72" && any(meta$CRC20_status == "CRC20")) stop("primary72 contains CRC_20")
  if (cohort == "sensitivity75" && !any(meta$CRC20_status == "CRC20")) stop("sensitivity75 lacks CRC_20")
  tab <- table(meta$patient_id, meta$tissue_group)
  complete <- apply(tab[, c("N","A","T"), drop=FALSE], 1, function(z) all(z == 1))
  if (!all(complete)) stop(sprintf("%s has incomplete N/A/T patient sets", cohort))
}

usage_matrix <- function(counts, map) {
  counts <- counts[map$isoform_id, , drop=FALSE]
  gene_counts <- rowsum(counts, group=map$gene_id, reorder=FALSE)
  denom <- gene_counts[map$gene_id, , drop=FALSE]
  u <- counts / pmax(denom, 1)
  u[denom == 0] <- NA_real_
  u
}

make_dtu_input <- function(cohort, counts, catalog) {
  map <- catalog[, c("isoform_id", "gene_id", "gene_name", "sqanti3_category")]
  map$gene_id <- as.character(map$gene_id)
  map <- map[map$isoform_id %in% rownames(counts), , drop=FALSE]
  n_isoforms_start <- nrow(counts)
  n_genes_start <- length(unique(map$gene_id[!is.na(map$gene_id) & map$gene_id != ""]))
  map <- map[!is.na(map$gene_id) & map$gene_id != "", , drop=FALSE]
  counts <- counts[map$isoform_id, , drop=FALSE]
  n_isoforms_with_gene_id <- nrow(map)
  n_genes_with_gene_id <- length(unique(map$gene_id))
  gene_iso_n <- table(map$gene_id)
  map <- map[map$gene_id %in% names(gene_iso_n)[gene_iso_n >= 2], , drop=FALSE]
  counts <- counts[map$isoform_id, , drop=FALSE]
  n_genes_with_2plus_isoforms <- length(unique(map$gene_id))
  n_isoforms_with_2plus_gene <- nrow(map)

  min_gene_expr <- 100
  min_feature_expr <- 10
  min_feature_prop <- 0.05
  min_gene_samples <- ceiling(ncol(counts) * 0.5)
  min_iso_samples <- ceiling(ncol(counts) / 6)

  gene_counts <- rowsum(counts, group=map$gene_id, reorder=FALSE)
  gene_keep <- rownames(gene_counts)[rowSums(gene_counts >= min_gene_expr) >= min_gene_samples]
  map <- map[map$gene_id %in% gene_keep, , drop=FALSE]
  counts <- counts[map$isoform_id, , drop=FALSE]
  u <- usage_matrix(counts, map)
  iso_keep <- rowSums(counts >= min_feature_expr) >= min_iso_samples |
    rowSums(u >= min_feature_prop, na.rm=TRUE) >= min_iso_samples
  map <- map[iso_keep, , drop=FALSE]
  counts <- counts[map$isoform_id, , drop=FALSE]
  tab <- table(map$gene_id)
  map <- map[map$gene_id %in% names(tab)[tab >= 2], , drop=FALSE]
  counts <- counts[map$isoform_id, , drop=FALSE]

  filt <- data.frame(
    analysis_cohort=cohort,
    n_isoforms_start=n_isoforms_start,
    n_genes_start=n_genes_start,
    n_isoforms_with_gene_id=n_isoforms_with_gene_id,
    n_genes_with_gene_id=n_genes_with_gene_id,
    n_genes_with_2plus_isoforms=n_genes_with_2plus_isoforms,
    n_isoforms_with_2plus_gene=n_isoforms_with_2plus_gene,
    min_gene_expr=min_gene_expr,
    min_gene_samples=min_gene_samples,
    min_feature_expr=min_feature_expr,
    min_feature_prop=min_feature_prop,
    min_iso_samples=min_iso_samples,
    n_genes_after_expression_filter=length(unique(map$gene_id)),
    n_isoforms_after_expression_filter=nrow(map)
  )
  list(counts=counts, map=map, filter=filt)
}

run_gene_de <- function(cohort, counts, meta, map) {
  gene_counts <- rowsum(counts[map$isoform_id, , drop=FALSE], group=map$gene_id, reorder=FALSE)
  meta2 <- meta
  meta2$tissue_group <- factor(meta2$tissue_group, levels=c("N","A","T"))
  meta2$patient_id <- factor(meta2$patient_id)
  rownames(meta2) <- meta2$sample_id
  dds <- DESeqDataSetFromMatrix(gene_counts, meta2, design=~ patient_id + tissue_group)
  dds <- estimateSizeFactors(dds, type="poscounts")
  dds <- DESeq(dds, quiet=FALSE)
  reslist <- list()
  for (i in seq_len(nrow(contrast_info))) {
    ci <- contrast_info[i,]
    res <- as.data.frame(results(dds, contrast=c("tissue_group", ci$group1, ci$group0), alpha=0.10))
    res$gene_id <- rownames(res)
    res$contrast <- ci$contrast
    res$analysis_cohort <- cohort
    res <- res[, c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","contrast","analysis_cohort")]
    write.table(res, file.path(out, "tables", paste0(cohort, "_gene_level_DE_", ci$contrast, ".tsv")),
                sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    reslist[[ci$contrast]] <- res
  }
  reslist
}

usage_means <- function(u, meta) {
  data.frame(
    isoform_id=rownames(u),
    mean_usage_N=rowMeans(u[, meta$sample_id[meta$tissue_group == "N"], drop=FALSE], na.rm=TRUE),
    mean_usage_A=rowMeans(u[, meta$sample_id[meta$tissue_group == "A"], drop=FALSE], na.rm=TRUE),
    mean_usage_T=rowMeans(u[, meta$sample_id[meta$tissue_group == "T"], drop=FALSE], na.rm=TRUE)
  )
}

stage_adjust_saturn <- function(pvals, map, alpha) {
  pvals[is.na(pvals)] <- 1
  split_idx <- split(seq_along(pvals), map$gene_id)
  p_gene <- vapply(split_idx, function(i) min(pvals[i], na.rm=TRUE), numeric(1))
  p_gene[!is.finite(p_gene)] <- 1
  q_screen <- p.adjust(p_gene, method="BH")
  p_conf <- matrix(pvals, ncol=1, dimnames=list(map$isoform_id, "transcript"))
  tx2gene <- data.frame(transcript=map$isoform_id, gene=map$gene_id)
  sr <- stageRTx(pScreen=q_screen, pConfirmation=p_conf, pScreenAdjusted=TRUE, tx2gene=tx2gene)
  sr <- stageWiseAdjustment(sr, method="dtu", alpha=alpha, allowNA=TRUE)
  adj <- getAdjustedPValues(sr, onlySignificantGenes=FALSE, order=FALSE)
  list(gene_p=p_gene, gene_q=q_screen, stage=adj)
}

extract_pvalue <- function(res) {
  for (nm in c("empirical_pval", "pval", "pvalue", "pValue", "p.value")) {
    if (nm %in% names(res)) return(as.numeric(res[[nm]]))
  }
  stop("No p-value column found in satuRn result")
}

sig_main_iso <- function(df) {
  !is.na(df$padj) & df$padj < 0.05 &
    !is.na(df$gene_level_padj) & df$gene_level_padj < 0.05 &
    abs(df$delta_usage) >= 0.10
}

sig_explore_iso <- function(df) {
  !is.na(df$padj_stageR_0.10) & df$padj_stageR_0.10 < 0.10 &
    !is.na(df$gene_level_padj) & df$gene_level_padj < 0.10 &
    abs(df$delta_usage) >= 0.05
}

dominant_switch_table <- function(cohort, u, meta, map) {
  rows <- list()
  for (i in seq_len(nrow(contrast_info))) {
    ci <- contrast_info[i,]
    for (gid in unique(map$gene_id)) {
      ids <- map$isoform_id[map$gene_id == gid]
      m1 <- rowMeans(u[ids, meta$sample_id[meta$tissue_group == ci$group1], drop=FALSE], na.rm=TRUE)
      m0 <- rowMeans(u[ids, meta$sample_id[meta$tissue_group == ci$group0], drop=FALSE], na.rm=TRUE)
      d1 <- names(which.max(m1)); d0 <- names(which.max(m0))
      rows[[length(rows)+1]] <- data.frame(
        analysis_cohort=cohort,
        gene_id=gid,
        gene_name=first_nonempty(map$gene_name[map$gene_id == gid]),
        contrast=ci$contrast,
        dominant_isoform_group1=d1,
        dominant_isoform_group2=d0,
        dominant_usage_group1=unname(m1[d1]),
        dominant_usage_group2=unname(m0[d0]),
        switch_status=ifelse(d1 != d0, "DOMINANT_SWITCH", "NO_DOMINANT_SWITCH"),
        max_abs_delta_usage=NA_real_
      )
    }
  }
  do.call(rbind, rows)
}

patient_consistency <- function(cohort, iso_df, u, meta) {
  sig <- iso_df[sig_main_iso(iso_df), , drop=FALSE]
  rows <- list()
  for (i in seq_len(nrow(sig))) {
    x <- sig[i,]
    ci <- contrast_info[contrast_info$contrast == x$contrast,]
    same <- 0; evaln <- 0
    for (pid in unique(meta$patient_id)) {
      s1 <- meta$sample_id[meta$patient_id == pid & meta$tissue_group == ci$group1]
      s0 <- meta$sample_id[meta$patient_id == pid & meta$tissue_group == ci$group0]
      if (length(s1) == 1 && length(s0) == 1) {
        diff <- u[x$isoform_id, s1] - u[x$isoform_id, s0]
        if (!is.na(diff) && diff != 0 && !is.na(x$delta_usage) && x$delta_usage != 0) {
          evaln <- evaln + 1
          same <- same + as.integer(sign(diff) == sign(x$delta_usage))
        }
      }
    }
    rows[[length(rows)+1]] <- data.frame(
      analysis_cohort=cohort, gene_id=x$gene_id, isoform_id=x$isoform_id,
      contrast=x$contrast, overall_delta_usage=x$delta_usage,
      n_patients_evaluable=evaln, n_patients_same_direction=same,
      patient_direction_consistency=ifelse(evaln > 0, same/evaln, NA_real_)
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

plot_usage_gene <- function(cohort, contrast, gid, u, meta, map, dom_df) {
  ids <- map$isoform_id[map$gene_id == gid]
  ids <- names(sort(rowMeans(u[ids, , drop=FALSE], na.rm=TRUE), decreasing=TRUE))[seq_len(min(6, length(ids)))]
  dat <- do.call(rbind, lapply(ids, function(id) {
    data.frame(sample_id=meta$sample_id, patient_id=meta$patient_id,
               tissue_group=meta$tissue_group, isoform_id=id, usage=as.numeric(u[id, meta$sample_id]))
  }))
  dom <- dom_df[dom_df$gene_id == gid & dom_df$contrast == contrast,]
  p <- ggplot(dat, aes(tissue_group, usage, group=patient_id)) +
    geom_line(alpha=0.25, color="grey55") +
    geom_point(aes(color=tissue_group), size=1.6, alpha=0.85) +
    facet_wrap(~ isoform_id, scales="free_y") +
    theme_bw(base_size=8) +
    labs(title=paste(cohort, contrast, gid, ifelse(nrow(dom), dom$switch_status[1], "")),
         x="", y="isoform usage")
  path <- file.path(out, "figures/top_DTU_usage_plots",
                    paste0(cohort, "_", contrast, "_", gsub("[^A-Za-z0-9_.-]", "_", gid), "_usage.pdf"))
  ggsave(path, p, width=8, height=5)
  path
}

run_cohort <- function(cohort, counts_path, meta_path, catalog) {
  message("START_COHORT ", cohort)
  counts <- read_counts(counts_path)
  meta <- read.delim(meta_path, check.names=FALSE)
  validate_inputs(cohort, counts, meta)
  meta$tissue_group <- factor(meta$tissue_group, levels=c("N","A","T"))
  meta$patient_id <- factor(meta$patient_id)
  dtu <- make_dtu_input(cohort, counts, catalog)
  write.table(dtu$filter, file.path(out, "tables", paste0(cohort, "_DTU_filter_summary.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(dtu$map, file.path(out, "matrices", paste0(cohort, "_DTU_isoform_gene_map.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  write.table(data.frame(isoform_id=rownames(dtu$counts), dtu$counts, check.names=FALSE),
              file.path(out, "matrices", paste0(cohort, "_DTU_counts_filtered.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  u <- usage_matrix(dtu$counts, dtu$map)
  write.table(data.frame(isoform_id=rownames(u), u, check.names=FALSE),
              file.path(out, "matrices", paste0(cohort, "_isoform_usage_matrix.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  mean_u <- usage_means(u, meta)
  gene_de <- run_gene_de(cohort, counts, meta, dtu$map)
  row_data <- data.frame(isoform_id=dtu$map$isoform_id, gene_id=dtu$map$gene_id)
  rownames(row_data) <- row_data$isoform_id
  col_data <- meta
  rownames(col_data) <- col_data$sample_id
  se <- SummarizedExperiment(assays=list(counts=dtu$counts[dtu$map$isoform_id, meta$sample_id, drop=FALSE]),
                             colData=col_data, rowData=row_data)
  message("SATURN_FIT ", cohort)
  se <- fitDTU(object=se, formula=~ patient_id + tissue_group, parallel=FALSE, verbose=TRUE)
  design <- model.matrix(~ patient_id + tissue_group, data=meta)
  L <- matrix(0, nrow=ncol(design), ncol=nrow(contrast_info), dimnames=list(colnames(design), contrast_info$contrast))
  L["tissue_groupA", "A_vs_N"] <- 1
  L["tissue_groupT", "T_vs_N"] <- 1
  L["tissue_groupT", "T_vs_A"] <- 1
  L["tissue_groupA", "T_vs_A"] <- -1
  message("SATURN_TEST ", cohort)
  se <- testDTU(object=se, contrasts=L, diagplot1=FALSE, diagplot2=FALSE, sort=FALSE)
  all_iso <- list(); all_gene <- list()
  dom <- dominant_switch_table(cohort, u, meta, dtu$map)
  for (i in seq_len(nrow(contrast_info))) {
    ci <- contrast_info[i,]
    res <- as.data.frame(rowData(se)[[paste0("fitDTUResult_", ci$contrast)]])
    pvals <- extract_pvalue(res)
    st05 <- stage_adjust_saturn(pvals, dtu$map, 0.05)
    st10 <- stage_adjust_saturn(pvals, dtu$map, 0.10)
    iso <- data.frame(gene_id=dtu$map$gene_id, gene_name=dtu$map$gene_name,
                      isoform_id=dtu$map$isoform_id, sqanti3_category=dtu$map$sqanti3_category,
                      contrast=ci$contrast, analysis_cohort=cohort, pvalue=pvals)
    iso <- merge(iso, mean_u, by="isoform_id", all.x=TRUE, sort=FALSE)
    iso$delta_usage <- iso[[paste0("mean_usage_", ci$group1)]] - iso[[paste0("mean_usage_", ci$group0)]]
    iso$abs_delta_usage <- abs(iso$delta_usage)
    iso$padj <- st05$stage$transcript[match(iso$isoform_id, st05$stage$txID)]
    iso$padj_stageR_0.10 <- st10$stage$transcript[match(iso$isoform_id, st10$stage$txID)]
    iso$padj[is.na(iso$padj)] <- 1
    iso$padj_stageR_0.10[is.na(iso$padj_stageR_0.10)] <- 1
    iso$gene_level_pvalue <- st05$gene_p[match(iso$gene_id, names(st05$gene_p))]
    iso$gene_level_padj <- st05$gene_q[match(iso$gene_id, names(st05$gene_q))]
    iso$dtu_status <- ifelse(sig_main_iso(iso), "DTU_FDR05_delta010",
                             ifelse(sig_explore_iso(iso), "DTU_FDR10_delta005", "not_significant"))
    max_delta <- aggregate(abs_delta_usage ~ gene_id, data=iso, max, na.rm=TRUE)
    names(max_delta)[names(max_delta) == "abs_delta_usage"] <- "max_abs_delta_usage"
    genes <- data.frame(gene_id=names(st05$gene_p), gene_level_pvalue=as.numeric(st05$gene_p),
                        gene_level_padj=as.numeric(st05$gene_q))
    genes$gene_name <- vapply(genes$gene_id, function(g) first_nonempty(dtu$map$gene_name[dtu$map$gene_id == g]), character(1))
    genes$contrast <- ci$contrast
    genes$analysis_cohort <- cohort
    genes$n_isoforms_tested <- as.integer(table(dtu$map$gene_id)[genes$gene_id])
    genes <- merge(genes, max_delta, by="gene_id", all.x=TRUE, sort=FALSE)
    de <- gene_de[[ci$contrast]]
    genes$gene_expression_padj <- de$padj[match(genes$gene_id, de$gene_id)]
    genes$gene_expression_DE_status <- ifelse(!is.na(genes$gene_expression_padj) & genes$gene_expression_padj < 0.10,
                                              "gene_DE_FDR10", "gene_expression_not_significant")
    genes$dtu_status <- ifelse(!is.na(genes$gene_level_padj) & genes$gene_level_padj < 0.05 &
                                 genes$gene_id %in% iso$gene_id[sig_main_iso(iso)],
                               "DTU_FDR05_delta010",
                               ifelse(!is.na(genes$gene_level_padj) & genes$gene_level_padj < 0.10 &
                                        genes$gene_id %in% iso$gene_id[sig_explore_iso(iso)],
                                      "DTU_FDR10_delta005", "not_significant"))
    genes$dominant_isoform_switch <- dom$switch_status[match(paste(genes$gene_id, genes$contrast), paste(dom$gene_id, dom$contrast))]
    genes <- genes[, c("gene_id","gene_name","contrast","analysis_cohort","n_isoforms_tested",
                       "gene_level_pvalue","gene_level_padj","max_abs_delta_usage",
                       "dominant_isoform_switch","gene_expression_DE_status","dtu_status","gene_expression_padj")]
    iso <- iso[, c("gene_id","gene_name","isoform_id","sqanti3_category","contrast","analysis_cohort",
                   "mean_usage_N","mean_usage_A","mean_usage_T","delta_usage","abs_delta_usage",
                   "pvalue","padj","gene_level_padj","dtu_status","padj_stageR_0.10")]
    write.table(genes, file.path(out, "tables", paste0(cohort, "_DTU_gene_", ci$contrast, ".tsv")),
                sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    write.table(iso, file.path(out, "tables", paste0(cohort, "_DTU_isoform_", ci$contrast, ".tsv")),
                sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    all_gene[[ci$contrast]] <- genes
    all_iso[[ci$contrast]] <- iso
  }
  gene_all <- do.call(rbind, all_gene)
  iso_all <- do.call(rbind, all_iso)
  dom$max_abs_delta_usage <- gene_all$max_abs_delta_usage[match(paste(dom$gene_id, dom$contrast), paste(gene_all$gene_id, gene_all$contrast))]
  write.table(dom, file.path(out, "tables", paste0(cohort, "_dominant_isoform_switch.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  silent <- gene_all[gene_all$dtu_status == "DTU_FDR05_delta010" &
                       (is.na(gene_all$gene_expression_padj) | gene_all$gene_expression_padj >= 0.10), , drop=FALSE]
  write.table(silent, file.path(out, "tables", paste0(cohort, "_gene_silent_DTU.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  cons <- patient_consistency(cohort, iso_all, u, meta)
  write.table(cons, file.path(out, "tables", paste0(cohort, "_DTU_patient_direction_consistency.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  list(cohort=cohort, filter=dtu$filter, map=dtu$map, usage=u, meta=meta,
       gene=gene_all, iso=iso_all, dominant=dom, silent=silent, consistency=cons)
}

catalog <- read.delim(catalog_path, check.names=FALSE)
primary <- run_cohort("primary72", counts_primary_path, meta_primary_path, catalog)
sensitivity <- run_cohort("sensitivity75", counts_sens_path, meta_sens_path, catalog)

filter_summary <- rbind(primary$filter, sensitivity$filter)
write.table(filter_summary, file.path(out, "tables", "DTU_filter_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
method_fallback <- data.frame(
  method_used="satuRn + stageR",
  reason="DRIMSeq + stageR was available but failed or became computationally impractical during genewise precision estimation: default DRIMSeq failed with non-finite Hessian values; add_uniform/nlminb remained stuck at genewise precision despite stricter filtering.",
  design="~ patient_id + tissue_group",
  gene_level_FDR="BH-adjusted gene screening p-values from the minimum satuRn isoform p-value per gene",
  isoform_level_FDR="stageR transcript confirmation adjusted p-values"
)
write.table(method_fallback, file.path(out, "tables", "DTU_method_fallbacks.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

summary_rows <- list(); cat_rows <- list()
for (obj in list(primary, sensitivity)) {
  for (cc in contrast_info$contrast) {
    gi <- obj$gene[obj$gene$contrast == cc,]
    ii <- obj$iso[obj$iso$contrast == cc,]
    dom <- obj$dominant[obj$dominant$contrast == cc,]
    silent <- obj$silent[obj$silent$contrast == cc,]
    main_dtu_genes <- unique(gi$gene_id[gi$dtu_status == "DTU_FDR05_delta010"])
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      analysis_cohort=obj$cohort, contrast=cc, n_genes_tested=nrow(gi), n_isoforms_tested=nrow(ii),
      n_DTU_genes_FDR05=sum(gi$dtu_status == "DTU_FDR05_delta010", na.rm=TRUE),
      n_DTU_isoforms_FDR05_delta010=sum(sig_main_iso(ii), na.rm=TRUE),
      n_DTU_genes_FDR10=sum(gi$dtu_status %in% c("DTU_FDR05_delta010","DTU_FDR10_delta005"), na.rm=TRUE),
      n_DTU_isoforms_FDR10_delta005=sum(sig_explore_iso(ii), na.rm=TRUE),
      n_gene_silent_DTU=nrow(silent),
      n_dominant_isoform_switch=sum(dom$switch_status == "DOMINANT_SWITCH" & dom$gene_id %in% main_dtu_genes, na.rm=TRUE)
    )
    for (cat in sort(unique(ii$sqanti3_category))) {
      sub <- ii[ii$sqanti3_category == cat,]
      ss <- sig_main_iso(sub)
      cat_rows[[length(cat_rows)+1]] <- data.frame(
        analysis_cohort=obj$cohort, contrast=cc, sqanti3_category=cat,
        n_isoforms_tested=nrow(sub), n_significant_isoforms=sum(ss, na.rm=TRUE),
        n_up_usage=sum(ss & sub$delta_usage > 0, na.rm=TRUE),
        n_down_usage=sum(ss & sub$delta_usage < 0, na.rm=TRUE)
      )
    }
  }
}
summary_df <- do.call(rbind, summary_rows)
cat_summary <- do.call(rbind, cat_rows)
write.table(summary_df, file.path(out, "tables", "DTU_summary_by_contrast.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(cat_summary, file.path(out, "tables", "DTU_summary_by_sqanti3_category.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

conc_rows <- list()
for (cc in contrast_info$contrast) {
  pi <- primary$iso[primary$iso$contrast == cc,]
  si <- sensitivity$iso[sensitivity$iso$contrast == cc,]
  pg <- primary$gene[primary$gene$contrast == cc,]
  sg <- sensitivity$gene[sensitivity$gene$contrast == cc,]
  m <- merge(pi[, c("gene_id","isoform_id","delta_usage","dtu_status")],
             si[, c("gene_id","isoform_id","delta_usage","dtu_status")],
             by=c("gene_id","isoform_id"), suffixes=c("_primary","_sensitivity"))
  pgenes <- unique(pg$gene_id[pg$dtu_status == "DTU_FDR05_delta010"])
  sgenes <- unique(sg$gene_id[sg$dtu_status == "DTU_FDR05_delta010"])
  overlap_genes <- intersect(pgenes, sgenes)
  sig_m <- m[m$gene_id %in% overlap_genes,]
  conc_rows[[length(conc_rows)+1]] <- data.frame(
    contrast=cc, n_primary_DTU_genes=length(pgenes), n_sensitivity_DTU_genes=length(sgenes),
    n_overlap_DTU_genes=length(overlap_genes),
    gene_overlap_fraction=ifelse(length(pgenes) > 0, length(overlap_genes)/length(pgenes), NA_real_),
    direction_concordance_isoforms=ifelse(nrow(sig_m) > 0, mean(sign(sig_m$delta_usage_primary) == sign(sig_m$delta_usage_sensitivity), na.rm=TRUE), NA_real_),
    cor_delta_usage_all_shared_isoforms=cor(m$delta_usage_primary, m$delta_usage_sensitivity, use="complete.obs"),
    n_primary_gene_silent_DTU=sum(primary$silent$contrast == cc),
    n_sensitivity_gene_silent_DTU=sum(sensitivity$silent$contrast == cc),
    notes=""
  )
  write.table(m, file.path(out, "figures", paste0("DTU_delta_usage_primary_vs_sensitivity_", cc, ".data.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  p <- ggplot(m, aes(delta_usage_primary, delta_usage_sensitivity)) +
    geom_point(alpha=0.35, size=0.35) + geom_hline(yintercept=0, linetype=2, color="grey50") +
    geom_vline(xintercept=0, linetype=2, color="grey50") + theme_bw(base_size=10) +
    labs(title=paste0("DTU delta usage primary vs sensitivity: ", cc),
         x="primary72 delta usage", y="sensitivity75 delta usage")
  ggsave(file.path(out, "figures", paste0("DTU_delta_usage_primary_vs_sensitivity_", cc, ".pdf")), p, width=5.5, height=5.5)
}
concordance <- do.call(rbind, conc_rows)
write.table(concordance, file.path(out, "tables", "DTU_primary72_vs_sensitivity75_concordance.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")

p <- ggplot(summary_df, aes(contrast, n_DTU_genes_FDR05, fill=analysis_cohort)) + geom_col(position="dodge") + theme_bw(base_size=10) + labs(title="Main-threshold DTU genes", x="", y="genes")
ggsave(file.path(out, "figures", "DTU_summary_barplot.pdf"), p, width=7, height=4.5)
p <- ggplot(summary_df, aes(contrast, n_gene_silent_DTU, fill=analysis_cohort)) + geom_col(position="dodge") + theme_bw(base_size=10) + labs(title="Gene-silent DTU", x="", y="genes")
ggsave(file.path(out, "figures", "DTU_gene_silent_summary.pdf"), p, width=7, height=4.5)
p <- ggplot(summary_df, aes(contrast, n_dominant_isoform_switch, fill=analysis_cohort)) + geom_col(position="dodge") + theme_bw(base_size=10) + labs(title="Dominant isoform switches", x="", y="genes")
ggsave(file.path(out, "figures", "DTU_dominant_switch_summary.pdf"), p, width=7, height=4.5)

iso_both <- rbind(primary$iso, sensitivity$iso)
for (cc in contrast_info$contrast) {
  dat <- iso_both[iso_both$contrast == cc,]
  dat$neglog10 <- -log10(pmax(dat$padj, .Machine$double.xmin))
  dat$sig <- sig_main_iso(dat)
  p <- ggplot(dat, aes(delta_usage, neglog10)) + geom_point(aes(color=sig), alpha=0.4, size=0.35) +
    facet_wrap(~ analysis_cohort) + scale_color_manual(values=c("FALSE"="grey70", "TRUE"="#B23A48")) +
    theme_bw(base_size=10) + labs(title=paste0("DTU delta usage volcano: ", cc), x="delta usage", y="-log10 stageR isoform FDR", color="main DTU")
  ggsave(file.path(out, "figures", paste0("DTU_delta_usage_volcano_", cc, ".pdf")), p, width=8, height=4.5)
}

top_pool <- primary$gene[primary$gene$dtu_status == "DTU_FDR05_delta010", , drop=FALSE]
if (nrow(top_pool) < 10) {
  top_pool <- rbind(top_pool, primary$gene[primary$gene$dtu_status == "DTU_FDR10_delta005", , drop=FALSE])
}
if (nrow(top_pool) == 0) {
  top_pool <- primary$gene
}
top_candidates <- top_pool[order(top_pool$gene_level_padj, -top_pool$max_abs_delta_usage), c("gene_id","contrast")]
top_candidates <- head(unique(top_candidates), 10)
for (i in seq_len(nrow(top_candidates))) {
  plot_usage_gene("primary72", top_candidates$contrast[i], top_candidates$gene_id[i],
                  primary$usage, primary$meta, primary$map, primary$dominant)
}

report <- paste0(
"# Phase 03 DTU / isoform usage switch report\n\n",
"## Inputs\n\n",
"- Primary72 counts: `", counts_primary_path, "`\n",
"- Primary72 metadata: `", meta_primary_path, "`\n",
"- Sensitivity75 counts: `", counts_sens_path, "`\n",
"- Sensitivity75 metadata: `", meta_sens_path, "`\n",
"- Master catalog: `", catalog_path, "`\n",
"- Phase02 DTE: `", phase02, "`\n",
"- Phase02B QC: `", phase02b, "`\n\n",
"## Methods\n\n",
"The final DTU analysis used `satuRn + stageR` with paired design `~ patient_id + tissue_group`. DRIMSeq + stageR was attempted first because it was available, but it failed or became computationally impractical during genewise precision estimation. The fallback reason is recorded in `tables/DTU_method_fallbacks.tsv`.\n\n",
"## Filtering Rules\n\n",
"Genes required at least two retained isoforms, gene count >=100 in at least half of samples, and isoform count >=10 or usage >=0.05 in at least one-sixth of samples.\n\n",
"```text\n", paste(capture.output(print(filter_summary)), collapse='\\n'), "\n```\n\n",
"## Contrast Directions\n\nA_vs_N = A - N; T_vs_A = T - A; T_vs_N = T - N.\n\n",
"## DTU Summary\n\n```text\n", paste(capture.output(print(summary_df)), collapse='\\n'), "\n```\n\n",
"## Primary72 vs Sensitivity75 Robustness\n\n```text\n", paste(capture.output(print(concordance)), collapse='\\n'), "\n```\n\n",
"## Gene-Silent DTU\n\nSee `tables/primary72_gene_silent_DTU.tsv` and `tables/sensitivity75_gene_silent_DTU.tsv`.\n\n",
"## Dominant Isoform Switch\n\nSee `tables/primary72_dominant_isoform_switch.tsv` and `tables/sensitivity75_dominant_isoform_switch.tsv`.\n\n",
"## Patient Direction Consistency\n\nSee `tables/primary72_DTU_patient_direction_consistency.tsv` and `tables/sensitivity75_DTU_patient_direction_consistency.tsv`.\n\n",
"## Top DTU Examples\n\nUsage plots are saved under `figures/top_DTU_usage_plots/`.\n\n",
"## Relationship With Phase02 DTE\n\nPhase03 tests within-gene usage changes and is distinct from Phase02 isoform-level DTE. Gene-level DESeq2 is auxiliary and used only for gene-silent labeling.\n\n",
"## Risks and Limitations\n\nDTU is statistical evidence of altered isoform usage, not ORF/NMD/domain or wet-lab validation. This phase did not run SUPPA2, ORF/NMD/domain annotation, or final candidate selection.\n\n",
"## Next Step\n\nReview DTU robustness and top usage plots before functional annotation or candidate prioritization.\n")
writeLines(report, file.path(out, "reports", "phase03_DTU_report.md"))

worklog <- paste0(
"# Phase 03 DTU worklog\n\n",
"## 日期\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
"## 任务目标\n\nPaired DTU / isoform usage switch analysis.\n\n",
"## 输入文件\n\nPhase01 frozen matrices and master catalog; Phase02 and Phase02B directories.\n\n",
"## 输出文件\n\nDTU gene/isoform tables, gene-silent tables, dominant switch tables, patient consistency tables, concordance, figures, usage plots, report and provenance under `03_DTU_paired`.\n\n",
"## 脚本\n\n`code/phase03_dtu_saturn.R`\n\n",
"## 命令\n\n`Rscript code/phase03_dtu_saturn.R 2>&1 | tee logs/phase03_dtu_saturn.log`\n\n",
"## 软件环境\n\nRscript: `Rscript`; packages satuRn, stageR, SummarizedExperiment, DESeq2, ggplot2.\n\n",
"## R sessionInfo\n\nSee `provenance/sessionInfo_phase03.txt`.\n\n",
"## 关键参数\n\nPaired design `~ patient_id + tissue_group`; main threshold gene FDR <0.05, isoform stageR FDR <0.05, |delta usage| >=0.10; exploratory threshold gene FDR <0.10, isoform stageR FDR <0.10, |delta usage| >=0.05.\n\n",
"## 运行状态\n\nSUCCESS\n\n",
"## 错误和修复\n\nDRIMSeq + stageR was attempted first and failed or became computationally impractical during genewise precision estimation. The final method is satuRn + stageR; no samples were removed and Phase01/02/02B outputs were not modified.\n\n",
"## 下一步\n\nStop here; do not enter SUPPA2, ORF/NMD/domain annotation, or final candidate selection until reviewed.\n")
writeLines(worklog, file.path(out, "reports", "phase03_DTU_worklog.md"))

capture.output(sessionInfo(), file=file.path(out, "provenance", "sessionInfo_phase03.txt"))
prov <- c(
  "{",
  paste0('  "timestamp": "', as.character(Sys.time()), '",'),
  '  "method": "satuRn + stageR",',
  '  "design": "~ patient_id + tissue_group",',
  paste0('  "output": "', out, '",'),
  '  "contrasts": ["A_vs_N", "T_vs_A", "T_vs_N"],',
  '  "main_threshold": "gene FDR < 0.05; isoform stageR FDR < 0.05; abs(delta usage) >= 0.10"',
  "}"
)
writeLines(prov, file.path(out, "provenance", "phase03_provenance.json"))
message("PHASE03_DTU_COMPLETE")
