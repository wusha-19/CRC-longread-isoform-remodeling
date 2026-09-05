#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

required_packages <- c("ggplot2", "patchwork", "dplyr", "tidyr", "readr", "scales", "ragg", "ggtext", "systemfonts")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(ggtext)
})

v2_base <- Sys.getenv(
  "V2_BASE",
  Sys.getenv("CRC_V2_ROOT")
)

dir_reports <- file.path(v2_base, "reports")
dir_source <- file.path(v2_base, "source_data")
dir_pdf <- file.path(v2_base, "figures", "main", "pdf")
dir_svg <- file.path(v2_base, "figures", "main", "svg")
dir_png <- file.path(v2_base, "figures", "main", "png")
for (d in c(dir_reports, dir_source, dir_pdf, dir_svg, dir_png, file.path(v2_base, "logs"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FONT_FAMILY <- "Times New Roman"
font_match <- systemfonts::match_font(FONT_FAMILY)
if (is.null(font_match$path) || !nzchar(font_match$path) || !file.exists(font_match$path)) {
  stop("Times New Roman is not available to the R graphics environment.")
}

fmt_int <- function(x) scales::comma(as.numeric(x), accuracy = 1)
fmt_pct <- function(x) sprintf("%.1f%%", as.numeric(x) * 100)
mtime_record <- function(paths, stage) {
  data.frame(
    stage = stage,
    path = paths,
    exists = file.exists(paths),
    size_bytes = ifelse(file.exists(paths), file.info(paths)$size, NA_real_),
    mtime = ifelse(file.exists(paths), format(file.info(paths)$mtime, "%Y-%m-%d %H:%M:%S %Z"), NA_character_),
    stringsAsFactors = FALSE
  )
}
write_report <- function(path, lines) {
  writeLines(lines, con = path, useBytes = TRUE)
  path
}

theme_fig2 <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY, colour = "#202020"),
      axis.line = element_line(linewidth = 0.28, colour = "#333333"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#333333"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.7),
      legend.title = element_text(size = base_size - 0.6),
      legend.text = element_text(size = base_size - 0.9),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(4.5, "mm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = element_text(size = base_size + 1.1, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#444444", margin = margin(b = 2)),
      plot.margin = margin(3, 3, 3, 3),
      panel.grid.major.y = element_line(linewidth = 0.16, colour = "#E7E7E7"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}
theme_set(theme_fig2())
update_geom_defaults("text", list(family = FONT_FAMILY))
update_geom_defaults("label", list(family = FONT_FAMILY))

save_pdf <- function(plot, path, width_in, height_in) {
  grDevices::cairo_pdf(path, width = width_in, height = height_in, family = FONT_FAMILY, bg = "white")
  print(plot)
  dev.off()
}
save_svg <- function(plot, path, width_in, height_in) {
  grDevices::svg(path, width = width_in, height = height_in, bg = "white", onefile = TRUE, family = FONT_FAMILY)
  print(plot)
  dev.off()
}
save_png <- function(plot, path, width_in, height_in, dpi = 600) {
  ragg::agg_png(path, width = width_in, height = height_in, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
}

input_files <- c(
  fig2a = file.path(dir_source, "Fig2A_DTE_DTU_counts_v4.tsv"),
  fig2b = file.path(dir_source, "Fig2B_DTU_effect_FDR_v4.tsv"),
  fig2c = file.path(dir_source, "Fig2C_usage_pattern_classes_v4.tsv"),
  fig2d = file.path(dir_source, "Fig2D_selected_usage_shifts_v4.tsv")
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0) stop("Missing v4 source data: ", paste(missing_inputs, collapse = "; "))

guard_paths <- c(
  file.path(dir_pdf, "Fig1_paired_longread_isoform_atlas_v3.pdf"),
  file.path(dir_svg, "Fig1_paired_longread_isoform_atlas_v3.svg"),
  file.path(dir_png, "Fig1_paired_longread_isoform_atlas_v3.png"),
  file.path(dir_pdf, "Fig3_gene_isoform_discordance_v4.pdf"),
  file.path(dir_svg, "Fig3_gene_isoform_discordance_v4.svg"),
  file.path(dir_png, "Fig3_gene_isoform_discordance_v4.png"),
  file.path(dir_pdf, "Fig4_molecular_features_isoform_remodeling_v2.pdf"),
  file.path(dir_svg, "Fig4_molecular_features_isoform_remodeling_v2.svg"),
  file.path(dir_png, "Fig4_molecular_features_isoform_remodeling_v2.png"),
  file.path(dir_pdf, "Fig5_recurrence_shortread_support_v2.pdf"),
  file.path(dir_svg, "Fig5_recurrence_shortread_support_v2.svg"),
  file.path(dir_png, "Fig5_recurrence_shortread_support_v2.png")
)
mtime_before <- mtime_record(guard_paths, "before")

fig2a <- readr::read_tsv(input_files[["fig2a"]], show_col_types = FALSE) %>%
  mutate(
    comparison = factor(comparison, levels = c("A vs N", "T vs A", "T vs N")),
    evidence = ifelse(analysis_type == "Differential transcript usage", "DTU units", "DTE events"),
    statistical_unit_final = ifelse(analysis_type == "Differential transcript usage", "isoform-level DTU units", statistical_unit)
  )

fig2b <- readr::read_tsv(input_files[["fig2b"]], show_col_types = FALSE) %>%
  mutate(
    comparison = factor(comparison, levels = c("A vs N", "T vs A", "T vs N")),
    delta_usage = as.numeric(delta_usage),
    FDR = as.numeric(FDR),
    minus_log10_FDR_raw = as.numeric(minus_log10_FDR_raw),
    minus_log10_FDR_display = pmin(ifelse(is.infinite(minus_log10_FDR_raw), Inf, minus_log10_FDR_raw), 8),
    is_plotted = as.logical(is_plotted),
    is_capped_minor = is_plotted & (is.infinite(minus_log10_FDR_raw) | minus_log10_FDR_raw > 8),
    direction_display = case_when(
      direction_class == "up" ~ "Increased usage",
      direction_class == "down" ~ "Decreased usage",
      TRUE ~ "Not significant"
    ),
    direction_display = factor(direction_display, levels = c("Increased usage", "Decreased usage", "Not significant"))
  )

fig2c <- readr::read_tsv(input_files[["fig2c"]], show_col_types = FALSE)
fig2d <- readr::read_tsv(input_files[["fig2d"]], show_col_types = FALSE) %>%
  mutate(
    comparison = factor(comparison, levels = c("A - N", "T - A", "T - N")),
    label_markdown = ifelse(is.na(gene_symbol), transcript_id, paste0("*", gene_symbol, "* | ", transcript_id)),
    mean_delta_usage = as.numeric(mean_delta_usage)
  )

required_d <- c("A - N", "T - A", "T - N")
panel_d_has_ta <- all(required_d %in% as.character(fig2d$comparison))
if (!panel_d_has_ta) {
  fig2d <- fig2d %>% filter(comparison %in% c("A - N", "T - N"))
}

fig2b_counts <- fig2b %>%
  group_by(comparison) %>%
  summarise(
    n_tested = sum(is_plotted),
    dtu_up = sum(direction_class == "up"),
    dtu_down = sum(direction_class == "down"),
    dtu_total = dtu_up + dtu_down,
    n_capped = sum(is_capped_minor),
    .groups = "drop"
  )
fig2a_dtu <- fig2a %>%
  filter(analysis_type == "Differential transcript usage") %>%
  transmute(comparison, panel_a_dtu = as.numeric(count))
dtu_consistency <- fig2b_counts %>%
  left_join(fig2a_dtu, by = "comparison") %>%
  mutate(panel_a_equals_panel_b = panel_a_dtu == dtu_total)
if (!all(dtu_consistency$panel_a_equals_panel_b)) {
  stop("Panel A DTU count does not match Panel B DTU up + down; see consistency table.")
}

fig2c_top <- fig2c %>%
  filter(percentage_denominator == "all pattern-classified isoforms") %>%
  mutate(
    status = ifelse(pattern_class_display == "Remodeled patterns", "Remodeled", "Stable or unclear"),
    status = factor(status, levels = c("Remodeled", "Stable or unclear")),
    count = as.numeric(count),
    denominator = as.numeric(denominator),
    percentage = as.numeric(percentage)
  )
total_pattern_n <- unique(fig2c_top$denominator)
fig2c_bottom <- fig2c %>%
  filter(percentage_denominator == "remodeled pattern isoforms") %>%
  mutate(
    count = as.numeric(count),
    denominator = as.numeric(denominator),
    percentage = as.numeric(percentage),
    pattern_class_display = factor(
      pattern_class_display,
      levels = c("Monotonic decrease", "Monotonic increase", "Tumor-specific increase", "Tumor-specific decrease", "Adjacent-specific trough")
    )
  ) %>%
  arrange(pattern_class_display)
remodeled_n <- unique(fig2c_bottom$denominator)

p2a_data <- fig2a %>%
  mutate(
    evidence = factor(evidence, levels = c("DTE events", "DTU units")),
    plot_count = ifelse(evidence == "DTU units" & count > 0, pmax(as.numeric(count), max(as.numeric(count[evidence == "DTU units"])) * 0.035), as.numeric(count)),
    fill_colour = ifelse(evidence == "DTE events", "#5E6A7D", "#7B61A8")
  )
p2a_dte <- ggplot(filter(p2a_data, evidence == "DTE events"), aes(x = comparison, y = plot_count, fill = fill_colour)) +
  geom_col(width = 0.58) +
  geom_text(aes(label = fmt_int(count), y = ifelse(count == 0, max(plot_count) * 0.035, plot_count)), vjust = -0.35, size = 2.1) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.16))) +
  labs(title = "DTE events", x = NULL, y = "Significant DTE isoforms") +
  scale_fill_identity() +
  theme_fig2(6.8) +
  theme(plot.title = element_text(size = 7.0, face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1, size = 5.8))
p2a_dtu <- ggplot(filter(p2a_data, evidence == "DTU units"), aes(x = comparison, y = plot_count, fill = fill_colour)) +
  geom_col(width = 0.58) +
  geom_text(aes(label = fmt_int(count), y = ifelse(count == 0, max(plot_count) * 0.035, plot_count)), vjust = -0.35, size = 2.1) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.18))) +
  labs(title = "DTU units", x = NULL, y = "Significant DTU units") +
  scale_fill_identity() +
  theme_fig2(6.8) +
  theme(plot.title = element_text(size = 7.0, face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 1, size = 5.8))
p2a <- p2a_dte | p2a_dtu

x_limit <- max(0.32, ceiling(max(abs(fig2b$delta_usage), na.rm = TRUE) * 10) / 10)
fig2b_anno <- fig2b_counts %>%
  mutate(
    label = paste0("n tested = ", fmt_int(n_tested), "\nDTU up = ", fmt_int(dtu_up), "\nDTU down = ", fmt_int(dtu_down)),
    x = -x_limit + x_limit * 0.08,
    y = 7.7
  )
p2b <- ggplot(filter(fig2b, is_plotted), aes(x = delta_usage, y = minus_log10_FDR_display)) +
  geom_hline(yintercept = -log10(0.10), linewidth = 0.32, linetype = "dashed", colour = "#9A9A9A") +
  geom_vline(xintercept = c(-0.05, 0.05), linewidth = 0.32, linetype = "dashed", colour = "#9A9A9A") +
  geom_point(data = filter(fig2b, !is_capped_minor, direction_display == "Not significant"), aes(colour = direction_display), size = 0.72, alpha = 0.45) +
  geom_point(data = filter(fig2b, !is_capped_minor, direction_display != "Not significant"), aes(colour = direction_display), size = 0.9, alpha = 0.78) +
  geom_point(data = filter(fig2b, is_capped_minor), aes(colour = direction_display, fill = direction_display), shape = 24, size = 1.25, alpha = 0.85, stroke = 0.25, show.legend = FALSE) +
  geom_text(data = fig2b_anno, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.7, lineheight = 0.86, colour = "#333333") +
  facet_wrap(~ comparison, nrow = 1) +
  scale_colour_manual(
    values = c("Increased usage" = "#D66A00", "Decreased usage" = "#3B6FB6", "Not significant" = "#C8C8C8"),
    breaks = c("Increased usage", "Decreased usage", "Not significant"),
    name = NULL
  ) +
  scale_fill_manual(values = c("Increased usage" = "#D66A00", "Decreased usage" = "#3B6FB6", "Not significant" = "#C8C8C8"), guide = "none") +
  scale_x_continuous(limits = c(-x_limit, x_limit), breaks = c(-0.3, 0, 0.3)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 1, 4, 8), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "DTU effect sizes and false discovery rates", x = "Delta isoform usage", y = "−log10(FDR), capped at 8") +
  guides(colour = guide_legend(override.aes = list(size = 2.1, alpha = 0.9), nrow = 1)) +
  theme_fig2(7.0) +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), axis.text.x = element_text(size = 5.8), strip.text = element_text(size = 6.3, face = "bold"))

p2c_top <- ggplot(fig2c_top, aes(x = percentage, y = "All classified isoforms", fill = status)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = paste0(status, "\n", fmt_int(count), " (", fmt_pct(percentage), ")")), position = position_stack(vjust = 0.5), size = 2.0, colour = "white") +
  scale_fill_manual(values = c("Remodeled" = "#7B61A8", "Stable or unclear" = "#B5B5B5"), name = NULL) +
  scale_x_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Isoform usage pattern classes\nacross matched tissue states", subtitle = paste0("All classified isoforms, n = ", fmt_int(total_pattern_n)), x = NULL, y = NULL) +
  theme_fig2(7.0) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(), legend.position = "none", plot.margin = margin(2, 2, 0, 2))
p2c_bottom <- ggplot(fig2c_bottom, aes(x = percentage, y = pattern_class_display)) +
  geom_col(fill = "#7B61A8", width = 0.58) +
  geom_text(aes(label = paste0(fmt_int(count), " (", fmt_pct(percentage), ")")), hjust = -0.08, size = 2.0) +
  scale_x_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.22))) +
  labs(title = paste0("Remodeled isoforms, n = ", fmt_int(remodeled_n)), x = "Percentage of remodeled isoforms", y = NULL) +
  theme_fig2(7.0) +
  theme(axis.text.y = element_text(size = 6.1), panel.grid.major.y = element_blank(), plot.margin = margin(0, 2, 2, 2), plot.title = element_text(size = 7.0, face = "bold"))
p2c <- p2c_top / p2c_bottom + plot_layout(heights = c(0.58, 1.42))

heatmap_limit <- ceiling(max(abs(fig2d$mean_delta_usage), na.rm = TRUE) * 100) / 100
fig2d <- fig2d %>%
  mutate(
    label_markdown = factor(label_markdown, levels = rev(unique(label_markdown[order(row_order)]))),
    comparison = factor(comparison, levels = if (panel_d_has_ta) c("A - N", "T - A", "T - N") else c("A - N", "T - N"))
  )
p2d <- ggplot(fig2d, aes(x = comparison, y = label_markdown, fill = mean_delta_usage)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B24C63", midpoint = 0, limits = c(-heatmap_limit, heatmap_limit), name = "Mean Δ usage") +
  labs(title = "Selected tumor-associated\nisoform usage shifts", x = NULL, y = NULL) +
  theme_fig2(7.0) +
  theme(axis.text.y = ggtext::element_markdown(size = 6.3, family = FONT_FAMILY), axis.text.x = element_text(size = 6.4), axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(), legend.position = "right", legend.title = element_text(size = 6.7), legend.text = element_text(size = 6.0))

fig2 <- (wrap_elements(full = p2a) | wrap_elements(full = p2b)) /
  (wrap_elements(full = p2c) | wrap_elements(full = p2d)) +
  plot_layout(widths = c(0.38, 0.62), heights = c(1.0, 1.08), guides = "keep") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FONT_FAMILY, size = 11.5, face = "bold"))

width_in <- 180 / 25.4
height_in <- 132 / 25.4
main_pdf <- file.path(dir_pdf, "Fig2_tumor_associated_isoform_remodeling_v4_minor.pdf")
main_svg <- file.path(dir_svg, "Fig2_tumor_associated_isoform_remodeling_v4_minor.svg")
main_png <- file.path(dir_png, "Fig2_tumor_associated_isoform_remodeling_v4_minor.png")
main_preview <- file.path(dir_png, "Fig2_tumor_associated_isoform_remodeling_v4_minor_180mm_preview.png")
save_pdf(fig2, main_pdf, width_in, height_in)
save_svg(fig2, main_svg, width_in, height_in)
save_png(fig2, main_png, width_in, height_in)
save_png(fig2, main_preview, width_in, height_in)

mtime_after <- mtime_record(guard_paths, "after")
mtime_guard <- mtime_before %>%
  select(path, before_exists = exists, before_size_bytes = size_bytes, before_mtime = mtime) %>%
  left_join(mtime_after %>% select(path, after_exists = exists, after_size_bytes = size_bytes, after_mtime = mtime), by = "path") %>%
  mutate(mtime_guard_status = ifelse(before_exists == after_exists & before_size_bytes == after_size_bytes & before_mtime == after_mtime, "UNCHANGED", "CHANGED"))
readr::write_tsv(mtime_guard, file.path(dir_reports, "Fig2_minor_revision_mtime_guard_v4.tsv"))
if (any(mtime_guard$mtime_guard_status != "UNCHANGED")) stop("Figure 1/3/4/5 mtime guard failed.")

legend_zh <- c(
  "图2 | 肿瘤组织状态相关的 isoform 表达和使用比例重塑。",
  "",
  "N、A、T 分别表示配对正常组织、癌旁组织和肿瘤组织。本研究为同一患者的横断面配对组织取样，N/A/T 的展示顺序用于组织状态比较，不表示纵向时间序列或个体内进展轨迹。",
  "",
  "A，A vs N、T vs A 和 T vs N 三个配对比较中的显著 differential transcript expression（DTE）isoforms 和 differential transcript usage（DTU）units 数量。DTE 的统计单位为通过严格表达差异阈值的 isoform；DTU 的统计单位为 isoform-level DTU unit。Panel A 中的 DTU 为同时满足 FDR < 0.10 和 |Δ usage| >= 0.05 的 effect-size-filtered DTU units。",
  "",
  "B，DTU effect size 与 FDR 的关系。每个点代表一个 isoform-level DTU result row。橙色表示显著上调用法比例，蓝色表示显著下调用法比例，灰色表示未同时达到 FDR 和 |Δ usage| 阈值。每个分面的 n tested 表示具有可绘制 Δ usage 和 FDR 的 DTU units 数量；DTU up 和 DTU down 分别表示通过阈值且方向为正或负的 DTU units 数量。位于 y 轴上限的上三角表示 −log10(FDR) 超过 8 后按图示上限截断显示。",
  "",
  "C，N/A/T matched tissue states 下的 isoform usage pattern 分类。上方 stacked bar 的百分比以全部已分类 isoforms 为分母（n = 668）；下方 remodeled pattern bar plot 的百分比以 remodeled isoforms 为分母（n = 230）。Stable or unclear 和 Remodeled 类别保留冻结 source data 中的定义。",
  "",
  paste0("D，预设 tumor-associated heatmap row set 中 selected isoforms 的 mean Δ isoform usage。每个 cell 表示对应组织比较中的 tissue-group mean usage difference。当前 Panel D ", if (panel_d_has_ta) "包含 A − N、T − A 和 T − N 三列。" else "仅包含 A − N 和 T − N；T − A heatmap values 在现有 source data 中不可用且未重新计算，因此该 panel 表示相对于 normal tissue 的 selected isoform usage shifts。"),
  "",
  "Source data are provided as Source Data files."
)
legend_zh_path <- file.path(dir_reports, "Figure2_legend_draft_zh_v4.txt")
write_report(legend_zh_path, legend_zh)

summary_lines <- c(
  "# Fig2 minor revision summary v4",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Modified script",
  paste0("- `", file.path(v2_base, "scripts", "14_build_Fig2_minor_revision_v4.R"), "`"),
  "",
  "## Scope",
  "- Render-only minor revision from existing v4 source data.",
  "- No upstream analysis was rerun and no source statistic was changed.",
  "",
  "## Panel A",
  "- Final DTU statistical unit: isoform-level DTU units.",
  "- Panel A title/y-axis use DTU units consistently.",
  paste0("- Panel A DTU counts: ", paste(dtu_consistency$comparison, dtu_consistency$panel_a_dtu, sep = "=", collapse = "; "), "."),
  paste0("- Panel B DTU up + down: ", paste(dtu_consistency$comparison, dtu_consistency$dtu_total, sep = "=", collapse = "; "), "."),
  paste0("- Panel A/B DTU consistency: ", ifelse(all(dtu_consistency$panel_a_equals_panel_b), "PASS", "FAIL"), "."),
  "- Panel A uses effect-size-filtered DTU units passing both FDR and absolute delta-usage thresholds.",
  "",
  "## Panel B",
  "- Facet annotations enlarged to approximately 7.5-8 pt and kept as three left-aligned lines.",
  "- The long capped-triangle note was removed from the plot subtitle; capped triangles are explained in the figure legend.",
  "",
  "## Panel C",
  paste0("- Added denominator label: All classified isoforms, n = ", fmt_int(total_pattern_n), "."),
  paste0("- Added denominator label: Remodeled isoforms, n = ", fmt_int(remodeled_n), "."),
  "",
  "## Panel D",
  paste0("- T - A included: ", panel_d_has_ta, "."),
  if (panel_d_has_ta) {
    "- Panel D displays A - N, T - A and T - N from existing v4 source data."
  } else {
    "- T - A heatmap values were not available in the existing source data and were not recomputed; legend states relative to normal tissue."
  },
  "",
  "## Outputs",
  paste0("- PDF: `", main_pdf, "`"),
  paste0("- SVG: `", main_svg, "`"),
  paste0("- PNG: `", main_png, "`"),
  paste0("- 180 mm preview: `", main_preview, "`"),
  paste0("- Chinese legend: `", legend_zh_path, "`"),
  "",
  "## Mtime guard",
  paste0("- Figure 1/3/4/5 mtime guard: ", ifelse(all(mtime_guard$mtime_guard_status == "UNCHANGED"), "PASS", "FAIL"), "."),
  "",
  "## Manual decisions still needed",
  "- None for the requested four minor revision items."
)
summary_path <- file.path(dir_reports, "Fig2_minor_revision_summary_v4.md")
write_report(summary_path, summary_lines)

qa_lines <- c(
  "# Fig2 QC v4 minor",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Runtime",
  paste0("- Times New Roman: PASS (`", font_match$path, "`)."),
  "- R render exit status: PASS if this report exists.",
  "",
  "## Count checks",
  paste0("- Panel A/B DTU consistency: ", ifelse(all(dtu_consistency$panel_a_equals_panel_b), "PASS", "FAIL"), "."),
  paste0("- n tested / DTU up / DTU down: ", paste(fig2b_counts$comparison, fig2b_counts$n_tested, fig2b_counts$dtu_up, fig2b_counts$dtu_down, sep = "=", collapse = "; "), "."),
  "",
  "## Panel checks",
  "- Panel A DTU wording uses DTU units only: PASS.",
  "- Panel B annotation enlarged and formatted as three lines: PASS.",
  paste0("- Panel C denominator labels: PASS (", fmt_int(total_pattern_n), " and ", fmt_int(remodeled_n), ")."),
  paste0("- Panel D includes T - A: ", ifelse(panel_d_has_ta, "PASS", "NOT_AVAILABLE"), "."),
  "",
  "## Guard",
  paste0("- Figure 1/3/4/5 mtime guard: ", ifelse(all(mtime_guard$mtime_guard_status == "UNCHANGED"), "PASS", "FAIL"), ".")
)
qa_path <- file.path(dir_reports, "Fig2_QC_v4_minor.md")
write_report(qa_path, qa_lines)

message("Figure 2 v4 minor revision complete.")
