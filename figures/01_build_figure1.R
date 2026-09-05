#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "scales", "patchwork",
  "cowplot", "circlize", "ggplotify", "svglite", "ragg", "systemfonts", "Cairo"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
})

BASE <- Sys.getenv(
  "V2_BASE",
  Sys.getenv("CRC_V2_ROOT")
)
MAIN_STUDY <- dirname(BASE)
FONT_FAMILY <- "Times New Roman"

DIR_REPORTS <- file.path(BASE, "reports")
DIR_LOGS <- file.path(BASE, "logs")
DIR_SOURCE <- file.path(BASE, "source_data")
DIR_PDF <- file.path(BASE, "figures", "main", "pdf")
DIR_SVG <- file.path(BASE, "figures", "main", "svg")
DIR_PNG <- file.path(BASE, "figures", "main", "png")
for (d in c(DIR_REPORTS, DIR_LOGS, DIR_SOURCE, DIR_PDF, DIR_SVG, DIR_PNG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

INPUT <- list(
  clinical = file.path(BASE, "clinical", "clinical_metadata_deidentified_v2.tsv"),
  age_bmi = file.path(BASE, "clinical", "Fig1G_age_BMI_from_source_excel_deidentified.tsv"),
  existing_patient_profile = file.path(BASE, "tables", "Fig1G_clinical_cohort_profile_source_data.tsv"),
  sample_metadata_75 = file.path(MAIN_STUDY, "01_master_catalog", "matrices", "metadata_sensitivity75.tsv"),
  sample_metrics_75 = file.path(MAIN_STUDY, "02_DTE_paired", "tables", "sensitivity75_DTE_preQC_sample_metrics.tsv"),
  atlas = file.path(BASE, "figures", "source_data", "Fig1C_v2_isoform_catalog_scale_source_data.tsv"),
  sqanti3 = file.path(BASE, "tables", "Fig1D_SQANTI3_structural_composition_source_data.tsv"),
  multiplicity = file.path(BASE, "tables", "Fig1E_isoforms_per_gene_distribution_source_data.tsv"),
  phase01_report = file.path(MAIN_STUDY, "01_master_catalog", "reports", "phase01_master_catalog_report.md"),
  mtime_baseline = file.path(DIR_LOGS, "Fig1_v3_pre_render_mtime_20260714.tsv")
)

input_status <- tibble(
  input = names(INPUT),
  path = unname(unlist(INPUT)),
  exists = file.exists(unname(unlist(INPUT))),
  size_bytes = ifelse(
    file.exists(unname(unlist(INPUT))),
    file.info(unname(unlist(INPUT)))$size,
    NA_real_
  )
)
if (any(!input_status$exists) || any(is.na(input_status$size_bytes) | input_status$size_bytes <= 0)) {
  readr::write_tsv(input_status, file.path(DIR_REPORTS, "Fig1_v3_input_status_failed.tsv"))
  stop("Figure 1 v3 required inputs are missing or empty.")
}

font_match <- systemfonts::match_fonts(FONT_FAMILY)
if (nrow(font_match) == 0 || is.na(font_match$path[1]) || !nzchar(font_match$path[1])) {
  stop("Times New Roman is not available in the rendering environment.")
}
Cairo::CairoFonts(
  regular = "Times New Roman:style=Regular",
  bold = "Times New Roman:style=Bold",
  italic = "Times New Roman:style=Italic",
  bolditalic = "Times New Roman:style=Bold Italic"
)

is_known_tnr_metric_warning <- function(w) {
  grepl(
    "font family 'Times New Roman' not found in PostScript font database",
    conditionMessage(w),
    fixed = TRUE
  )
}

globalCallingHandlers(
  warning = function(w) {
    if (is_known_tnr_metric_warning(w)) invokeRestart("muffleWarning")
  }
)

muffle_known_tnr_metric_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (is_known_tnr_metric_warning(w)) invokeRestart("muffleWarning")
    }
  )
}

read_tsv_char <- function(path) {
  readr::read_tsv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | toupper(x) %in% c("NA", "N/A", "NULL", ".")] <- NA_character_
  x
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
format_n <- function(x) scales::comma(as.numeric(x), accuracy = 1)
format_pct1 <- function(x) scales::percent(as.numeric(x), accuracy = 0.1)
pt_size <- function(x) x / ggplot2::.pt
normalize_patient_id <- function(x) {
  x <- gsub("_", "-", clean_text(x))
  sub("^CRC-0+([0-9]+)$", "CRC-\\1", x)
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- suppressWarnings(system2("sha256sum", path, stdout = TRUE, stderr = TRUE))
  if (length(out) == 0 || !grepl("^[0-9a-fA-F]{64}", out[1])) return(NA_character_)
  sub("^([0-9a-fA-F]{64}).*$", "\\1", out[1])
}

stop_equal <- function(actual, expected, label, tolerance = 1e-8) {
  if (length(actual) != length(expected) ||
      !isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = tolerance))) {
    stop(
      label, " differs from the frozen result. Actual: ", paste(actual, collapse = ", "),
      "; expected: ", paste(expected, collapse = ", ")
    )
  }
}

write_tsv_v3 <- function(x, filename) {
  path <- file.path(DIR_SOURCE, filename)
  readr::write_tsv(x, path, na = "NA")
  path
}

write_report_tsv <- function(x, filename) {
  path <- file.path(DIR_REPORTS, filename)
  readr::write_tsv(x, path, na = "NA")
  path
}

TISSUE_COLOURS <- c(
  Normal = "#3F77B5",
  Adjacent = "#45A99A",
  Tumor = "#D95F02"
)
AGE_COLOUR <- "#7561A8"
BMI_COLOUR <- "#C46D5E"
MISSING_COLOUR <- "#D2D2D2"
CLINICAL_COLOURS <- c(
  "Sex: Female" = "#C9859E",
  "Sex: Male" = "#6F82A7",
  "Morphology: Ulcerative" = "#9B7E61",
  "Morphology: Protruding" = "#D0AA70",
  "Morphology: Other/unspecified" = "#BFC1C4",
  "Site: Colon" = "#7A6F9B",
  "Site: Rectum" = "#5E8F78",
  "Site: Unspecified" = "#BFC1C4"
)
SQANTI_COLOURS <- c(
  FSM = "#526D82",
  ISM = "#8B8B8B",
  NIC = "#5B8E8A",
  NNC = "#B07A4A",
  Other = "#D8D8D8"
)
MULTIPLICITY_COLOUR <- "#657F96"

theme_fig1 <- function(base_size = 8.0) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(linewidth = 0.32, colour = "#2B2B2B"),
      axis.ticks = element_line(linewidth = 0.32, colour = "#2B2B2B"),
      axis.title = element_text(size = 8.4, colour = "#222222"),
      axis.text = element_text(size = 7.8, colour = "#333333"),
      legend.title = element_text(size = 7.4, face = "bold"),
      legend.text = element_text(size = 7.0),
      panel.grid = element_blank(),
      plot.margin = margin(2, 3, 2, 2)
    )
}
theme_set(theme_fig1())
update_geom_defaults("text", list(family = FONT_FAMILY))
update_geom_defaults("label", list(family = FONT_FAMILY))

# -----------------------------------------------------------------------------
# Input freeze and fail-fast audit
# -----------------------------------------------------------------------------

clinical_raw <- read_tsv_char(INPUT$clinical)
age_bmi_raw <- read_tsv_char(INPUT$age_bmi)
existing_profile <- read_tsv_char(INPUT$existing_patient_profile)
sample_meta_raw <- read_tsv_char(INPUT$sample_metadata_75)
sample_metrics_raw <- read_tsv_char(INPUT$sample_metrics_75)
atlas_raw <- read_tsv_char(INPUT$atlas)
sqanti_raw <- read_tsv_char(INPUT$sqanti3)
multiplicity_raw <- read_tsv_char(INPUT$multiplicity)

clinical_required <- c("patient_id", "sex", "site_group", "tumor_form", "manual_review_flags")
age_bmi_required <- c(
  "patient_id", "patient_number", "age_years", "age_status",
  "BMI_source_cached", "BMI_recomputed", "BMI_final", "BMI_status", "audit_flags"
)
sample_meta_required <- c("sample_id", "patient_id", "tissue_group")
sample_metric_required <- c("sample_id", "patient_id", "tissue_group", "detected_isoforms", "library_size")
if (length(setdiff(clinical_required, names(clinical_raw))) > 0 ||
    length(setdiff(age_bmi_required, names(age_bmi_raw))) > 0 ||
    length(setdiff(sample_meta_required, names(sample_meta_raw))) > 0 ||
    length(setdiff(sample_metric_required, names(sample_metrics_raw))) > 0) {
  stop("One or more Figure 1 v3 inputs lack required columns.")
}

patient_profile <- clinical_raw %>%
  transmute(
    patient_id = clean_text(patient_id),
    sex_raw = clean_text(sex),
    site_raw = clean_text(site_group),
    morphology_raw = clean_text(tumor_form),
    manual_review_flags = clean_text(manual_review_flags)
  ) %>%
  left_join(
    age_bmi_raw %>%
      transmute(
        patient_id = clean_text(patient_id),
        patient_number = as.integer(as_num(patient_number)),
        age = as_num(age_years),
        age_status_source = clean_text(age_status),
        BMI_source_cached = as_num(BMI_source_cached),
        BMI_recomputed = as_num(BMI_recomputed),
        BMI = as_num(BMI_final),
        BMI_status_source = clean_text(BMI_status),
        age_BMI_audit_flags = clean_text(audit_flags)
      ),
    by = "patient_id"
  ) %>%
  mutate(
    age = ifelse(age_status_source == "available" & is.finite(age) & between(age, 18, 100), age, NA_real_),
    BMI = ifelse(BMI_status_source == "available_recomputed" & is.finite(BMI), BMI, NA_real_),
    sex = case_when(
      sex_raw == "Female" ~ "Female",
      sex_raw == "Male" ~ "Male",
      TRUE ~ "Unknown/Missing"
    ),
    tumor_site = case_when(
      site_raw == "Colon" ~ "Colon",
      site_raw == "Rectum" ~ "Rectum",
      site_raw == "Colorectum_unspecified" ~ "Unspecified",
      TRUE ~ "Unknown/Missing"
    ),
    tumor_morphology = case_when(
      grepl("溃疡", morphology_raw) ~ "Ulcerative",
      grepl("隆起", morphology_raw) ~ "Protruding",
      TRUE ~ "Other/unspecified"
    )
  ) %>%
  arrange(patient_number) %>%
  mutate(
    display_patient_id = sprintf("P%02d", patient_number),
    patient_order = row_number(),
    patient_order_rule = "Natural CRC patient-number order; no clinical-variable sorting"
  )

expected_ids <- paste0("CRC-", seq_len(25))
if (nrow(patient_profile) != 25 || n_distinct(patient_profile$patient_id) != 25 ||
    !identical(patient_profile$patient_id, expected_ids)) {
  stop("Patient profile must contain CRC-1 through CRC-25 exactly once in natural order.")
}
if (sum(is.finite(patient_profile$age)) != 25 || sum(is.finite(patient_profile$BMI)) != 24) {
  stop("Age/BMI availability differs from the frozen 25/24 audit.")
}
if (!is.na(patient_profile$BMI[patient_profile$patient_id == "CRC-8"])) {
  stop("CRC-8 BMI must remain missing and must not be imputed.")
}
if (!isTRUE(all.equal(
  patient_profile$BMI[patient_profile$patient_id == "CRC-7"],
  26.638918,
  tolerance = 1e-6
))) {
  stop("CRC-7 audited BMI differs from the frozen value.")
}

sample_meta <- sample_meta_raw %>%
  transmute(
    sample_id = clean_text(sample_id),
    patient_id = normalize_patient_id(patient_id),
    tissue_group = clean_text(tissue_group)
  )
sample_metrics <- sample_metrics_raw %>%
  transmute(
    sample_id = clean_text(sample_id),
    metrics_patient_id = normalize_patient_id(patient_id),
    metrics_tissue_group = clean_text(tissue_group),
    detected_isoforms = as_num(detected_isoforms),
    library_size = as_num(library_size)
  )

sample_detection <- sample_meta %>%
  inner_join(sample_metrics, by = "sample_id") %>%
  mutate(
    patient_match = patient_id == metrics_patient_id,
    tissue_match = tissue_group == metrics_tissue_group,
    tissue_label = recode(tissue_group, N = "Normal", A = "Adjacent", T = "Tumor")
  ) %>%
  select(sample_id, patient_id, tissue_group, tissue_label, detected_isoforms, library_size,
         patient_match, tissue_match)

if (nrow(sample_detection) != 75 || n_distinct(sample_detection$sample_id) != 75 ||
    any(!sample_detection$patient_match) || any(!sample_detection$tissue_match) ||
    any(!sample_detection$tissue_group %in% c("N", "A", "T"))) {
  stop("Official sensitivity75 metadata and pre-QC sample metrics do not form a valid 75-sample join.")
}

pair_audit <- sample_detection %>%
  count(patient_id, tissue_group, name = "n_samples") %>%
  complete(patient_id = expected_ids, tissue_group = c("N", "A", "T"), fill = list(n_samples = 0))
if (nrow(pair_audit) != 75 || any(pair_audit$n_samples != 1)) {
  stop("Every patient must have exactly one Normal, Adjacent and Tumor sample.")
}

sample_wide <- sample_detection %>%
  select(patient_id, tissue_group, sample_id) %>%
  pivot_wider(names_from = tissue_group, values_from = sample_id, names_prefix = "sample_")

patient_profile <- patient_profile %>%
  left_join(sample_wide, by = "patient_id") %>%
  rename(
    normal_sample_id = sample_N,
    adjacent_sample_id = sample_A,
    tumor_sample_id = sample_T
  ) %>%
  mutate(
    complete_NAT_pair = !is.na(normal_sample_id) & !is.na(adjacent_sample_id) & !is.na(tumor_sample_id),
    missing_fields = case_when(
      is.na(BMI) ~ "BMI",
      TRUE ~ "none"
    )
  )
if (!all(patient_profile$complete_NAT_pair)) stop("One or more patients lack a complete N/A/T pair.")

continuous_summary <- bind_rows(
  tibble(variable = "Age", unit = "years", values = list(patient_profile$age)),
  tibble(variable = "BMI", unit = "kg/m²", values = list(patient_profile$BMI))
) %>%
  rowwise() %>%
  mutate(
    n = sum(is.finite(values)),
    minimum = min(values, na.rm = TRUE),
    q1 = as.numeric(quantile(values, 0.25, na.rm = TRUE, names = FALSE)),
    median = median(values, na.rm = TRUE),
    q3 = as.numeric(quantile(values, 0.75, na.rm = TRUE, names = FALSE)),
    maximum = max(values, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(-values)

categorical_summary <- bind_rows(
  patient_profile %>% count(variable = "Sex", category = sex, name = "n"),
  patient_profile %>% count(variable = "Tumor morphology", category = tumor_morphology, name = "n"),
  patient_profile %>% count(variable = "Tumor site", category = tumor_site, name = "n")
) %>%
  group_by(variable) %>%
  mutate(percentage = n / sum(n)) %>%
  ungroup()

expected_categorical <- tribble(
  ~variable, ~category, ~expected_n,
  "Sex", "Female", 12,
  "Sex", "Male", 13,
  "Tumor morphology", "Ulcerative", 16,
  "Tumor morphology", "Protruding", 8,
  "Tumor morphology", "Other/unspecified", 1,
  "Tumor site", "Colon", 12,
  "Tumor site", "Rectum", 12,
  "Tumor site", "Unspecified", 1
)
categorical_check <- expected_categorical %>%
  left_join(categorical_summary, by = c("variable", "category"))
if (any(is.na(categorical_check$n)) || any(categorical_check$n != categorical_check$expected_n)) {
  stop("Clinical categorical distributions differ from the frozen audit.")
}
closure_check <- categorical_summary %>%
  group_by(variable) %>%
  summarise(n = sum(n), percentage = sum(percentage), .groups = "drop")
if (any(closure_check$n != 25) || any(abs(closure_check$percentage - 1) > 1e-10)) {
  stop("One or more clinical categorical summaries do not close to 25 patients and 100%.")
}

atlas <- atlas_raw %>%
  transmute(
    metric = clean_text(metric),
    value = as_num(value),
    source_metric = clean_text(source_metric),
    note = clean_text(note)
  ) %>%
  mutate(
    display_order = match(metric, c("Master catalog isoforms", "Quantified isoforms", "Genes represented")),
    display_label = recode(
      metric,
      `Master catalog isoforms` = "Catalog isoforms",
      `Quantified isoforms` = "Quantified isoforms",
      `Genes represented` = "Represented genes"
    ),
    unit = ifelse(metric == "Genes represented", "genes", "isoforms")
  ) %>%
  arrange(display_order)
stop_equal(atlas$value, c(66805, 63623, 16530), "Atlas overview values")

sqanti <- sqanti_raw %>%
  transmute(
    structural_class = clean_text(structural_class),
    n_isoforms = as_num(n_isoforms),
    percentage = as_num(percent)
  ) %>%
  mutate(
    display_order = match(structural_class, c("FSM", "ISM", "NIC", "NNC", "Other"))
  ) %>%
  arrange(display_order)
if (any(is.na(sqanti$display_order))) stop("Unexpected SQANTI3 category.")
stop_equal(sqanti$n_isoforms, c(31143, 26460, 5998, 2289, 915), "SQANTI3 counts")
if (sum(sqanti$n_isoforms) != 66805 || abs(sum(sqanti$percentage) - 1) > 1e-10) {
  stop("SQANTI3 counts or percentages do not close to the 66,805-isoform catalog.")
}
sqanti <- sqanti %>%
  mutate(
    xmin = lag(cumsum(percentage), default = 0),
    xmax = cumsum(percentage),
    display_label = paste0(structural_class, "  ", format_n(n_isoforms), " (", format_pct1(percentage), ")")
  )

multiplicity <- multiplicity_raw %>%
  transmute(
    source_category = clean_text(diversity_bin),
    n_genes = as_num(n_genes),
    percentage = as_num(percent_genes),
    denominator = as_num(total_genes_counted)
  ) %>%
  mutate(
    display_order = match(source_category, c("1 isoform", "2-3 isoforms", "4-5 isoforms", ">=6 isoforms")),
    display_category = recode(
      source_category,
      `1 isoform` = "1 isoform",
      `2-3 isoforms` = "2–3 isoforms",
      `4-5 isoforms` = "4–5 isoforms",
      `>=6 isoforms` = "≥6 isoforms"
    )
  ) %>%
  arrange(display_order)
if (any(is.na(multiplicity$display_order))) stop("Unexpected isoform-multiplicity category.")
stop_equal(multiplicity$n_genes, c(4952, 4932, 2841, 3805), "Isoform multiplicity counts")
if (sum(multiplicity$n_genes) != 16530 || abs(sum(multiplicity$percentage) - 1) > 1e-10 ||
    any(multiplicity$denominator != 16530)) {
  stop("Isoform multiplicity does not close to 16,530 represented genes and 100%.")
}

# Verify the rebuilt clinical profile against the previously frozen display table.
existing_profile_check <- existing_profile %>%
  transmute(
    patient_id = clean_text(patient_id),
    old_age = as_num(age_display),
    old_BMI = as_num(BMI_display),
    old_sex = clean_text(sex),
    old_morphology = clean_text(tumor_morphology),
    old_site = clean_text(site)
  ) %>%
  right_join(
    patient_profile %>%
      select(patient_id, age, BMI, sex, tumor_morphology, tumor_site),
    by = "patient_id"
  )
if (nrow(existing_profile_check) != 25 ||
    !isTRUE(all.equal(existing_profile_check$old_age, existing_profile_check$age, tolerance = 1e-8)) ||
    !isTRUE(all.equal(existing_profile_check$old_BMI, existing_profile_check$BMI, tolerance = 1e-8)) ||
    any(existing_profile_check$old_sex != existing_profile_check$sex) ||
    any(existing_profile_check$old_morphology != existing_profile_check$tumor_morphology) ||
    any(existing_profile_check$old_site != existing_profile_check$tumor_site)) {
  stop("Rebuilt patient profile differs from the frozen Figure 1G source table.")
}

# -----------------------------------------------------------------------------
# Source-data and reproducibility mapping
# -----------------------------------------------------------------------------

patient_check <- patient_profile %>%
  transmute(
    patient_id,
    display_patient_id,
    patient_order,
    age,
    BMI,
    sex,
    morphology_or_differentiation = tumor_morphology,
    tumor_site,
    normal_sample_id,
    adjacent_sample_id,
    tumor_sample_id,
    complete_NAT_pair,
    missing_fields
  )

track_definition <- tribble(
  ~track_order, ~track_position, ~track_name, ~variable, ~source_file, ~source_column, ~data_type, ~encoding, ~scale_or_categories, ~color, ~missing_definition, ~legend_text, ~patient_order_rule,
  1, "outer", "Patient ID", "patient_id", INPUT$age_bmi, "patient_id; patient_number", "identifier", "Tangential P01-P25 display label linked to the frozen CRC patient ID", "P01 to P25", "#2B2B2B", "Not applicable", "P01-P25 are compact de-identified display labels; the source-data table preserves the corresponding CRC patient IDs.", "Natural CRC patient-number order; no clinical-variable sorting",
  2, "outer-middle", "Age", "age", INPUT$age_bmi, "age_years", "continuous", "Radial lollipop; line and filled point", "30-90 years; actual values", AGE_COLOUR, "No missing age values", "Purple radial lollipop indicates age in years.", "Natural CRC patient-number order; no clinical-variable sorting",
  3, "middle", "BMI", "BMI", INPUT$age_bmi, "BMI_final", "continuous", "Radial lollipop; line and filled point; open grey point for missing", "15-30 kg/m²; actual values", BMI_COLOUR, "CRC-8 is NA; open grey point; never zero or imputed", "Coral radial lollipop indicates BMI in kg/m²; an open grey point indicates missing BMI.", "Natural CRC patient-number order; no clinical-variable sorting",
  4, "inner-middle", "Sex", "sex", INPUT$clinical, "sex", "categorical", "Equal-width categorical ring", "Female; Male", paste(CLINICAL_COLOURS[c("Sex: Female", "Sex: Male")], collapse = ";"), "Unknown/Missing would be light grey", "Sex categories are shown by the first categorical ring.", "Natural CRC patient-number order; no clinical-variable sorting",
  5, "inner", "Tumor morphology", "tumor_form", INPUT$clinical, "tumor_form", "categorical", "Equal-width categorical ring", "Ulcerative; Protruding; Other/unspecified", paste(CLINICAL_COLOURS[c("Morphology: Ulcerative", "Morphology: Protruding", "Morphology: Other/unspecified")], collapse = ";"), "Other/unspecified is retained as a separate category", "Tumor morphology is derived only from the audited tumor_form field.", "Natural CRC patient-number order; no clinical-variable sorting",
  6, "innermost", "Tumor site", "site_group", INPUT$clinical, "site_group", "categorical", "Equal-width categorical ring", "Colon; Rectum; Unspecified", paste(CLINICAL_COLOURS[c("Site: Colon", "Site: Rectum", "Site: Unspecified")], collapse = ";"), "Unspecified is retained as a separate category", "Tumor site categories are shown by the innermost ring.", "Natural CRC patient-number order; no clinical-variable sorting"
)

fig1a_source <- patient_check %>%
  mutate(
    age_unit = "years",
    BMI_unit = "kg/m²",
    patient_order_rule = "Natural CRC patient-number order; no clinical-variable sorting"
  )
fig1b_source <- atlas %>%
  transmute(display_order, metric = display_label, value, unit, source_metric, note)
fig1c_source <- sqanti %>%
  transmute(display_order, structural_class, count = n_isoforms, percentage, denominator = 66805)
fig1d_source <- multiplicity %>%
  transmute(display_order, category = display_category, count = n_genes, percentage, denominator)
fig1e_source <- sample_detection %>%
  transmute(
    patient_id,
    sample_id,
    tissue_group,
    tissue_label,
    detected_isoforms,
    library_size,
    detection_definition = "Number of isoforms with raw count > 0 in the locked complete 75-sample matrix",
    pairing_key = patient_id
  ) %>%
  arrange(match(patient_id, expected_ids), match(tissue_group, c("N", "A", "T")))

path_a <- write_tsv_v3(fig1a_source, "Fig1A_patient_profile_v3.tsv")
path_b <- write_tsv_v3(fig1b_source, "Fig1B_atlas_overview_v3.tsv")
path_c <- write_tsv_v3(fig1c_source, "Fig1C_SQANTI3_composition_v3.tsv")
path_d <- write_tsv_v3(fig1d_source, "Fig1D_isoform_multiplicity_v3.tsv")
path_e <- write_tsv_v3(fig1e_source, "Fig1E_sample_detection_v3.tsv")
write_report_tsv(track_definition, "Fig1A_circular_track_definition_v3.tsv")
write_report_tsv(patient_check, "Fig1A_patient_profile_data_check_v3.tsv")

panel_a_long <- bind_rows(
  patient_profile %>% transmute(patient_id, patient_order, variable = "age", category = NA_character_, value = age, unit = "years"),
  patient_profile %>% transmute(patient_id, patient_order, variable = "BMI", category = NA_character_, value = BMI, unit = "kg/m²"),
  patient_profile %>% transmute(patient_id, patient_order, variable = "sex", category = sex, value = NA_real_, unit = NA_character_),
  patient_profile %>% transmute(patient_id, patient_order, variable = "tumor_morphology", category = tumor_morphology, value = NA_real_, unit = NA_character_),
  patient_profile %>% transmute(patient_id, patient_order, variable = "tumor_site", category = tumor_site, value = NA_real_, unit = NA_character_)
) %>%
  transmute(
    panel = "A",
    display_element = "patient-level profile",
    patient_id,
    sample_id = NA_character_,
    tissue_group = NA_character_,
    variable,
    category,
    value,
    unit,
    count = NA_real_,
    percentage = NA_real_,
    display_order = patient_order,
    notes = "One row per patient and displayed clinical variable"
  )

combined_source <- bind_rows(
  panel_a_long,
  fig1b_source %>%
    transmute(
      panel = "B", display_element = metric, patient_id = NA_character_, sample_id = NA_character_,
      tissue_group = NA_character_, variable = "atlas_overview", category = metric, value,
      unit, count = value, percentage = NA_real_, display_order, notes = note
    ),
  fig1c_source %>%
    transmute(
      panel = "C", display_element = "SQANTI3 composition", patient_id = NA_character_, sample_id = NA_character_,
      tissue_group = NA_character_, variable = "structural_class", category = structural_class,
      value = percentage, unit = "catalog fraction", count, percentage, display_order,
      notes = paste0("Denominator = ", denominator, " catalog isoforms")
    ),
  fig1d_source %>%
    transmute(
      panel = "D", display_element = "isoform multiplicity", patient_id = NA_character_, sample_id = NA_character_,
      tissue_group = NA_character_, variable = "isoforms_per_gene_bin", category,
      value = percentage, unit = "fraction of represented genes", count, percentage, display_order,
      notes = paste0("Denominator = ", denominator, " represented genes")
    ),
  fig1e_source %>%
    transmute(
      panel = "E", display_element = "sample-level detection", patient_id, sample_id,
      tissue_group, variable = "detected_isoforms", category = tissue_label,
      value = detected_isoforms, unit = "isoforms per sample", count = detected_isoforms,
      percentage = NA_real_, display_order = match(tissue_group, c("N", "A", "T")),
      notes = detection_definition
    )
)
combined_path <- write_tsv_v3(combined_source, "Fig1_source_data_v3.tsv")

source_mapping <- tribble(
  ~panel, ~display_element, ~source_file, ~source_sheet_or_table, ~source_columns, ~filter_or_definition, ~aggregation, ~n_observations, ~missing_value_handling, ~plot_type, ~notes,
  "A", "Patient ID and patient order", INPUT$age_bmi, "de-identified audited age/BMI freeze", "patient_id; patient_number", "P01-P25 display labels mapped one-to-one to CRC-1 through CRC-25", "None; natural numeric order", 25, "No missing IDs", "open circular sectors", "No ordering by clinical variables",
  "A", "Age", INPUT$age_bmi, "de-identified audited age/BMI freeze", "age_years; age_status", "available and 18-100 years", "Patient-level values; median/IQR summary", 25, "No missing age", "boxplot + points and radial lollipop", "Actual years; no 0-1 normalization",
  "A", "BMI", INPUT$age_bmi, "de-identified audited age/BMI freeze", "BMI_final; BMI_status", "available_recomputed", "Patient-level values; median/IQR summary", 24, "CRC-8 retained as NA; open grey point in circular track", "boxplot + points and radial lollipop", "Actual kg/m²; never imputed or plotted as zero",
  "A", "Sex", INPUT$clinical, "clinical_metadata_deidentified_v2.tsv", "sex", "Frozen de-identified categories", "Count and percentage", 25, "Unknown/Missing would be separate light grey", "100% stacked bar + categorical ring", "Female 12; Male 13",
  "A", "Tumor morphology", INPUT$clinical, "clinical_metadata_deidentified_v2.tsv", "tumor_form", "Ulcerative if source contains ulcerative term; Protruding if source contains protruding term; otherwise Other/unspecified", "Count and percentage", 25, "Other/unspecified retained", "100% stacked bar + categorical ring", "This is morphology, not differentiation",
  "A", "Tumor site", INPUT$clinical, "clinical_metadata_deidentified_v2.tsv", "site_group", "Colon; Rectum; Colorectum_unspecified to Unspecified", "Count and percentage", 25, "Unspecified retained", "100% stacked bar + categorical ring", "No left/right inference",
  "A", "N/A/T sample completeness", INPUT$sample_metadata_75, "metadata_sensitivity75.tsv", "sample_id; patient_id; tissue_group", "Exactly one N, A and T per patient", "Patient-level wide audit", 75, "Any missing tissue would fail-fast", "paired-study design annotation", "25 complete triplets",
  "B", "Atlas overview", INPUT$atlas, "Fig1C_v2_isoform_catalog_scale_source_data.tsv", "metric; value; source_metric", "Frozen formal metrics", "None", 3, "Missing metric fails", "three equal-width numeric modules", "Genes and isoforms are not compared by bar length",
  "C", "SQANTI3 structural composition", INPUT$sqanti3, "Fig1D_SQANTI3_structural_composition_source_data.tsv", "structural_class; n_isoforms; percent", "FSM/ISM/NIC/NNC/Other in fixed order", "Catalog composition", 5, "Missing category fails", "100% stacked horizontal composition bar", "Denominator 66,805 catalog isoforms",
  "D", "Isoform multiplicity", INPUT$multiplicity, "Fig1E_isoforms_per_gene_distribution_source_data.tsv", "diversity_bin; n_genes; percent_genes", "Four frozen bins", "Percentage of represented genes", 4, "Missing category fails", "vertical bar chart", "Denominator 16,530 represented genes",
  "E", "Detected isoforms per paired sample", INPUT$sample_metrics_75, "sensitivity75_DTE_preQC_sample_metrics.tsv", "sample_id; patient_id; tissue_group; detected_isoforms", "detected_isoforms = number of raw-count-positive isoforms per sample", "None; descriptive sample-level values", 75, "No samples removed", "paired lines + points + boxplots", "Pairing key is patient_id; no inferential test"
)
mapping_path <- write_report_tsv(source_mapping, "Fig1_source_data_mapping_v3.tsv")

# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------

panel_shell <- function(plot, tag, title, title_size = 9.5, title_x = 0.075) {
  cowplot::ggdraw(xlim = c(0, 1), ylim = c(0, 1), clip = "on") +
    cowplot::draw_label(
      tag, x = 0.005, y = 0.995, hjust = 0, vjust = 1,
      size = 13, fontface = "bold", fontfamily = FONT_FAMILY, colour = "#111111"
    ) +
    cowplot::draw_label(
      title, x = title_x, y = 0.987, hjust = 0, vjust = 1,
      size = title_size, fontface = "bold", fontfamily = FONT_FAMILY,
      colour = "#171717", lineheight = 0.95
    ) +
    cowplot::draw_plot(plot, x = 0.0, y = 0.0, width = 1.0, height = 0.88)
}

make_study_design <- function() {
  design_nodes <- tibble(
    x = 1:3,
    tissue = factor(c("Normal", "Adjacent", "Tumor"), levels = c("Normal", "Adjacent", "Tumor")),
    label = c("Normal\nn=25", "Adjacent\nn=25", "Tumor\nn=25")
  )
  ggplot(design_nodes, aes(x = x, y = 0.28)) +
    annotate("segment", x = 1, xend = 3, y = 0.28, yend = 0.28, linewidth = 0.38, colour = "#A8A8A8") +
    geom_point(aes(colour = tissue), size = 3.8) +
    geom_text(aes(y = 0.04, label = label), size = pt_size(7.3), lineheight = 0.90, colour = "#222222") +
    annotate(
      "text", x = 2, y = 0.95,
      label = "25 patients with colorectal cancer",
      size = pt_size(8.4), fontface = "bold", family = FONT_FAMILY, colour = "#202020"
    ) +
    annotate(
      "text", x = 2, y = 0.69,
      label = "25 complete N/A/T triplets  |  75 long-read transcriptomes",
      size = pt_size(7.5), family = FONT_FAMILY, colour = "#4D4D4D"
    ) +
    scale_colour_manual(values = TISSUE_COLOURS, guide = "none") +
    coord_cartesian(xlim = c(0.55, 3.45), ylim = c(-0.15, 1.03), clip = "off") +
    theme_void(base_family = FONT_FAMILY) +
    theme(plot.margin = margin(0, 2, 0, 2))
}

make_continuous_summary <- function(data, value_col, title, unit, colour, limits, breaks) {
  values <- data[[value_col]]
  stats <- tibble(
    n = sum(is.finite(values)),
    q1 = quantile(values, 0.25, na.rm = TRUE, names = FALSE),
    median = median(values, na.rm = TRUE),
    q3 = quantile(values, 0.75, na.rm = TRUE, names = FALSE)
  )
  plot_data <- data %>% filter(is.finite(.data[[value_col]])) %>%
    mutate(y_offset = seq(-0.085, 0.085, length.out = n()))
  ggplot(plot_data, aes(x = .data[[value_col]], y = 1 + y_offset)) +
    geom_boxplot(
      data = plot_data,
      mapping = aes(x = .data[[value_col]], y = 1, group = 1),
      inherit.aes = FALSE,
      orientation = "y", width = 0.34, outlier.shape = NA,
      fill = "white", colour = "#4A4A4A", linewidth = 0.34
    ) +
    geom_point(size = 1.25, alpha = 0.82, colour = colour) +
    scale_x_continuous(limits = limits, breaks = breaks, expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(breaks = NULL, limits = c(0.70, 1.30)) +
    labs(
      title = paste0(title, " (", unit, ")"),
      subtitle = paste0(
        "n=", stats$n, "; median ", format(round(stats$median, 1), nsmall = 1),
        "\nQ1-Q3 ", format(round(stats$q1, 1), nsmall = 1), "-",
        format(round(stats$q3, 1), nsmall = 1)
      ),
      x = NULL, y = NULL
    ) +
    theme_fig1(7.4) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 8.0, face = "bold", hjust = 0, lineheight = 0.88),
      plot.subtitle = element_text(size = 7.0, hjust = 0, colour = "#555555", margin = margin(b = 1), lineheight = 0.92),
      plot.margin = margin(1, 4, 2, 2)
    )
}

categorical_plot_data <- categorical_summary %>%
  mutate(
    fill_key = case_when(
      variable == "Sex" ~ paste0("Sex: ", category),
      variable == "Tumor morphology" ~ paste0("Morphology: ", category),
      variable == "Tumor site" ~ paste0("Site: ", category),
      TRUE ~ NA_character_
    ),
    bar_label = case_when(
      percentage < 0.10 ~ "",
      variable == "Sex" ~ paste0(category, "\nn=", n, " (", format_pct1(percentage), ")"),
      variable == "Tumor morphology" ~ paste0(category, "\nn=", n, " (", format_pct1(percentage), ")"),
      variable == "Tumor site" ~ paste0(category, "\nn=", n, " (", format_pct1(percentage), ")"),
      TRUE ~ ""
    ),
    text_colour = ifelse(fill_key == "Morphology: Protruding", "#2A2A2A", "white"),
    category_order = case_when(
      variable == "Sex" & category == "Female" ~ 1,
      variable == "Sex" & category == "Male" ~ 2,
      variable == "Tumor morphology" & category == "Ulcerative" ~ 1,
      variable == "Tumor morphology" & category == "Protruding" ~ 2,
      variable == "Tumor morphology" & category == "Other/unspecified" ~ 3,
      variable == "Tumor site" & category == "Colon" ~ 1,
      variable == "Tumor site" & category == "Rectum" ~ 2,
      variable == "Tumor site" & category == "Unspecified" ~ 3,
      TRUE ~ NA_real_
    ),
    variable_y = recode(variable, `Tumor site` = 1, `Tumor morphology` = 2, Sex = 3)
  ) %>%
  group_by(variable) %>%
  arrange(category_order, .by_group = TRUE) %>%
  mutate(
    xmin = lag(cumsum(percentage), default = 0),
    xmax = cumsum(percentage),
    xmid = (xmin + xmax) / 2,
    fill_key = factor(
      fill_key,
      levels = c(
        "Sex: Female", "Sex: Male",
        "Morphology: Ulcerative", "Morphology: Protruding", "Morphology: Other/unspecified",
        "Site: Colon", "Site: Rectum", "Site: Unspecified"
      )
    )
  ) %>%
  ungroup()
if (any(is.na(categorical_plot_data$fill_key)) ||
    any(!categorical_plot_data$fill_key %in% names(CLINICAL_COLOURS))) {
  stop("Clinical category colour mapping is incomplete.")
}

make_categorical_summary <- function() {
  ggplot(categorical_plot_data) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = variable_y - 0.34, ymax = variable_y + 0.34, fill = fill_key),
      colour = "white", linewidth = 0.45
    ) +
    geom_text(
      aes(x = xmid, y = variable_y, label = bar_label, colour = text_colour),
      size = pt_size(7.0), fontface = "bold", lineheight = 0.88, show.legend = FALSE
    ) +
    scale_fill_manual(values = CLINICAL_COLOURS, guide = "none") +
    scale_colour_identity() +
    scale_x_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = 1:3,
      labels = c("Tumor site", "Tumor\nmorphology", "Sex"),
      limits = c(0.55, 3.45),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL, y = NULL,
      caption = "Grey 4.0% segments: morphology other/unspecified, n=1; site unspecified, n=1."
    ) +
    theme_fig1(7.4) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = 7.0, colour = "#555555"),
      axis.text.y = element_text(size = 7.2, face = "bold", lineheight = 0.9),
      plot.caption = element_text(size = 7.0, hjust = 0, colour = "#555555", margin = margin(t = 2), lineheight = 0.9),
      plot.margin = margin(2, 4, 2, 2)
    )
}

draw_open_patient_circle <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  on.exit(circlize::circos.clear(), add = TRUE)
  par(mar = c(0.1, 0.1, 0.1, 0.1), family = FONT_FAMILY, ps = 8.5, xpd = NA)
  circlize::circos.clear()

  sectors <- factor(patient_profile$display_patient_id, levels = patient_profile$display_patient_id)
  n_sector <- nrow(patient_profile)
  circlize::circos.par(
    start.degree = 135,
    clock.wise = TRUE,
    gap.after = c(rep(0.8, n_sector - 1), 86),
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.0025, 0.0025),
    points.overflow.warning = FALSE,
    canvas.xlim = c(-1.30, 1.30),
    canvas.ylim = c(-1.30, 1.30)
  )
  circlize::circos.initialize(factors = sectors, x = rep(1, n_sector), xlim = c(0.5, 1.5))

  circlize::circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.20, bg.border = NA,
    panel.fun = function(x, y) {
      sector_name <- circlize::get.cell.meta.data("sector.index")
      circlize::circos.segments(1, 0.02, 1, 0.20, col = "#B7B7B7", lwd = 0.45)
      circlize::circos.text(
        1, 0.63, sector_name,
        facing = "clockwise", niceFacing = TRUE,
        cex = 0.82, col = "#292929", family = FONT_FAMILY
      )
    }
  )

  draw_continuous_track <- function(values, limits, ticks, colour, missing_open = FALSE) {
    circlize::circos.trackPlotRegion(
      ylim = limits, track.height = 0.15, bg.border = "#D0D0D0", bg.lwd = 0.45,
      panel.fun = function(x, y) {
        i <- circlize::get.cell.meta.data("sector.numeric.index")
        value <- values[i]
        if (is.finite(value)) {
          circlize::circos.segments(1, limits[1], 1, value, col = colour, lwd = 0.85)
          circlize::circos.points(1, value, pch = 16, cex = 0.52, col = colour)
        } else if (missing_open) {
          circlize::circos.points(
            1, limits[1] + 0.10 * diff(limits),
            pch = 1, cex = 0.64, col = "#7F7F7F"
          )
        }
      }
    )
  }

  draw_categorical_track <- function(values, prefix) {
    keys <- paste0(prefix, values)
    colours <- unname(CLINICAL_COLOURS[keys])
    colours[is.na(colours)] <- MISSING_COLOUR
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = 0.082, bg.border = NA,
      panel.fun = function(x, y) {
        i <- circlize::get.cell.meta.data("sector.numeric.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        circlize::circos.rect(
          xleft = xlim[1], ybottom = 0,
          xright = xlim[2], ytop = 1,
          col = colours[i], border = "white", lwd = 0.45
        )
      }
    )
  }

  draw_continuous_track(patient_profile$age, c(30, 90), c(30, 60, 90), AGE_COLOUR)
  draw_continuous_track(patient_profile$BMI, c(15, 30), c(15, 22.5, 30), BMI_COLOUR, missing_open = TRUE)
  draw_categorical_track(patient_profile$sex, "Sex: ")
  draw_categorical_track(patient_profile$tumor_morphology, "Morphology: ")
  draw_categorical_track(patient_profile$tumor_site, "Site: ")

}

p_age <- make_continuous_summary(patient_profile, "age", "Age", "years", AGE_COLOUR, c(30, 90), c(30, 60, 90))
p_bmi <- make_continuous_summary(patient_profile, "BMI", "BMI", "kg/m²", BMI_COLOUR, c(15, 30), c(15, 22.5, 30))
p_cat <- make_categorical_summary()
p_summary <- cowplot::plot_grid(
  p_age,
  p_bmi,
  p_cat,
  ncol = 1,
  rel_heights = c(0.28, 0.28, 0.44),
  align = "v"
)
p_circle <- muffle_known_tnr_metric_warning(
  ggplotify::as.ggplot(~draw_open_patient_circle())
) +
  theme(plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1)
p_design <- make_study_design()

p_a_body <- cowplot::ggdraw(xlim = c(0, 1), ylim = c(0, 1), clip = "on") +
  cowplot::draw_label(
    "Cohort summary", x = 0.012, y = 0.955, hjust = 0, vjust = 1,
    size = 8.2, fontface = "bold", fontfamily = FONT_FAMILY, colour = "#262626"
  ) +
  cowplot::draw_label(
    "Open circular patient-level profile", x = 0.415, y = 0.955, hjust = 0, vjust = 1,
    size = 8.2, fontface = "bold", fontfamily = FONT_FAMILY, colour = "#262626"
  ) +
  cowplot::draw_plot(p_summary, x = 0.005, y = 0.105, width = 0.390, height = 0.81) +
  cowplot::draw_plot(p_circle, x = 0.385, y = 0.075, width = 0.610, height = 0.86) +
  cowplot::draw_label(
    "Age\n(years)\n30 | 60 | 90", x = 0.405, y = 0.690, hjust = 0, vjust = 0.5,
    size = 7.0, fontfamily = FONT_FAMILY, colour = AGE_COLOUR, lineheight = 0.90
  ) +
  cowplot::draw_label(
    "BMI\n(kg/m²)\n15 | 22.5 | 30", x = 0.405, y = 0.505, hjust = 0, vjust = 0.5,
    size = 7.0, fontfamily = FONT_FAMILY, colour = BMI_COLOUR, lineheight = 0.90
  ) +
  cowplot::draw_label(
    "Sex", x = 0.425, y = 0.390, hjust = 0, vjust = 0.5,
    size = 7.0, fontfamily = FONT_FAMILY, colour = "#444444"
  ) +
  cowplot::draw_label(
    "Morphology", x = 0.425, y = 0.345, hjust = 0, vjust = 0.5,
    size = 7.0, fontfamily = FONT_FAMILY, colour = "#444444"
  ) +
  cowplot::draw_label(
    "Tumor site", x = 0.425, y = 0.300, hjust = 0, vjust = 0.5,
    size = 7.0, fontfamily = FONT_FAMILY, colour = "#444444"
  ) +
  cowplot::draw_label(
    "P01-P25: de-identified patients  |  tracks outer to inner: ID, age, BMI, sex, morphology, site",
    x = 0.012, y = 0.074, hjust = 0, vjust = 1,
    size = 7.0, fontfamily = FONT_FAMILY, colour = "#4A4A4A"
  ) +
  cowplot::draw_label(
    "Purple: age  |  coral: BMI  |  open grey point: missing  |  ring colours match cohort summary",
    x = 0.012, y = 0.025, hjust = 0, vjust = 1,
    size = 7.0, fontfamily = FONT_FAMILY, colour = "#5A5A5A"
  )

p_a <- cowplot::plot_grid(
  p_design,
  p_a_body,
  ncol = 1,
  rel_heights = c(0.17, 0.83),
  align = "v"
)

make_atlas_overview <- function() {
  xs <- c(0.17, 0.50, 0.83)
  p <- cowplot::ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
    cowplot::draw_line(x = c(0.335, 0.335), y = c(0.19, 0.82), colour = "#D4D4D4", size = 0.35) +
    cowplot::draw_line(x = c(0.665, 0.665), y = c(0.19, 0.82), colour = "#D4D4D4", size = 0.35)
  for (i in seq_len(nrow(atlas))) {
    p <- p +
      cowplot::draw_label(
        format_n(atlas$value[i]), x = xs[i], y = 0.62,
        size = 12.5, fontface = "bold", fontfamily = FONT_FAMILY, colour = "#42576A"
      ) +
      cowplot::draw_label(
        sub(" ", "\n", atlas$display_label[i]), x = xs[i], y = 0.31,
        size = 7.6, fontfamily = FONT_FAMILY, colour = "#303030", lineheight = 0.92
      )
  }
  p
}

p_c_bar <- ggplot(sqanti) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0.33, ymax = 0.72, fill = structural_class),
    colour = "white", linewidth = 0.4
  ) +
  annotate(
    "text", x = 0, y = 0.91, hjust = 0,
    label = "Catalog isoforms, n=66,805",
    size = pt_size(7.5), fontface = "bold", family = FONT_FAMILY, colour = "#333333"
  ) +
  scale_fill_manual(values = SQANTI_COLOURS, guide = "none") +
  scale_x_continuous(labels = label_percent(accuracy = 1), breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.12, 1.0), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(
    axis.text.x = element_text(size = 7.0, colour = "#555555"),
    axis.ticks.x = element_line(linewidth = 0.3, colour = "#555555"),
    plot.margin = margin(1, 4, 0, 4)
  )

p_c_key <- ggplot(sqanti, aes(y = reorder(structural_class, -display_order))) +
  geom_point(aes(x = 0, colour = structural_class), shape = 15, size = 3.0) +
  geom_text(
    aes(x = 0.06, label = display_label), hjust = 0,
    size = pt_size(7.0), colour = "#2F2F2F"
  ) +
  scale_colour_manual(values = SQANTI_COLOURS, guide = "none") +
  coord_cartesian(xlim = c(-0.01, 0.74), clip = "off") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.margin = margin(0, 4, 1, 4))

p_c <- cowplot::plot_grid(p_c_bar, p_c_key, ncol = 1, rel_heights = c(0.40, 0.60), align = "v")

p_d <- ggplot(
  multiplicity,
  aes(x = factor(display_category, levels = display_category), y = percentage)
) +
  geom_col(width = 0.62, fill = MULTIPLICITY_COLOUR) +
  geom_text(
    aes(label = paste0(format_pct1(percentage), "\nn=", format_n(n_genes))),
    vjust = -0.28, size = pt_size(7.1), lineheight = 0.92, colour = "#252525"
  ) +
  scale_x_discrete(labels = c("1\nisoform", "2–3\nisoforms", "4–5\nisoforms", "≥6\nisoforms")) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    breaks = c(0, 0.1, 0.2, 0.3),
    limits = c(0, 0.37),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Percentage of represented genes") +
  theme_fig1(8.0) +
  theme(
    axis.text.x = element_text(size = 7.8, lineheight = 0.92),
    panel.grid.major.y = element_line(linewidth = 0.22, colour = "#E6E6E6"),
    axis.line.x = element_line(linewidth = 0.32),
    plot.margin = margin(3, 5, 3, 5)
  )

sample_plot <- sample_detection %>%
  mutate(
    patient_number = as.integer(sub("^CRC-", "", patient_id)),
    tissue_label = factor(tissue_label, levels = c("Normal", "Adjacent", "Tumor")),
    tissue_x = as.numeric(tissue_label),
    patient_offset = scales::rescale(patient_number, to = c(-0.075, 0.075)),
    x_position = tissue_x + patient_offset
  )

p_e <- ggplot(sample_plot, aes(y = detected_isoforms)) +
  geom_line(
    aes(x = x_position, group = patient_id),
    colour = "#B5B5B5", linewidth = 0.34, alpha = 0.48
  ) +
  geom_boxplot(
    aes(x = tissue_x, group = tissue_label, colour = tissue_label, fill = tissue_label),
    width = 0.45, alpha = 0.13, outlier.shape = NA, linewidth = 0.45
  ) +
  geom_point(
    aes(x = x_position, colour = tissue_label),
    size = 1.35, alpha = 0.80
  ) +
  annotate(
    "text", x = 1:3, y = max(sample_plot$detected_isoforms) + 650,
    label = "n=25", size = pt_size(7.0), family = FONT_FAMILY, colour = "#4A4A4A"
  ) +
  scale_colour_manual(values = TISSUE_COLOURS, guide = "none") +
  scale_fill_manual(values = TISSUE_COLOURS, guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = c("Normal", "Adjacent", "Tumor"), limits = c(0.62, 3.38)) +
  scale_y_continuous(
    labels = label_comma(),
    breaks = seq(35000, 55000, by = 5000),
    limits = c(min(sample_plot$detected_isoforms) - 900, max(sample_plot$detected_isoforms) + 1500),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Detected isoforms per sample") +
  theme_fig1(8.0) +
  theme(
    panel.grid.major.y = element_line(linewidth = 0.22, colour = "#E6E6E6"),
    axis.line.x = element_line(linewidth = 0.32),
    axis.text.x = element_text(size = 8.0),
    plot.margin = margin(3, 5, 3, 5)
  )

p_a_shell <- panel_shell(p_a, "A", "Paired study design and patient-level cohort profile", title_size = 9.5, title_x = 0.055)
p_b_shell <- panel_shell(make_atlas_overview(), "B", "Atlas overview", title_size = 9.5, title_x = 0.16)
p_c_shell <- panel_shell(p_c, "C", "SQANTI3 structural composition\nof catalog isoforms", title_size = 8.2, title_x = 0.16)
p_d_shell <- panel_shell(p_d, "D", "Isoform multiplicity across\nrepresented genes", title_size = 8.8, title_x = 0.10)
p_e_shell <- panel_shell(p_e, "E", "Sample-level isoform detection across paired tissues", title_size = 9.3, title_x = 0.075)

fig1 <- cowplot::ggdraw(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  cowplot::draw_plot(p_a_shell, x = 0.012, y = 0.425, width = 0.675, height = 0.565) +
  cowplot::draw_plot(p_b_shell, x = 0.705, y = 0.720, width = 0.283, height = 0.270) +
  cowplot::draw_plot(p_c_shell, x = 0.705, y = 0.425, width = 0.283, height = 0.280) +
  cowplot::draw_plot(p_d_shell, x = 0.012, y = 0.018, width = 0.405, height = 0.385) +
  cowplot::draw_plot(p_e_shell, x = 0.440, y = 0.018, width = 0.548, height = 0.385)

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------

open_svg <- function(path, width_in, height_in) {
  svglite::svglite(path, width = width_in, height = height_in, bg = "white")
}

save_vector_and_png <- function(plot, basename, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  pdf_path <- file.path(DIR_PDF, paste0(basename, ".pdf"))
  svg_path <- file.path(DIR_SVG, paste0(basename, ".svg"))
  png_path <- file.path(DIR_PNG, paste0(basename, ".png"))

  Cairo::CairoPDF(pdf_path, width = width_in, height = height_in, family = FONT_FAMILY, bg = "white")
  muffle_known_tnr_metric_warning(print(plot))
  grDevices::dev.off()

  open_svg(svg_path, width_in, height_in)
  muffle_known_tnr_metric_warning(print(plot))
  grDevices::dev.off()

  ragg::agg_png(
    png_path, width = width_in, height = height_in, units = "in",
    res = dpi, background = "white", scaling = 1
  )
  muffle_known_tnr_metric_warning(print(plot))
  grDevices::dev.off()

  tibble(format = c("pdf", "svg", "png"), path = c(pdf_path, svg_path, png_path))
}

save_preview <- function(plot, filename, width_mm, height_mm, dpi = 600) {
  path <- file.path(DIR_PNG, filename)
  ragg::agg_png(
    path,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    units = "in",
    res = dpi,
    background = "white",
    scaling = 1
  )
  muffle_known_tnr_metric_warning(print(plot))
  grDevices::dev.off()
  path
}

fig1_files <- save_vector_and_png(
  fig1,
  "Fig1_paired_longread_isoform_atlas_v3",
  width_mm = 180,
  height_mm = 170,
  dpi = 600
)
fig1a_files <- save_vector_and_png(
  p_a_shell,
  "Fig1A_patient_profile_open_circular_v3",
  width_mm = 180,
  height_mm = 115,
  dpi = 600
)
fig1_preview <- save_preview(
  fig1,
  "Fig1_paired_longread_isoform_atlas_v3_180mm_preview.png",
  width_mm = 180,
  height_mm = 170,
  dpi = 600
)
fig1a_preview <- save_preview(
  p_a_shell,
  "Fig1A_patient_profile_open_circular_v3_180mm_preview.png",
  width_mm = 180,
  height_mm = 115,
  dpi = 600
)

all_outputs <- c(fig1_files$path, fig1a_files$path, fig1_preview, fig1a_preview,
                 path_a, path_b, path_c, path_d, path_e, combined_path, mapping_path)
if (any(!file.exists(all_outputs)) || any(file.info(all_outputs)$size <= 0)) {
  stop("One or more Figure 1 v3 outputs are missing or empty.")
}

# -----------------------------------------------------------------------------
# Figure 2-5 integrity guard
# -----------------------------------------------------------------------------

baseline <- read_tsv_char(INPUT$mtime_baseline) %>%
  mutate(
    size_bytes = as_num(size_bytes),
    mtime_epoch = as_num(mtime_epoch),
    basename = basename(path)
  )
guard <- baseline %>%
  filter(grepl("^(Fig|Figure)[2-5]", basename)) %>%
  rowwise() %>%
  mutate(
    current_exists = file.exists(path),
    current_size_bytes = ifelse(current_exists, file.info(path)$size, NA_real_),
    current_mtime_epoch = ifelse(current_exists, as.numeric(file.info(path)$mtime), NA_real_),
    current_sha256 = ifelse(current_exists, sha256_file(path), NA_character_),
    unchanged = current_exists &
      current_size_bytes == size_bytes &
      floor(current_mtime_epoch) == mtime_epoch &
      current_sha256 == sha256
  ) %>%
  ungroup() %>%
  select(scope, path, baseline_size_bytes = size_bytes, current_size_bytes,
         baseline_mtime_epoch = mtime_epoch, current_mtime_epoch,
         baseline_sha256 = sha256, current_sha256, unchanged)
guard_path <- write_report_tsv(guard, "Fig1_Fig2_to_Fig5_mtime_guard_v3.tsv")
if (nrow(guard) == 0 || any(!guard$unchanged)) {
  stop("Figure 2-5 mtime/SHA256 guard failed. See ", guard_path)
}

# -----------------------------------------------------------------------------
# Legend and reports
# -----------------------------------------------------------------------------

age_s <- continuous_summary %>% filter(variable == "Age")
bmi_s <- continuous_summary %>% filter(variable == "BMI")
legend_lines <- c(
  "Figure 1. A paired long-read isoform atlas across the normal-adjacent-tumor continuum in colorectal cancer.",
  "",
  paste0(
    "(A) Study design and patient-level clinical profile of 25 patients with complete matched normal (N), adjacent (A) and tumor (T) tissues, yielding 75 long-read transcriptomes. ",
    "The cohort summary shows age, body mass index (BMI), sex, tumor morphology and tumor site. Age and BMI points represent individual patients; boxes show the median and interquartile range, with whiskers extending to 1.5 times the interquartile range. ",
    "Age was available for all 25 patients (median ", format(round(age_s$median, 1), nsmall = 1), " years; Q1-Q3 ",
    format(round(age_s$q1, 1), nsmall = 1), "-", format(round(age_s$q3, 1), nsmall = 1), "), and BMI was available for 24 patients (median ",
    format(round(bmi_s$median, 1), nsmall = 1), " kg/m²; Q1-Q3 ", format(round(bmi_s$q1, 1), nsmall = 1), "-",
    format(round(bmi_s$q3, 1), nsmall = 1), "). Each sector in the open circular profile represents one de-identified patient (P01 to P25) in natural patient-number order; the source-data table preserves the one-to-one mapping to CRC-1 through CRC-25. ",
    "Tracks from outer to inner show patient ID, age, BMI, sex, tumor morphology and tumor site. Radial lines extend from the displayed scale minimum to the observed age or BMI, and filled points mark the observed values. An open grey point denotes the missing BMI value for CRC-8; missing values were not imputed or plotted as zero. Grey denotes missing or unspecified information."
  ),
  "",
  paste0(
    "(B) Overview of the long-read isoform atlas: ", format_n(atlas$value[1]), " catalog isoforms, ",
    format_n(atlas$value[2]), " quantified isoforms and ", format_n(atlas$value[3]), " represented genes. ",
    "Catalog isoforms are the final SQANTI3-annotated master catalog. Quantified isoforms are the 63,623 isoforms retained in the locked primary expression-analysis matrix under the inherited support rule (raw count at least 2 in at least 8 samples and total count at least 20). Represented genes are unique gene identifiers represented in the master catalog."
  ),
  "",
  paste0(
    "(C) SQANTI3 structural classification of catalog isoforms. Counts and percentages are calculated relative to the full catalog of ",
    format_n(sum(sqanti$n_isoforms)), " isoforms. FSM, full-splice match; ISM, incomplete-splice match; NIC, novel in catalog; NNC, novel not in catalog."
  ),
  "",
  paste0(
    "(D) Distribution of isoform multiplicity across ", format_n(sum(multiplicity$n_genes)),
    " represented genes. Bars show the percentage of represented genes; labels provide percentages and gene counts."
  ),
  "",
  paste0(
    "(E) Sample-level numbers of detected isoforms across paired normal, adjacent and tumor tissues (n=25 samples per tissue). ",
    "Each point represents one sample, and lines connect samples from the same patient using patient ID as the pairing key. Boxes show the median and interquartile range, with whiskers extending to 1.5 times the interquartile range. ",
    "Detected isoforms are defined as isoforms with a raw count greater than zero in each sample of the locked complete 75-sample matrix. This panel is a descriptive dataset overview: detected numbers may reflect sequencing depth, sample quality and transcript complexity and are not interpreted as evidence of tissue-specific biological remodeling."
  ),
  "",
  "Normal denotes histologically normal tissue, adjacent denotes paired tissue adjacent to the tumor, and tumor denotes the matched colorectal tumor. No inferential test was added for this figure."
)
legend_path <- file.path(DIR_REPORTS, "Figure1_legend_draft_v3.txt")
writeLines(legend_lines, legend_path, useBytes = TRUE)

output_info <- tibble(path = c(fig1_files$path, fig1a_files$path, fig1_preview, fig1a_preview)) %>%
  mutate(
    size_bytes = file.info(path)$size,
    mtime = format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S %z"),
    sha256 = vapply(path, sha256_file, character(1))
  )

redesign_lines <- c(
  "# Figure 1 v3 redesign summary",
  "",
  "## Execution boundary",
  "",
  "- A new Figure 1-only render script was created; the existing combined Figure 1-3 script was not modified or executed.",
  "- No FASTQ preprocessing, Pychopper, minimap2, FLAIR, SQANTI3, DTE/DTU, SUPPA2, ORF/NMD or Pfam analysis was rerun.",
  "- No raw data or frozen statistical result was modified.",
  "- Figure 2-5 files and source data passed the pre/post size, mtime and SHA256 guard.",
  "",
  "## Script and inputs",
  "",
  paste0("- Render script: `", file.path(BASE, "scripts", "10_build_Fig1_paired_longread_isoform_atlas_v3.R"), "`."),
  paste0("- Panel A: `", INPUT$age_bmi, "`, `", INPUT$clinical, "` and `", INPUT$sample_metadata_75, "`."),
  paste0("- Panel B: `", INPUT$atlas, "`."),
  paste0("- Panel C: `", INPUT$sqanti3, "`."),
  paste0("- Panel D: `", INPUT$multiplicity, "`."),
  paste0("- Panel E: `", INPUT$sample_metrics_75, "` joined to `", INPUT$sample_metadata_75, "` by sample ID."),
  "",
  "## Panel A: open circular patient profile",
  "",
  "- The former small closed circle was replaced by a 25-sector open circular layout with an 86-degree left-side gap.",
  "- The upper strip explicitly shows 25 patients, complete matched Normal/Adjacent/Tumor collection, n=25 per tissue and 75 long-read transcriptomes.",
  "- The left column reports cohort-level age, BMI, sex, tumor morphology and tumor site using the same values and colors as the circular profile.",
  "- Tracks from outer to inner are: patient ID, age, BMI, sex, tumor morphology and tumor site.",
  "- Compact display IDs P01-P25 are mapped one-to-one to CRC-1 through CRC-25 and ordered by natural patient number; no clinical-variable sorting was used.",
  paste0("- Age: n=", age_s$n, ", median ", round(age_s$median, 1), ", Q1 ", round(age_s$q1, 1), ", Q3 ", round(age_s$q3, 1), ", range ", round(age_s$minimum, 1), "-", round(age_s$maximum, 1), " years."),
  paste0("- BMI: n=", bmi_s$n, ", median ", round(bmi_s$median, 2), ", Q1 ", round(bmi_s$q1, 2), ", Q3 ", round(bmi_s$q3, 2), ", range ", round(bmi_s$minimum, 2), "-", round(bmi_s$maximum, 2), " kg/m²."),
  "- CRC-8 is the only missing BMI value; it remains NA and is shown as an open grey point, never as zero or an imputed value.",
  "- The audited CRC-7 BMI recomputation (26.638918 kg/m²) is retained.",
  "- The source field tumor_form is displayed as tumor morphology, not differentiation.",
  "- Sex: Female 12 (48.0%), Male 13 (52.0%).",
  "- Tumor morphology: Ulcerative 16 (64.0%), Protruding 8 (32.0%), Other/unspecified 1 (4.0%).",
  "- Tumor site: Colon 12 (48.0%), Rectum 12 (48.0%), Unspecified 1 (4.0%). No left/right site was inferred.",
  "",
  "## Panels B-E",
  "",
  "- Panel B uses three equal-width numerical modules because gene and isoform units should not be compared by a shared bar length.",
  paste0("- Formal atlas values are ", format_n(atlas$value[1]), " catalog isoforms, ", format_n(atlas$value[2]), " quantified isoforms and ", format_n(atlas$value[3]), " represented genes."),
  paste0("- Panel C uses one 100% composition bar: ", paste(sqanti$display_label, collapse = "; "), "."),
  paste0("- Panel D uses percentage of represented genes with denominator ", format_n(sum(multiplicity$n_genes)), ": ", paste0(multiplicity$display_category, " ", format_pct1(multiplicity$percentage), " (n=", format_n(multiplicity$n_genes), ")", collapse = "; "), "."),
  "- Panel E uses the existing official complete 75-sample metrics table. Lines use patient_id, not row position, and all 25 N/A/T triplets are retained.",
  "- Panel E remains descriptive and contains no p values, significance stars or biological-remodeling interpretation.",
  "",
  "## Numeric closure and legacy issues",
  "",
  "- SQANTI3 counts close to 66,805 and percentages close to 100.0%.",
  "- Isoform-multiplicity counts close to 16,530 and percentages close to 100.0%.",
  "- All clinical categorical distributions close to 25 patients and 100.0%.",
  "- The old visual could leave the 4% tumor-site category implicit; v3 explicitly labels Unspecified 1 (4.0%).",
  "- The previous Figure 1 sample-detection table contained 72 samples/24 patients. V3 uses the existing official sensitivity75 pre-QC metrics to display all 25 complete pairs without recomputation.",
  "- The normalized clinical metadata contains an invalid CRC-1 age entry; the already audited age/BMI freeze remains authoritative and records age 50.",
  "",
  "## Remaining limitation",
  "",
  "- The 63,623 quantified-isoform atlas metric belongs to the locked primary expression-analysis matrix, whereas Panel E intentionally uses the existing complete 75-sample descriptive metrics to retain all 25 matched triplets. This cohort-boundary distinction is stated in the legend and must not be interpreted as a new analysis.",
  "",
  "## Outputs",
  "",
  paste0("- `", output_info$path, "` (", format_n(output_info$size_bytes), " bytes; SHA256 ", output_info$sha256, ")"),
  paste0("- Combined source data: `", combined_path, "`."),
  paste0("- Source-data mapping: `", mapping_path, "`."),
  paste0("- Legend draft: `", legend_path, "`."),
  paste0("- Figure 2-5 guard: `", guard_path, "`.")
)
redesign_path <- file.path(DIR_REPORTS, "Fig1_redesign_summary_v3.md")
writeLines(redesign_lines, redesign_path, useBytes = TRUE)

qc_lines <- c(
  "# Figure 1 v3 QC report",
  "",
  paste0("Render date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "",
  "## Data and scope QA",
  "",
  "- PASS: Figure 1 is restricted to the paired cohort and long-read isoform-atlas overview.",
  "- PASS: No DTE, DTU, functional consequence, clinical association, recurrence, prognosis or representative-gene result is shown.",
  "- PASS: No upstream analysis was rerun and no existing Figure 1-5 file was overwritten.",
  "- PASS: 25 patients, 75 samples and exactly one N/A/T sample per patient.",
  "- PASS: patient labels include all P01-P25 in natural numeric order, with one-to-one CRC ID mapping preserved in source data.",
  "- PASS: age n=25 and BMI n=24; CRC-8 BMI remains NA and is not encoded as zero.",
  "- PASS: all clinical categorical summaries close to n=25 and 100.0%.",
  "- PASS: SQANTI3 closes to 66,805 and 100.0%; multiplicity closes to 16,530 and 100.0%.",
  "- PASS: Panel E uses patient_id for lines, retains all samples/outliers and adds no inferential test.",
  "",
  "## Rendering QA",
  "",
  paste0("- PASS: Times New Roman resolved to `", font_match$path[1], "`."),
  "- PASS: PDF and SVG were rendered through R vector devices.",
  "- PASS: PNG and 180 mm previews were rendered at 600 dpi with a white background.",
  "- PASS: full figure width is 180 mm; standalone Panel A width is 180 mm.",
  "- PENDING MANUAL VISUAL QA: clipping, patient-label readability, composition labels, Panel E pair lines and overall balance must be reviewed from the downloaded 180 mm PNG.",
  "",
  "## Figure 2-5 guard",
  "",
  paste0("- PASS: ", nrow(guard), " pre-existing Figure 2-5 figure/source-data files retained identical size, mtime and SHA256."),
  paste0("- Detailed guard: `", guard_path, "`."),
  "",
  "## Output inventory",
  "",
  paste0("- ", output_info$path, " | ", format_n(output_info$size_bytes), " bytes | ", output_info$mtime),
  paste0("- `", path_a, "`"),
  paste0("- `", path_b, "`"),
  paste0("- `", path_c, "`"),
  paste0("- `", path_d, "`"),
  paste0("- `", path_e, "`"),
  paste0("- `", combined_path, "`"),
  paste0("- `", mapping_path, "`"),
  paste0("- `", legend_path, "`")
)
qc_path <- file.path(DIR_REPORTS, "Fig1_QC_v3.md")
writeLines(qc_lines, qc_path, useBytes = TRUE)

message("Figure 1 v3 render completed.")
message("Main PNG: ", fig1_files$path[fig1_files$format == "png"])
message("Panel A PNG: ", fig1a_files$path[fig1a_files$format == "png"])
message("Figure 2-5 guard: PASS (", nrow(guard), " files)")
