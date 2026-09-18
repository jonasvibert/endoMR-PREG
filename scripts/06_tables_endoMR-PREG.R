#!/usr/bin/env Rscript

###############################################################################
# 06_tables_endoMR-PREG.R
#
# Generate main and supplementary tables for:
#   Genetic liability to endometriosis → pregnancy outcomes (endoMR-PREG)
###############################################################################
# REVISION LOG
#   [R3.2] Trio-based MR de-emphasis (Reviewer #3, "Trio-based analysis and
#          title framing"). Supplementary Table S4 ("trio-based MR estimates":
#          maternal/fetal/paternal DONUTS estimates) has been REMOVED from the
#          exported tables and replaced by two new tables built from the main
#          (non-trio) harmonised dataset:
#            - Table S4: heterogeneity statistics (Cochran's Q, IVW only)
#            - Table S5: MR-Egger intercept (directional pleiotropy)
#          Confirmed rationale: the trio table simply duplicated the content
#          of Figure 3 ("Trios overview", script 07) - removed as redundant,
#          not solely to reduce prominence per the reviewer's wording.
#   [STYLE] US spelling applied throughout outcome labels (haemorrhage ->
#          hemorrhage, anaemia -> anemia). Confirmed: author preference, for
#          consistency with the rest of the manuscript (not a journal
#          requirement).
#   [FIX]  McBride et al. (MR-PREG) publication year corrected 2026 -> 2025 in
#          Table 1 (table1_gwas).
#   [REVERTED] domain_meta: order_within swap for ga_all/rup_memb within the
#          "Pregnancy timing" domain was unintentional - reverted to original
#          order (ga_all before rup_memb).
#   [EDIT] Table numbering/naming reworked across the board (e.g. Supp_Table_S1
#          -> Supp_Table_1A/B/C, new Table_S1_harmonised_summary, Supp_Table_S2
#          -> Supplementary_Table_2_Endometriosis_Instruments, Supp_Table_S3
#          -> Table_S3_all_MR_methods) to align with final manuscript table
#          numbering.
#   [R2.1] Reviewer #2, comment 1 (colocalisation): added Supplementary Table
#          S6, the locus-by-locus coloc.abf results (41 loci) referenced in
#          the response letter, from results/coloc_endo_praevia_by_locus.csv
#          (script 05.3). Previously computed but not exported as a table.
#   [R2.5] Reviewer #2, comment 5 (adenomyosis misclassification): added
#          Supplementary Table S7, the exploratory adenomyosis -> placenta
#          praevia MR (IVW/MR-Egger/Weighted median, Cochran's Q, Egger
#          intercept), from results/koller_adenomyosis_*.csv (script 05.1).
#          Includes an explicit sample-overlap caveat, since (unlike the
#          primary endometriosis analysis) no MRlap/LDSC overlap check was
#          run for the Koller adenomyosis GWAS vs. FinnGen R12.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(here)
  library(data.table)
  library(tidyr)
})

source(here::here("config", "config.R"))
source(here::here("config", "utils.R"))

setup_logging("5_tables_MR_analysis")
log_info("=== Starting table generation for endoMR-PREG study ===")

load_packages(c("openxlsx", "knitr", "tibble", "purrr", "gt", "TwoSampleMR"))

results_dir <- PATHS$results
tables_dir  <- file.path(results_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

harm_file <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")
if (!file.exists(harm_file)) {
  log_error("Harmonised data not found: ", harm_file)
  stop("Harmonised data not found.")
}

###############################################################################
# 1) LOAD HARMONISED DATA & DEFINE 30 PRIMARY OUTCOMES
###############################################################################

dat <- data.table::fread(harm_file)
log_info("Loaded harmonised data with ", nrow(dat), " rows")

stopifnot(all(c("id.exposure", "beta.exposure", "beta.outcome") %in% colnames(dat)))
stopifnot("outcome" %in% colnames(dat))

# ---- 30 retained outcomes (aligned across scripts) ---------------------------

vars_keep <- c(
  # Placenta & bleeding (7)
  "Antepartum_bleeding",
  "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta",
  "finngen_R12_O15_PLAC_PRAEVIA",
  "finngen_R12_O15_PLAC_DISORD",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",
  
  # Membranes (1)
  "rup_memb",
  
  # Preterm birth (2)
  "pretb_all",
  "vpretb_all",
  
  # Growth / GA / weight (6)
  "ga_all",
  "sga",
  "lbw_all",
  "hbw_all",
  "lga",
  "zbw_all",
  
  # Neonatal / Apgar / perinatal death (4)
  "lowapgar1",
  "lowapgar5",
  "nicu",
  "sb_subsamp",

  # Maternal complications (6)
  "anaemia_preg_all",
  "gdm_subsamp",
  "gh_subsamp",
  "hdp_subsamp",
  "pe_subsamp",
  "depr_subsamp",
  
  # Other obstetric (4)
  "induction",
  "posttb_all",
  "el_cs",
  "em_cs"
)

# --- [STYLE] US spelling applied to outcome labels (see REVISION LOG) ---
outcome_labels <- c(
    # Bleeding (4)
    Antepartum_bleeding                       = "Antepartum bleeding",
    Postpartum_hemorrhage                     = "Postpartum hemorrhage (any)",
    Postpartum_hemorrhage_due_to_atony        = "PPH due to atony",
    Postpartum_hemorrhage_due_to_retained_placenta = "PPH due to retained placenta",
    
    # Placenta (3)
    finngen_R12_O15_PLAC_PRAEVIA              = "Placenta praevia",
    finngen_R12_O15_PLAC_DISORD               = "Placental disorders",
    finngen_R12_O15_PLAC_PREMAT_SEPAR         = "Premature placental separation",
    
    # Membranes (1)
    rup_memb                                  = "Premature rupture of membranes",
    
    # Birth timing (4)
    pretb_all                                 = "Preterm birth <37 weeks (any)",
    vpretb_all                                = "Very preterm birth <34 weeks",
    posttb_all                                = "Post-term birth",
    ga_all                                    = "Gestational age",
    
    # Fetal Growth (5)
    sga                                       = "Small for gestational age",
    lbw_all                                   = "Low birthweight <2500g",
    hbw_all                                   = "High birthweight >4000g",
    lga                                       = "Large for gestational age",
    zbw_all                                   = "Z-score birthweight",
    
    # Neonatal adaptation / perinatal death (4)
    lowapgar1                                 = "Low Apgar score at 1 min",
    lowapgar5                                 = "Low Apgar score at 5 min",
    nicu                                      = "NICU admission",
    sb_subsamp                                = "Stillbirth",

    # Maternal complications (6)
    anaemia_preg_all                          = "Pregnancy anemia",
    gdm_subsamp                               = "Gestational diabetes",
    gh_subsamp                                = "Gestational hypertension",
    hdp_subsamp                               = "Hypertensive disorders of pregnancy",
    pe_subsamp                                = "Preeclampsia",
    depr_subsamp                              = "Postpartum depression",
    
    # Caesarean section (2)
    el_cs                                     = "Elective caesarean section",
    em_cs                                     = "Emergency caesarean section",
    
    # Other obstetric timing (1)
    induction                                 = "Labour induction"
  )

# Restrict to the 30 outcomes
dat <- dat[dat$outcome %in% vars_keep, , drop = FALSE]

labels_df <- data.frame(
  outcome      = names(outcome_labels),
  outcome_full = unname(outcome_labels),
  stringsAsFactors = FALSE
)

dat <- merge(dat, labels_df, by = "outcome", all.x = TRUE, sort = FALSE)
dat$id.outcome <- dat$outcome

dat <- dat %>%
  dplyr::filter(
    !is.na(beta.outcome),
    !is.na(se.outcome),
    se.outcome > 0,
    !is.na(pval.outcome)
  )

log_info("Outcomes included in tables (n = ", length(unique(dat$outcome)), "):")
print(sort(unique(dat$outcome_full)))

# Continuous outcomes (aligned with sensitivity script)
continuous_outcomes <- c("ga_all", "zbw_all")

###############################################################################
# 2) HELPERS
###############################################################################

label_source <- function(x) {
  dplyr::case_when(
    grepl("^finngen_R12_", x) ~ "FinnGen R12",
    grepl("^Postpartum_hemorrhage", x) ~ "Westergaard (PPH)",
    grepl("^Antepartum_bleeding|^Early_bleeding", x) ~ "Westergaard (PPH)",
    TRUE ~ "MR-PREG"
  )
}

classify_outcome_type <- function(df) {
  df %>%
    dplyr::group_by(outcome) %>%
    dplyr::summarise(
      has_cases = any(!is.na(ncase.outcome)),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      Type = if_else(has_cases, "Binary", "Continuous")
    ) %>%
    dplyr::select(outcome, Type)
}

load_mr_results <- function(stem) {
  path <- file.path(results_dir, paste0(stem, ".csv"))
  if (!file.exists(path)) {
    stop("MR results file not found: ", path)
  }
  df <- readr::read_csv(path, show_col_types = FALSE)
  df <- readr::type_convert(df)
  df
}

fmt_p <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ formatC(p, format = "e", digits = 2),
    p < 0.01  ~ sprintf("%.3f", p),
    TRUE      ~ sprintf("%.2f", p)
  )
}

label_group <- function(outcome, outcome_label) {
  placental_disorders <- c(
    "finngen_R12_O15_PLAC_PRAEVIA",
    "finngen_R12_O15_PLAC_DISORD",
    "finngen_R12_O15_PLAC_PREMAT_SEPAR"
  )
  hdp_outcomes <- c("hdp_subsamp", "gh_subsamp", "pe_subsamp")
  pregnancy_timing <- c("ga_all", "pretb_all", "vpretb_all", "posttb_all")
  fetal_growth_bw <- c("hbw_all", "lbw_all", "sga", "lga", "zbw_all")
  labour_delivery <- c("induction", "rup_memb", "cs", "el_cs", "em_cs")
  bleeding_haem <- c(
    "Antepartum_bleeding", "Early_bleeding_with_any_outcome",
    "Early_bleeding_ending_in_live_birth", "Postpartum_hemorrhage",
    "Postpartum_hemorrhage_due_to_atony", "Postpartum_hemorrhage_due_to_retained_placenta",
    "antepartum_bleeding", "ap_bleeding", "pph_all", "pph_atony",
    "pph_retained_placenta", "pph_other"
  )
  maternal_metab_haem <- c("gdm_subsamp", "anaemia_preg_all")
  maternal_mental     <- c("depr_subsamp")
  neonatal_condition  <- c("lowapgar1", "lowapgar5", "nicu", "sb_subsamp")

  dplyr::case_when(
    outcome %in% placental_disorders ~ "Placental disorders",
    outcome %in% hdp_outcomes ~ "Hypertensive disorders of pregnancy",
    outcome %in% pregnancy_timing |
      grepl("^ga", outcome, ignore.case = TRUE) ~ "Pregnancy timing",
    outcome %in% fetal_growth_bw |
      outcome_label %in% c("Z-score birthweight") ~ "Fetal growth and birthweight",
    outcome %in% labour_delivery ~ "Labour and delivery complications",
    outcome %in% bleeding_haem |
      outcome_label %in% c("Antepartum bleeding", "PPH due to atony",
                           "PPH due to retained placenta",
                           "Postpartum hemorrhage") ~ "Bleeding and haemorrhage",
    outcome %in% maternal_metab_haem ~ "Maternal metabolic/haematologic complications",
    outcome %in% maternal_mental |
      outcome_label %in% c("Postpartum depression", "Perinatal depression") ~
      "Maternal mental health",
    outcome %in% neonatal_condition ~ "Neonatal condition at birth",
    TRUE ~ "Other"
  )
}

snps_file <- file.path(results_dir, "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv")
if (file.exists(snps_file)) {
  clumped_snps  <- data.table::fread(snps_file)
  snp_whitelist <- unique(clumped_snps$SNP)
} else {
  snp_whitelist <- unique(dat$SNP)
  log_warn("Clumped SNP file not found. Using SNPs from harmonised data (n = ",
           length(snp_whitelist), ").")
}
log_info("Using ", length(snp_whitelist), " SNP instruments")

# Outcome types for main MR tables
type_df <- tibble::tibble(
  outcome = vars_keep,
  Type    = dplyr::if_else(outcome %in% continuous_outcomes, "Continuous", "Binary")
)

###############################################################################
# 3) MAIN TABLE 1 - DESCRIPTION OF GWAS DATASETS
###############################################################################

log_info("Creating Table 1: description of GWAS datasets...")

# --- [FIX] McBride et al. (MR-PREG) year corrected 2026 -> 2025 (see REVISION LOG) ---
# --- [R3.3] Koller et al. (2026) row added - sensitivity exposure GWAS used
#     in 05.1_koller_sensitivity_endoMR-PREG.R. Instrument is the 86 EUR loci
#     of the paper's own Supplementary Table 4 (published beta/SE, not
#     re-derived); EUR sample size per the paper's Supplementary Table 1. ---
# --- [EDIT] Table 1 simplified to 7 columns (Study, Year, Dataset/consortium,
#     Phenotype Group, Ancestry, Sample size, Main adjustments); Sex and
#     Notes/accession columns dropped per author request. ---
table1_gwas <- tibble::tribble(
  ~Study,                      ~Year, ~Dataset_or_consortium,                                           ~Phenotype_Group,                                    ~Ancestry,                                                                                                          ~Sample_size,                                                                                                                                          ~Main_adjustments,
  "Rahmioglu et al.",          2023,  "International Endometriosis Genetics Consortium + UK Biobank",   "Endometriosis (exposure)",                          "Predominantly European (~98%) + Japanese (~2%)",                                                                 "60,674 cases; 701,926 controls",                                                                                                                     "Age, batch, principal components, study-specific covariates",
  "McBride et al. (MR-PREG)",  2025,  "MR-PREG collaboration (ALSPAC, BiB, MoBa, UKB, FinnGen + GWAS)", "Adverse pregnancy and perinatal outcomes",          "Predominantly European genetic ancestry; separate South Asian ancestry analyses available for some Born in Bradford outcomes", "Up to 678,001 women",                                                                                                                                "Age, study/centre, principal components",
  "FinnGen R12",               2024,  "FinnGen (release 12)",                                           "Placental phenotypes",                              "Predominantly Finnish/European genetic ancestry",                                                                 "Up to 223,001 women",                                                                                                                                "Age, batch, principal components",
  "Westergaard et al.",        2024,  "Nordic registry-based GWAS (6 cohorts)",                          "Bleeding in pregnancy and postpartum haemorrhage",  "Northern European (registry-based)",                                                                              "Up to 331,792 women",                                                                                                                                "Age, parity, calendar year, cohort",
  "Koller et al.",             2026,  "Multi-ancestry endometriosis GWAS",                              "Endometriosis, sensitivity exposure",               "European (EUR-specific effect estimates, Supplementary Table 4; full GWAS spans 6 ancestries)",                  "99,407 cases; 1,093,534 controls (EUR); 1,388,600 women in the full multi-ancestry GWAS", "Study-specific covariates"
)

table1_file <- file.path(tables_dir, "Table1_GWAS_sources.csv")
write.csv(table1_gwas, table1_file, row.names = FALSE)
log_info("Table 1 (GWAS sources) saved: ", table1_file)

table1_gwas_gt <- table1_gwas %>%
  gt::gt() %>%
  gt::tab_header(
    title = gt::md("**Table 1. Exposure and outcome GWAS datasets included in the MR analyses**")
  ) %>%
  gt::cols_label(
    Study                 = "Study",
    Year                  = "Year",
    Dataset_or_consortium = "Dataset or consortium",
    Phenotype_Group       = "Phenotype Group",
    Ancestry              = "Ancestry",
    Sample_size           = "Sample size",
    Main_adjustments      = "Main adjustments"
  ) %>%
  gt::tab_options(
    table.font.size = 11,
    heading.align   = "left"
  )

table1_gwas_gt

###############################################################################
# 4) MAIN TABLE 2 - SUMMARY OF 30 OUTCOMES & SAMPLE SIZES
###############################################################################

log_info("Creating Table 2: pregnancy outcome sample sizes and SNP counts...")

preg_summary <- dat %>%
  dplyr::filter(SNP %in% snp_whitelist) %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarise(
    `No. SNPs` = dplyr::n_distinct(SNP),
    `N total`  = suppressWarnings(max(samplesize.outcome, na.rm = TRUE)),
    Cases      = suppressWarnings(max(ncase.outcome, na.rm = TRUE)),
    Controls   = suppressWarnings(max(ncontrol.outcome, na.rm = TRUE)),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    `N total` = ifelse(is.infinite(`N total`), NA, `N total`),
    Cases     = ifelse(is.infinite(Cases), NA, Cases),
    Controls  = ifelse(is.infinite(Controls), NA, Controls),
    Source    = label_source(outcome)
  ) %>%
  dplyr::mutate(
    `N total` = ifelse(
      is.na(`N total`) & !is.na(Cases) & !is.na(Controls),
      Cases + Controls,
      `N total`
    )
  )

table1_preg <- preg_summary %>%
  dplyr::mutate(
    Outcome_label = dplyr::recode(outcome, !!!outcome_labels, .default = outcome)
  ) %>%
  dplyr::select(
    Source,
    Outcome = Outcome_label,
    `No. SNPs`,
    `N total`,
    Cases,
    Controls
  ) %>%
  dplyr::arrange(Source, Outcome)

# Manual sizes (kept as in original code; only rows that match an Outcome will be used)
manual_sizes <- tibble::tribble(
  ~Source,             ~Outcome,                             ~`N total`, ~Cases,  ~Controls,
  # FinnGen R12 - sample sizes from FinnGen R12 phenotype manifest
  "FinnGen R12",       "Placenta praevia",                    223001L,     1815L,   221186L,
  "FinnGen R12",       "Placental disorders",                 221519L,      333L,   221186L,
  "FinnGen R12",       "Premature placental separation",      222061L,      875L,   221186L,
  # Westergaard (PPH)
  "Westergaard (PPH)", "Antepartum bleeding",                  331792L,     3236L,   328556L,
  "Westergaard (PPH)", "Early bleeding (live birth)",          331792L,     6356L,   325436L,
  "Westergaard (PPH)", "Early bleeding (any outcome)",         331792L,    28898L,   302894L,
  "Westergaard (PPH)", "PPH due to atony",                     274857L,    13048L,   261809L,
  "Westergaard (PPH)", "PPH due to retained placenta",         272683L,     6256L,   266427L,
  "Westergaard (PPH)", "Postpartum hemorrhage (any)",          331792L,    21521L,   310271L
)

table1_preg <- table1_preg %>%
  dplyr::left_join(
    manual_sizes,
    by = c("Source", "Outcome"),
    suffix = c("", ".manual")
  ) %>%
  dplyr::mutate(
    `N total` = dplyr::coalesce(`N total`, `N total.manual`),
    Cases     = dplyr::coalesce(Cases,     Cases.manual),
    Controls  = dplyr::coalesce(Controls,  Controls.manual)
  ) %>%
  dplyr::select(Source, Outcome, `No. SNPs`, `N total`, Cases, Controls) %>%
  dplyr::arrange(Source, Outcome)

table1_preg_file <- file.path(tables_dir, "Table2_pregnancy_sample_sizes.csv")
write.csv(table1_preg, table1_preg_file, row.names = FALSE)
log_info("Table 2 saved: ", table1_preg_file)

if (requireNamespace("knitr", quietly = TRUE)) {
  cat("\n=== TABLE 2: SAMPLE SIZES ===\n")
  print(knitr::kable(
    table1_preg,
    align   = "llrrrr",
    caption = "Table 2. Sample sizes and SNP counts for pregnancy outcomes."
  ))
}

###############################################################################
# 5) MAIN TABLE 3 - PRIMARY IVW MR ESTIMATES (FDR-CORRECTED)
#    Ordered by domain (same order as Figure 2)
###############################################################################

log_info("Creating Table 3: primary IVW MR estimates with FDR correction...")

# Domain ordering - mirrors outcome_meta in script 07
domain_levels_t3 <- c(
  "Placental disorders",
  "Bleeding & haemorrhage",
  "Pregnancy timing",
  "Labour & delivery",
  "Hypertensive disorders",
  "Fetal growth & birthweight",
  "Maternal metabolic/haematologic",
  "Maternal mental health",
  "Neonatal condition"
)

# --- [REVERTED] "Pregnancy timing" order restored to original (ga_all before rup_memb) ---
domain_meta <- tibble::tribble(
  ~outcome,                                          ~domain,                           ~order_within,
  "finngen_R12_O15_PLAC_PRAEVIA",                    "Placental disorders",             1,
  "finngen_R12_O15_PLAC_DISORD",                     "Placental disorders",             2,
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",               "Placental disorders",             3,
  "Antepartum_bleeding",                             "Bleeding & haemorrhage",          1,
  "Postpartum_hemorrhage",                           "Bleeding & haemorrhage",          2,
  "Postpartum_hemorrhage_due_to_atony",              "Bleeding & haemorrhage",          3,
  "Postpartum_hemorrhage_due_to_retained_placenta",  "Bleeding & haemorrhage",          4,
  "pretb_all",                                       "Pregnancy timing",                1,
  "vpretb_all",                                      "Pregnancy timing",                2,
  "posttb_all",                                      "Pregnancy timing",                3,
  "ga_all",                                          "Pregnancy timing",                4,
  "rup_memb",                                        "Pregnancy timing",                5,
  "induction",                                       "Labour & delivery",               1,
  "el_cs",                                           "Labour & delivery",               2,
  "em_cs",                                           "Labour & delivery",               3,
  "gh_subsamp",                                      "Hypertensive disorders",          1,
  "hdp_subsamp",                                     "Hypertensive disorders",          2,
  "pe_subsamp",                                      "Hypertensive disorders",          3,
  "sga",                                             "Fetal growth & birthweight",      1,
  "lga",                                             "Fetal growth & birthweight",      2,
  "lbw_all",                                         "Fetal growth & birthweight",      3,
  "hbw_all",                                         "Fetal growth & birthweight",      4,
  "zbw_all",                                         "Fetal growth & birthweight",      5,
  "gdm_subsamp",                                     "Maternal metabolic/haematologic", 1,
  "anaemia_preg_all",                                "Maternal metabolic/haematologic", 2,
  "depr_subsamp",                                    "Maternal mental health",           1,
  "lowapgar1",                                       "Neonatal condition",              1,
  "lowapgar5",                                       "Neonatal condition",              2,
  "nicu",                                            "Neonatal condition",              3,
  "sb_subsamp",                                      "Neonatal condition",              4
) %>%
  dplyr::mutate(domain = factor(domain, levels = domain_levels_t3))

ivw_res <- load_mr_results("ivw_results") %>%
  dplyr::filter(
    outcome %in% vars_keep,
    method == "Inverse variance weighted"
  ) %>%
  dplyr::left_join(type_df,     by = "outcome") %>%
  dplyr::left_join(domain_meta, by = "outcome") %>%
  dplyr::mutate(
    Outcome_label = dplyr::recode(outcome, !!!outcome_labels, .default = outcome),
    `Data source` = label_source(outcome),
    q_fdr         = p.adjust(pval, method = "fdr")
  )

table3_main <- ivw_res %>%
  dplyr::mutate(
    Effect_scale = if_else(Type == "Binary", "OR", "beta"),
    Estimate_val = if_else(Type == "Binary", exp(b), b),
    CI_low_val   = if_else(Type == "Binary", exp(b - 1.96 * se), b - 1.96 * se),
    CI_high_val  = if_else(Type == "Binary", exp(b + 1.96 * se), b + 1.96 * se),
    `Estimate`   = if_else(
      Type == "Binary",
      sprintf("%.2f", Estimate_val),
      sprintf("%.3f", Estimate_val)
    ),
    `95% CI`     = if_else(
      Type == "Binary",
      sprintf("(%.2f-%.2f)", CI_low_val, CI_high_val),
      sprintf("(%.3f-%.3f)", CI_low_val, CI_high_val)
    ),
    `P-value` = fmt_p(pval),
    `q (FDR)` = sprintf("%.3f", q_fdr)
  ) %>%
  dplyr::arrange(domain, order_within) %>%
  dplyr::select(
    Domain       = domain,
    Outcome      = Outcome_label,
    Effect_scale,
    `Estimate`,
    `95% CI`,
    `P-value`,
    `q (FDR)`,
    `No. SNPs` = nsnp,
    `Data source`
  )

table3_file <- file.path(tables_dir, "Table3_IVW_FDR_main_results.csv")
write.csv(table3_main, table3_file, row.names = FALSE)
log_info("Table 3 (IVW main FDR) saved: ", table3_file)

if (requireNamespace("knitr", quietly = TRUE)) {
  cat("\n=== TABLE 3: PRIMARY IVW MR ESTIMATES WITH FDR CORRECTION ===\n")
  print(knitr::kable(
    table3_main,
    align = c("l","l","l","l","l","r","r","r","r","r"),
    caption = "Table 3. IVW Mendelian randomization estimates for genetic liability to endometriosis and pregnancy outcomes (FDR-corrected)."
  ))
}

###############################################################################
# 6) SUPPLEMENTARY TABLES 1A-C - DEFINITIONS OF PERINATAL OUTCOMES
###############################################################################

log_info("Creating Supplementary Tables 1A-C (definitions of perinatal outcomes)...")

Supp_1A <- tibble::tribble(
  ~`Binary outcomes`,                ~`Case definition`,                                                ~`Control definition`,                                        ~`Exclusion criteria`,                                                                                         ~`Contributing studies`,
  "Pregnancy loss outcomes",         NA_character_,                                                    NA_character_,                                                NA_character_,                                                                                                 NA_character_,
  "Miscarriage",                     ">=1 pregnancy loss before 20 gestational weeks",                  "No pregnancy loss",                                          "Multiple births",                                                                                             "ALSPAC, MoBa, UKB, FinnGen",
  "Stillbirth",                      ">=1 pregnancy loss at or after 20 gestational weeks",             "No pregnancy loss",                                          "Multiple births",                                                                                             "ALSPAC, BiB, MoBa, UKB",
  "Maternal morbidity outcomes",     NA_character_,                                                    NA_character_,                                                NA_character_,                                                                                                 NA_character_,
  "GDM",                             "Diabetes mellitus diagnosed in pregnancy",                       "No GDM or pre-existing diabetes mellitus",                   "Pre-existing diabetes, multiple births, non-live births",                                                     "BiB, MoBa, UKB, GenDIP",
  "Perinatal depression",            "Maternal depression during pregnancy and up to one year after birth", "No maternal depression",                               "Pre-existing depression, multiple births, non-live births",                                                   "ALSPAC, MoBa, UKB, PGC",
  "Labour outcomes",                 NA_character_,                                                    NA_character_,                                                NA_character_,                                                                                                 NA_character_,
  "Induction of labour",             "Artificial stimulation of uterine contractions",                 "No artificial stimulation of uterine contractions",         "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "Prelabour rupture of membranes",  "Rupture of the amniotic sac before 37 gestational weeks",       "No rupture of the amniotic sac before 37 gestational weeks","Multiple births, non-live births",                                                                           "ALSPAC, MoBa, FinnGen",
  "Caesarean section",               "Delivery by caesarean section",                                 "No delivery by caesarean section",                          "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB, FinnGen",
  "Offspring birth outcomes",        NA_character_,                                                    NA_character_,                                                NA_character_,                                                                                                 NA_character_,
  "LBW",                             "<2,500 g",                                                       ">=2,500 to 4,500 g",                                          "Multiple births, non-live births, PTBs (GA <37 weeks)",                                                      "ALSPAC, MoBa, UKB",
  "HBW",                             ">4,500 g",                                                       ">=2,500 to 4,500 g",                                          "Multiple births, non-live births, PTBs (GA <37 weeks)",                                                      "ALSPAC, MoBa, UKB",
  "Pre-term birth",                  "GA <37 weeks",                                                   "GA >=37 to <42 weeks",                                       "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB, FinnGen, EGG",
  "Post-term birth",                 "GA >=42 weeks",                                                   "GA >=37 to <42 weeks",                                       "Multiple births, non-live births, elective caesarean section, physician-induced labour",                      "ALSPAC, BiB, MoBa, UKB, FinnGen, EGG",
  "SGA",                             "Birth weight <10th percentile for GA",                          "Birth weight >=10th percentile for GA",                      "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "LGA",                             "Birth weight >90th percentile for GA",                          "Birth weight <=90th percentile for GA",                      "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "Low Apgar score at 1 minute",     "<7 points",                                                      ">=7 points",                                                  "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa",
  "Low Apgar score at 5 minutes",    "<7 points",                                                      ">=7 points",                                                  "Multiple births, non-live births",                                                                           "ALSPAC, MoBa",
  "NICU admission",                  "Offspring admitted to the NICU",                                "Offspring not admitted to the NICU",                        "Multiple births, non-live births",                                                                           "ALSPAC, MoBa"
)

s1A_file <- file.path(tables_dir, "Supp_Table_1A_primary_perinatal_definitions.csv")
write.csv(Supp_1A, s1A_file, row.names = FALSE)
log_info("Supplementary Table 1A saved: ", s1A_file)

Supp_1B <- tibble::tribble(
  ~`Binary outcomes`,          ~`Case definition`,                                      ~`Control definition`,                   ~`Exclusion criteria`,                                                                                         ~`Contributing studies`,
  "Pregnancy loss outcomes",   NA_character_,                                          NA_character_,                          NA_character_,                                                                                                 NA_character_,
  "Sporadic miscarriage",      "1-2 pregnancy loss before 20 gestational weeks",       "No pregnancy loss",                    "Age at menarche <9 or >17, underlying conditions leading to miscarriage (see supplement), multiple births",   "ALSPAC, MoBa, UKB",
  "Recurrent miscarriage",     ">=3 pregnancy loss before 20 gestational weeks",        "No pregnancy loss",                    "Age at menarche <9 or >17, underlying conditions leading to miscarriage (see supplement), multiple births",   "ALSPAC, MoBa, UKB, FinnGen",
  "Labour outcomes",           NA_character_,                                          NA_character_,                          NA_character_,                                                                                                 NA_character_,
  "Emergency caesarean section","Delivery by emergency caesarean section",            "No delivery by caesarean section",     "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "Elective caesarean section","Delivery by elective caesarean section",               "No delivery by caesarean section",     "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "Offspring birth outcomes",  NA_character_,                                          NA_character_,                          NA_character_,                                                                                                 NA_character_,
  "Very PTB",                  "GA <34 weeks",                                         "GA >=37 to <42 weeks",                  "Multiple births, non-live births, elective caesarean section, physician-induced labour",                      "ALSPAC, BiB, MoBa, UKB",
  "Spontaneous PTB",           "GA <37 weeks",                                         "GA >=37 to <42 weeks",                  "Multiple births, non-live births, elective caesarean section, physician-induced labour",                      "ALSPAC, BiB, MoBa, UKB"
)

s1B_file <- file.path(tables_dir, "Supp_Table_1B_secondary_binary_perinatal_definitions.csv")
write.csv(Supp_1B, s1B_file, row.names = FALSE)
log_info("Supplementary Table 1B saved: ", s1B_file)

Supp_1C <- tibble::tribble(
  ~`Continuous outcomes`,      ~Units,       ~`Exclusion criteria`,                                                                                                    ~`Contributing studies`,
  "Offspring birth outcomes",  NA_character_,  NA_character_,                                                                                                         NA_character_,
  "Birth weight*",             "SD",         "Multiple births, non-live births, extreme birth weight (>5 SD from sex-specific study mean), PTBs (GA <37 weeks)",      "ALSPAC, BiB, MoBa, UKB",
  "Gestational age*",          "Weeks",      "Multiple births, non-live births",                                                                                      "ALSPAC, BiB, MoBa, UKB, EGG"
)

s1C_file <- file.path(tables_dir, "Supp_Table_1C_continuous_perinatal_definitions.csv")
write.csv(Supp_1C, s1C_file, row.names = FALSE)
log_info("Supplementary Table 1C saved: ", s1C_file)

###############################################################################
# 7) TABLE S1 - HARMONISED DATASET SUMMARY (SNPs PER OUTCOME)
###############################################################################

log_info("Creating Table S1: harmonised dataset summary...")

Table_S1 <- dat %>%
  dplyr::filter(SNP %in% snp_whitelist) %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarise(
    `No. SNPs` = dplyr::n_distinct(SNP),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    Outcome_label = dplyr::recode(outcome, !!!outcome_labels, .default = outcome),
    `Data source` = label_source(outcome)
  ) %>%
  dplyr::left_join(type_df, by = "outcome") %>%
  dplyr::mutate(
    Type = dplyr::if_else(is.na(Type), "Binary", Type)
  ) %>%
  dplyr::left_join(domain_meta, by = "outcome") %>%
  dplyr::mutate(
    domain = factor(domain, levels = domain_levels_t3)
  ) %>%
  dplyr::arrange(domain, order_within) %>%
  dplyr::select(
    Domain        = domain,
    `Data source`,
    Outcome       = outcome,
    Outcome_label,
    Type,
    `No. SNPs`
  )

s1_file <- file.path(tables_dir, "Table_S1_harmonised_summary.csv")
write.csv(Table_S1, s1_file, row.names = FALSE)
log_info("Table S1 (harmonised summary) saved: ", s1_file)

###############################################################################
# 8) SUPPLEMENTARY TABLE 2 - GENETIC INSTRUMENTS FOR ENDOMETRIOSIS
###############################################################################

log_info("Creating Supplementary Table 2: endometriosis instruments...")

if (!exists("coalesce_into")) {
  coalesce_into <- function(df, new_name, candidates) {
    df %>%
      dplyr::mutate(
        !!new_name := dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(candidates)))
      )
  }
}

clumped2_file <- file.path(results_dir, "clumped2_with_stats.rds")
if (!file.exists(clumped2_file)) {
  stop("File 'clumped2_with_stats.rds' not found. Please run 01_select_instruments_endoMR-PREG.R first.")
}
clumped2 <- readRDS(clumped2_file)

clumped2 <- clumped2 %>%
  coalesce_into("beta.exposure",       c("beta.exposure",       "beta.exposure.y",       "beta.exposure.x")) %>%
  coalesce_into("se.exposure",         c("se.exposure",         "se.exposure.y",         "se.exposure.x")) %>%
  coalesce_into("eaf.exposure",        c("eaf.exposure",        "eaf.exposure.y",        "eaf.exposure.x")) %>%
  coalesce_into("pval.exposure",       c("pval.exposure",       "pval.exposure.y",       "pval.exposure.x")) %>%
  coalesce_into("samplesize.exposure", c("samplesize.exposure", "samplesize.exposure.y", "samplesize.exposure.x")) %>%
  coalesce_into("effect_allele.exposure",
                c("effect_allele.exposure", "effect_allele.exposure.y", "effect_allele.exposure.x", "effect_allele")) %>%
  coalesce_into("other_allele.exposure",
                c("other_allele.exposure",  "other_allele.exposure.y",  "other_allele.exposure.x",  "other_allele"))

if (!all(c("R2_i", "F_i") %in% names(clumped2))) {
  clumped2 <- clumped2 %>%
    dplyr::mutate(
      R2_i = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
      F_i  = R2_i * (samplesize.exposure - 2) / (1 - R2_i)
    )
}

supp_table2 <- clumped2 %>%
  dplyr::select(
    SNP,
    Beta            = beta.exposure,
    SE              = se.exposure,
    `P-value`       = pval.exposure,
    `Effect allele` = effect_allele.exposure,
    `Other allele`  = other_allele.exposure,
    EAF             = eaf.exposure,
    r2              = R2_i,
    F               = F_i
  ) %>%
  dplyr::mutate(
    Beta    = round(Beta, 4),
    SE      = round(SE, 4),
    EAF     = round(EAF, 4),
    r2      = signif(r2, 3),
    F       = round(F, 2),
    `P-value` = formatC(`P-value`, format = "e", digits = 2)
  )

supp2_file <- file.path(tables_dir, "Supplementary_Table_2_Endometriosis_Instruments.xlsx")
openxlsx::write.xlsx(
  supp_table2,
  supp2_file,
  rowNames = FALSE
)
log_info("Supplementary Table 2 saved: ", supp2_file)

###############################################################################
# 9) SUPPLEMENTARY TABLE S3 - SENSITIVITY ANALYSES (WIDE FORMAT, eBioMedicine)
#
#   One row per outcome, ordered by the 9 clinical domains (same as Figure 2 /
#   Table 3). Columns: Domain, Scale (OR / beta), IVW, MR-Egger, Weighted Median,
#   Weighted Mode effect estimates, Egger intercept, Cochran's Q p-value, and
#   MR-PRESSO as a diagnostic block.
#   Effect estimates are ORs for binary outcomes and beta for continuous outcomes;
#   all column headers use "effect (95% CI)" to avoid misleading OR labelling
#   for continuous outcomes.
###############################################################################

log_info("Creating Table S3: wide-format sensitivity analyses (eBioMedicine)...")

# ── Helper functions ──────────────────────────────────────────────────────────

# Scalar: returns "X.XX (X.XX-X.XX)" for OR (binary) or "X.XXX (X.XXX-X.XXX)" for beta
fmt_est_ci <- function(b, se, type) {
  if (is.na(b) || is.na(se)) return(NA_character_)
  if (type == "binary") {
    sprintf("%.2f (%.2f\u2013%.2f)", exp(b), exp(b - 1.96*se), exp(b + 1.96*se))
  } else {
    sprintf("%.3f (%.3f\u2013%.3f)", b, b - 1.96*se, b + 1.96*se)
  }
}

# Vectorized p-value formatter (1 decimal scientific notation below 0.001)
fmt_p_s3 <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ NA_character_,
    p < 0.001 ~ formatC(p, format = "e", digits = 1),
    TRUE      ~ sprintf("%.3f", p)
  )
}

# ── Load all-methods results ──────────────────────────────────────────────────
all_res_s3 <- load_mr_results("all_mr_methods") %>%
  dplyr::filter(
    outcome %in% vars_keep,
    method  %in% c("Inverse variance weighted", "MR Egger",
                   "Weighted median", "Weighted mode")
  ) %>%
  dplyr::mutate(
    type       = dplyr::if_else(outcome %in% continuous_outcomes, "continuous", "binary"),
    method_key = dplyr::recode(method,
      "Inverse variance weighted" = "IVW",
      "MR Egger"                  = "Egger",
      "Weighted median"           = "WM",
      "Weighted mode"             = "WMode"
    )
  )

# ── Pivot to wide: one row per outcome ───────────────────────────────────────
s3_wide <- all_res_s3 %>%
  dplyr::select(outcome, type, method_key, b, se, pval) %>%
  tidyr::pivot_wider(
    names_from  = method_key,
    values_from = c(b, se, pval),
    names_glue  = "{method_key}_{.value}"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    # "effect (95% CI)" headers - generic label covers both OR and beta
    `IVW effect (95% CI)`             = fmt_est_ci(IVW_b,   IVW_se,   type),
    `IVW p-value`                     = fmt_p_s3(IVW_pval),
    IVW_p_raw                         = IVW_pval,   # kept for Excel bold logic
    `MR-Egger effect (95% CI)`        = fmt_est_ci(Egger_b,  Egger_se,  type),
    `Weighted median effect (95% CI)` = fmt_est_ci(WM_b,    WM_se,    type),
    `Weighted mode effect (95% CI)`   = fmt_est_ci(WMode_b, WMode_se, type),
    # Scale column: tells readers which metric each row uses
    Scale = dplyr::if_else(type == "binary", "OR", "\u03b2")
  ) %>%
  dplyr::ungroup()

# ── Egger intercept ───────────────────────────────────────────────────────────
# Format: "beta = X.XXX (p = X.XXX)" - explicit beta symbol clarifies it is always
# a regression intercept on the log-OR / beta scale, regardless of outcome type.
pleio_file <- file.path(results_dir, "mr_pleiotropy.csv")
if (file.exists(pleio_file)) {
  pleio_raw <- readr::read_csv(pleio_file, show_col_types = FALSE) %>%
    dplyr::filter(outcome %in% vars_keep) %>%
    dplyr::mutate(
      `Egger intercept (beta), p-value` = sprintf(
        "\u03b2\u00a0=\u00a0%.4f (p\u00a0=\u00a0%s)",
        as.numeric(egger_intercept),
        fmt_p_s3(as.numeric(pval))
      )
    ) %>%
    dplyr::select(outcome, `Egger intercept (beta), p-value`)
} else {
  pleio_raw <- tibble::tibble(
    outcome                           = character(),
    `Egger intercept (beta), p-value` = character()
  )
  log_warn("mr_pleiotropy.csv not found - Egger intercept column will be NA.")
}

# ── MR-RAPS (winner's-curse-robust estimator) ────────────────────────────────
# --- [R4.6] Added as an extra column rather than a standalone table, since it
#     was run across all 30 outcomes exactly like the other sensitivity
#     estimators already in this table. Source: results/mr_raps.csv (script
#     05, section 3b). ---
raps_file <- file.path(results_dir, "mr_raps.csv")
if (file.exists(raps_file)) {
  raps_raw <- readr::read_csv(raps_file, show_col_types = FALSE) %>%
    dplyr::filter(outcome %in% vars_keep) %>%
    dplyr::mutate(
      type_raps = dplyr::if_else(outcome %in% continuous_outcomes, "continuous", "binary"),
      `MR-RAPS effect (95% CI)` = purrr::pmap_chr(
        list(b, se, type_raps), ~ fmt_est_ci(..1, ..2, ..3)
      ),
      `MR-RAPS p-value` = fmt_p_s3(as.numeric(pval))
    ) %>%
    dplyr::select(outcome, `MR-RAPS effect (95% CI)`, `MR-RAPS p-value`)
} else {
  raps_raw <- tibble::tibble(
    outcome                    = character(),
    `MR-RAPS effect (95% CI)`  = character(),
    `MR-RAPS p-value`          = character()
  )
  log_warn("mr_raps.csv not found - MR-RAPS columns will be NA.")
}

# ── Cochran's Q (IVW heterogeneity) ──────────────────────────────────────────
het_file <- file.path(results_dir, "mr_heterogeneity.csv")
if (file.exists(het_file)) {
  het_raw <- readr::read_csv(het_file, show_col_types = FALSE) %>%
    dplyr::filter(outcome %in% vars_keep,
                  method  == "Inverse variance weighted") %>%
    dplyr::mutate(`Cochran's Q p-value` = fmt_p_s3(as.numeric(Q_pval))) %>%
    dplyr::select(outcome, `Cochran's Q p-value`)
} else {
  # Compute on the fly from harmonised data
  dat_s3_het <- data.table::fread(harm_file)
  dat_s3_het <- dat_s3_het[dat_s3_het$outcome %in% vars_keep, ]
  het_raw <- TwoSampleMR::mr_heterogeneity(dat_s3_het) %>%
    dplyr::filter(method == "Inverse variance weighted") %>%
    dplyr::mutate(`Cochran's Q p-value` = fmt_p_s3(Q_pval)) %>%
    dplyr::select(outcome, `Cochran's Q p-value`)
  log_warn("mr_heterogeneity.csv not found - computed on the fly.")
}

# ── MR-PRESSO (diagnostic columns) ───────────────────────────────────────────
presso_file <- file.path(results_dir, "mr_presso_results.csv")
if (file.exists(presso_file)) {
  presso_raw_s3 <- readr::read_csv(presso_file, show_col_types = FALSE) %>%
    dplyr::filter(outcome_id %in% vars_keep) %>%
    dplyr::rename(outcome = outcome_id) %>%
    dplyr::mutate(
      type_presso = dplyr::if_else(outcome %in% continuous_outcomes, "continuous", "binary"),
      `MR-PRESSO outliers (n)` = as.integer(n_outliers),
      `MR-PRESSO corrected effect (95% CI)` = dplyr::if_else(
        !is.na(n_outliers) & n_outliers >= 1,
        purrr::pmap_chr(
          list(corrected_estimate, se_corrected, type_presso),
          ~ fmt_est_ci(..1, ..2, ..3)
        ),
        "\u2014"   # em-dash when no outliers detected
      ),
      `MR-PRESSO distortion p-value` = dplyr::if_else(
        !is.na(n_outliers) & n_outliers >= 1,
        fmt_p_s3(pval_corrected),
        "\u2014"
      )
    ) %>%
    dplyr::select(
      outcome,
      `MR-PRESSO outliers (n)`,
      `MR-PRESSO corrected effect (95% CI)`,
      `MR-PRESSO distortion p-value`
    )
} else {
  presso_raw_s3 <- tibble::tibble(
    outcome                                = character(),
    `MR-PRESSO outliers (n)`               = integer(),
    `MR-PRESSO corrected effect (95% CI)`  = character(),
    `MR-PRESSO distortion p-value`         = character()
  )
  log_warn("mr_presso_results.csv not found - PRESSO columns will be NA.")
}

# ── Join all components + domain ordering ─────────────────────────────────────
Table_S3_wide <- s3_wide %>%
  dplyr::left_join(pleio_raw,     by = "outcome") %>%
  dplyr::left_join(raps_raw,      by = "outcome") %>%
  dplyr::left_join(het_raw,       by = "outcome") %>%
  dplyr::left_join(presso_raw_s3, by = "outcome") %>%
  dplyr::left_join(
    domain_meta %>% dplyr::mutate(outcome_id = outcome),
    by = c("outcome" = "outcome_id")
  ) %>%
  dplyr::mutate(
    Outcome = dplyr::recode(outcome, !!!outcome_labels, .default = outcome),
    Domain  = factor(domain, levels = domain_levels_t3),
    # Replace NA PRESSO cells with em-dash
    `MR-PRESSO outliers (n)` = dplyr::coalesce(`MR-PRESSO outliers (n)`, 0L),
    `MR-PRESSO corrected effect (95% CI)` = dplyr::if_else(
      is.na(`MR-PRESSO corrected effect (95% CI)`), "\u2014",
      `MR-PRESSO corrected effect (95% CI)`
    ),
    `MR-PRESSO distortion p-value` = dplyr::if_else(
      is.na(`MR-PRESSO distortion p-value`), "\u2014",
      `MR-PRESSO distortion p-value`
    )
  ) %>%
  dplyr::arrange(Domain, order_within) %>%
  dplyr::select(
    Domain,
    Outcome,
    Scale,                                  # OR or beta - clarifies metric per row
    `IVW effect (95% CI)`,
    `IVW p-value`,
    `MR-Egger effect (95% CI)`,
    `Egger intercept (beta), p-value`,
    `Weighted median effect (95% CI)`,
    `Weighted mode effect (95% CI)`,
    `MR-RAPS effect (95% CI)`,
    `MR-RAPS p-value`,
    `Cochran's Q p-value`,
    `MR-PRESSO outliers (n)`,
    `MR-PRESSO corrected effect (95% CI)`,
    `MR-PRESSO distortion p-value`,
    IVW_p_raw   # used for Excel bold; removed before export
  )

# ── CSV export ────────────────────────────────────────────────────────────────
Table_S3_csv <- Table_S3_wide %>% dplyr::select(-IVW_p_raw)
readr::write_csv(
  Table_S3_csv,
  file.path(tables_dir, "Supplementary_Table_S3_sensitivity_analyses.csv")
)
log_info("Table S3 (wide format, CSV) saved.")

# ── Excel export with bold for IVW p < 0.05 ──────────────────────────────────
if (requireNamespace("openxlsx", quietly = TRUE)) {

  wb_s3 <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb_s3, "S3 Sensitivity analyses")

  # Caption at row 1
  caption_text <- paste0(
    "Supplementary Table S3. Sensitivity analyses for the two-sample Mendelian ",
    "randomization of genetic liability to endometriosis on pregnancy outcomes."
  )
  openxlsx::writeData(wb_s3, "S3 Sensitivity analyses",
                      x = caption_text, startRow = 1, startCol = 1)
  openxlsx::addStyle(wb_s3, "S3 Sensitivity analyses",
                     style = openxlsx::createStyle(
                       fontName = "Arial", fontSize = 10, textDecoration = "bold",
                       wrapText = TRUE
                     ),
                     rows = 1, cols = 1)

  # Table starts at row 3
  tbl_s3_export <- Table_S3_wide %>% dplyr::select(-IVW_p_raw)
  openxlsx::writeDataTable(
    wb_s3, "S3 Sensitivity analyses",
    x          = tbl_s3_export,
    startRow   = 3, startCol = 1,
    tableStyle = "TableStyleLight9",
    withFilter = TRUE
  )

  # Bold IVW effect (95% CI) cell (col 4) where IVW p < 0.05
  # Offset: +3 rows (caption row 1 + blank row 2 + header row 3)
  bold_rows <- which(!is.na(Table_S3_wide$IVW_p_raw) &
                       Table_S3_wide$IVW_p_raw < 0.05) + 3
  if (length(bold_rows) > 0) {
    openxlsx::addStyle(wb_s3, "S3 Sensitivity analyses",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 9,
                         textDecoration = "bold"
                       ),
                       rows = bold_rows, cols = 4,   # col 4 = "IVW effect (95% CI)"
                       gridExpand = TRUE)
  }

  # Column widths: Domain, Outcome, Scale, IVW eff, IVW p, Egger eff,
  #                Egger intercept, WM eff, WMode eff, RAPS eff, RAPS p,
  #                Q p, PRESSO n, PRESSO corr, PRESSO distort
  openxlsx::setColWidths(
    wb_s3, "S3 Sensitivity analyses",
    cols   = seq_len(ncol(tbl_s3_export)),
    widths = c(26, 32, 7, 22, 12, 22, 30, 22, 22, 22, 12, 14, 16, 26, 22)
  )

  # Header style
  openxlsx::addStyle(wb_s3, "S3 Sensitivity analyses",
                     style = openxlsx::createStyle(
                       fontName = "Arial", fontSize = 10, textDecoration = "bold",
                       fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                       wrapText = TRUE, valign = "top"
                     ),
                     rows = 3, cols = seq_len(ncol(tbl_s3_export)),
                     gridExpand = TRUE)

  # Alternate domain shading: light grey for odd-numbered domains
  domain_vals <- as.integer(Table_S3_wide$Domain)
  shade_rows  <- which(domain_vals %% 2 == 1) + 3   # odd domains get shading
  if (length(shade_rows) > 0) {
    openxlsx::addStyle(wb_s3, "S3 Sensitivity analyses",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 9, fgFill = "#F5F7FB"
                       ),
                       rows = shade_rows,
                       cols = seq_len(ncol(tbl_s3_export)),
                       gridExpand = TRUE, stack = TRUE)
  }

  # Footnote at bottom
  footnote_row <- nrow(tbl_s3_export) + 5
  footnote_text <- paste0(
    "Abbreviations: IVW, inverse variance weighted; MR, Mendelian randomization; ",
    "OR, odds ratio; \u03b2, beta coefficient; SNP, single-nucleotide polymorphism. ",
    "Scale column indicates the effect metric used: OR for binary outcomes and \u03b2 for continuous outcomes ",
    "(gestational age and birthweight z-score). ",
    "All effect estimates and 95% confidence intervals in the \u2018effect (95% CI)\u2019 columns are expressed ",
    "on the OR scale for binary outcomes and as \u03b2 coefficients for continuous outcomes; ",
    "no log-transformation has been applied to continuous outcomes. ",
    "Cochran\u2019s Q p-value tests heterogeneity of SNP-specific causal estimates under the IVW model; ",
    "a low p-value indicates evidence of heterogeneity, which may reflect pleiotropy. ",
    "The Egger intercept (\u03b2) tests for directional (unbalanced) horizontal pleiotropy; ",
    "departure from zero suggests that InSIDE assumption violations may bias IVW estimates. ",
    "MR-RAPS (Robust Adjusted Profile Score) is a winner's-curse-robust estimator run across ",
    "all 30 outcomes for consistency with the other sensitivity methods above. ",
    "MR-PRESSO outlier-corrected estimates and distortion p-values are shown only for outcomes ",
    "where at least one outlier SNP was detected; \u2014 indicates no outliers. ",
    "The distortion test evaluates whether removal of outlier SNPs significantly changes the causal estimate. ",
    "Bold IVW effect estimates indicate p\u00a0<\u00a00.05 (nominal significance, uncorrected for multiple testing)."
  )
  openxlsx::writeData(wb_s3, "S3 Sensitivity analyses",
                      x = footnote_text,
                      startRow = footnote_row, startCol = 1)
  openxlsx::addStyle(wb_s3, "S3 Sensitivity analyses",
                     style = openxlsx::createStyle(
                       fontName = "Arial", fontSize = 8, fontColour = "#595959",
                       wrapText = TRUE, valign = "top"
                     ),
                     rows = footnote_row, cols = 1)
  openxlsx::setRowHeights(wb_s3, "S3 Sensitivity analyses",
                          rows = footnote_row, heights = 90)
  openxlsx::mergeCells(wb_s3, "S3 Sensitivity analyses",
                       cols = 1:ncol(tbl_s3_export), rows = footnote_row)

  openxlsx::saveWorkbook(
    wb_s3,
    file.path(tables_dir, "Supplementary_Table_S3_sensitivity_analyses.xlsx"),
    overwrite = TRUE
  )
  log_info("Table S3 (wide format, Excel) saved.")
}

message("Supplementary Table S3 done.")

# Keep s3_file reference for downstream compatibility
s3_file <- file.path(tables_dir, "Supplementary_Table_S3_sensitivity_analyses.csv")

###############################################################################
# 10) TABLE S4 - TRIO-BASED MR ESTIMATES (MATERNAL/FETAL/PATERNAL)
#
# --- [R3.2, REVERTED] The [R3.2] change that replaced this table with
#     heterogeneity statistics (Cochran's Q) was itself reverted: the
#     manuscript's own "TABLES AND FIGURES" list and the Results text
#     ("Intergenerational (trio-based) analyses", citing "Supplementary
#     Table S4") both expect S4 to be the trio table, and the
#     heterogeneity/Cochran's Q content is already fully covered by the
#     "Cochran's Q p-value" column already present in Table S3 - so nothing
#     is lost by dropping the standalone heterogeneity table.
#     BLOCKER FIX: the version of this table previously pasted into the
#     manuscript did not match the Results text at all for Fetal/Paternal
#     (its values traced to no file anywhere in results/ - apparently a
#     lost/stale independent-single-genotype model). This version is built
#     directly from results/trios_adj_mr_results_by_outcome_long.csv
#     (script 04.2, mutually-adjusted DONUTS model), restricted to the 3
#     "conditioned" estimates the Methods promise ("Maternal-Fetal Genotype
#     Correction": maternal adj. for fetal = primary, fetal adj. for
#     maternal = sensitivity, paternal = negative control) - verified to
#     reproduce the exact numbers quoted in Results (e.g. pretb_all maternal
#     adj. OR=0.71/q=0.007, paternal OR=0.78/p=0.118). The unadjusted
#     maternal estimate (quoted once in Results, as an illustrative
#     before/after-adjustment contrast for preterm birth only) is not a
#     standing row here to avoid re-introducing a 4th ambiguous category. ---
###############################################################################

log_info("Creating Table S4: trio-based MR estimates (maternal/fetal/paternal)...")

trios_long_file <- file.path(results_dir, "trios_adj_mr_results_by_outcome_long.csv")
if (!file.exists(trios_long_file)) {
  stop("File 'trios_adj_mr_results_by_outcome_long.csv' not found. Run 04.2_fetal_effect_endoMR-PREG.R first.")
}

trios_long <- readr::read_csv(trios_long_file, show_col_types = FALSE)

Table_S4 <- trios_long %>%
  dplyr::filter(origin %in% c("Maternal (adj. fetal)", "Fetal (adj.)", "Paternal (adj.)")) %>%
  dplyr::mutate(
    Outcome = dplyr::recode(outcome_id, !!!outcome_labels, .default = outcome_full),
    origin  = dplyr::recode(origin,
      "Maternal (adj. fetal)" = "Maternal (adjusted for fetal genotype)",
      "Fetal (adj.)"          = "Fetal (adjusted for maternal genotype)",
      "Paternal (adj.)"       = "Paternal (adjusted, negative control)"
    ),
    b   = round(b, 3),
    se  = round(se, 3),
    OR  = dplyr::if_else(scale == "binary", round(OR, 3), NA_real_),
    LCL = dplyr::if_else(scale == "binary", round(LCL, 3), NA_real_),
    UCL = dplyr::if_else(scale == "binary", round(UCL, 3), NA_real_),
    pval = signif(pval, 3),
    qval = round(qval, 3)
  ) %>%
  dplyr::select(Outcome, origin, nsnp, b, se, OR, LCL, UCL, pval, qval) %>%
  dplyr::arrange(Outcome, origin)

s4_file <- file.path(tables_dir, "Supplementary_Table_S4_trio_based_MR.csv")
write.csv(Table_S4, s4_file, row.names = FALSE)
log_info("Table S4 (trio-based MR) saved: ", s4_file)

log_info("=== Table generation for endoMR-PREG completed successfully ===")

###############################################################################
# 12) EXPORT ALL TABLES TO HTML
###############################################################################

log_info("Exporting all tables to HTML...")

suppressPackageStartupMessages(library(gt))

html_dir <- file.path(tables_dir, "html")
dir.create(html_dir, showWarnings = FALSE, recursive = TRUE)

# --- [EDIT] Table names/numbering reworked to match final manuscript numbering (see REVISION LOG) ---
tables_to_export <- list(
  Table1_GWAS_sources                      = table1_gwas,
  Table2_pregnancy_sample_sizes            = table1_preg,
  Table3_IVW_FDR_main_results              = table3_main,
  Supp_Table_1A_primary_perinatal_def      = Supp_1A,
  Supp_Table_1B_secondary_binary_perinatal = Supp_1B,
  Supp_Table_1C_continuous_perinatal       = Supp_1C,
  Table_S1_harmonised_summary              = Table_S1,
  Supplementary_Table_2_Endometriosis_Instruments = supp_table2,
  Table_S3_all_MR_methods                  = Table_S3_csv,
  Table_S4_trio_based_MR                   = Table_S4
)

table_titles <- list(
  Table1_GWAS_sources =
    "Table 1. Exposure and outcome GWAS datasets included in the MR analyses",
  
  Table2_pregnancy_sample_sizes =
    "Table 2. Pregnancy outcomes and sample sizes in MR-PREG and related cohorts",
  
  Table3_IVW_FDR_main_results =
    "Table 3. Main inverse-variance weighted Mendelian randomisation results with FDR correction",
  
  Supp_Table_1A_primary_perinatal_def =
    "Supplementary Table 1A. Definitions and sources of primary perinatal outcomes",
  
  Supp_Table_1B_secondary_binary_perinatal =
    "Supplementary Table 1B. Definitions and sources of secondary binary perinatal outcomes",
  
  Supp_Table_1C_continuous_perinatal =
    "Supplementary Table 1C. Definitions and sources of continuous perinatal outcomes",
  
  Table_S1_harmonised_summary =
    "Table S1. Summary of harmonised exposure-outcome datasets used in MR analyses",
  
  Supplementary_Table_2_Endometriosis_Instruments =
    "Supplementary Table 2. Genetic instruments for endometriosis liability (Rahmioglu et al.)",
  
  Table_S3_all_MR_methods =
    "Supplementary Table S3. Sensitivity analyses across all MR methods (IVW, MR-Egger, Weighted Median, Weighted Mode, MR-PRESSO) for genetic liability to endometriosis on pregnancy outcomes",
  
  Table_S4_trio_based_MR =
    "Supplementary Table S4. Trio-based Mendelian randomization estimates for the effect of genetic liability to endometriosis on pregnancy and perinatal outcomes by maternal, fetal, and paternal origin"
)

for (nm in names(tables_to_export)) {
  df <- tables_to_export[[nm]]

  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    next
  }

  title_text <- table_titles[[nm]]
  if (is.null(title_text)) {
    title_text <- nm
  }

  gt_tbl <- df %>%
    gt::gt() %>%
    gt::tab_header(
      title = gt::md(paste0("**", title_text, "**"))
    )

  out_file <- file.path(html_dir, paste0(nm, ".html"))

  gt::gtsave(gt_tbl, out_file)
  log_info("Saved HTML table: ", out_file)
}

###############################################################################
# SUPPLEMENTARY TABLE S5 - LEAVE-ONE-COHORT SENSITIVITY ANALYSIS
#
#   IVW estimates across all outcomes with each contributing study
#   sequentially excluded (leave-one-cohort-out design).
#   Source: results/leave_one_cohort_results.csv (produced by script 05).
#   Outcomes and domains ordered as in Table 3 / Figure 2.
###############################################################################

loc_path <- file.path(results_dir, "leave_one_cohort_results.csv")

if (!file.exists(loc_path)) {
  warning("leave_one_cohort_results.csv not found - skipping Table S5.")
} else {

  loc_raw <- readr::read_csv(loc_path, show_col_types = FALSE)

  # Cohort labels: make them consistent and reader-friendly
  cohort_label_map <- c(
    "Main analysis"   = "Main analysis (all cohorts)",
    "Without ALSPAC"  = "Excluding ALSPAC",
    "Without BIB-SA"  = "Excluding BiB (South Asian)",
    "Without BIB-WE"  = "Excluding BiB (White European)",
    "Without FinnGen" = "Excluding FinnGen",
    "Without MOBA"    = "Excluding MoBa",
    "Without Public"  = "Excluding public GWAS",
    "Without UKB"     = "Excluding UK Biobank"
  )

  loc_fmt <- loc_raw %>%
    dplyr::filter(method == "Inverse variance weighted") %>%
    dplyr::left_join(
      domain_meta %>%
        dplyr::mutate(outcome_id = outcome) %>%
        dplyr::select(outcome_id, domain, order_within),
      by = c("id.outcome" = "outcome_id")
    ) %>%
    dplyr::mutate(
      Outcome = dplyr::recode(id.outcome, !!!outcome_labels, .default = id.outcome),
      Domain  = factor(domain, levels = domain_levels_t3),
      # Effect estimates: OR for binary, beta for continuous
      est = dplyr::if_else(outcome_type == "binary", exp(b),             b),
      lci = dplyr::if_else(outcome_type == "binary", exp(b - 1.96 * se), b - 1.96 * se),
      uci = dplyr::if_else(outcome_type == "binary", exp(b + 1.96 * se), b + 1.96 * se),
      Scale = dplyr::if_else(outcome_type == "binary", "OR", "\u03b2"),
      `Estimate (95% CI)` = dplyr::if_else(
        outcome_type == "binary",
        sprintf("%.2f (%.2f\u2013%.2f)", est, lci, uci),
        sprintf("%.3f (%.3f\u2013%.3f)", est, lci, uci)
      ),
      `P value` = dplyr::case_when(
        pval < 0.001 ~ formatC(pval, format = "e", digits = 1),
        TRUE         ~ sprintf("%.3f", pval)
      ),
      Cohort = dplyr::recode(analysis_type, !!!cohort_label_map, .default = analysis_type),
      # Order: Main analysis first, then cohorts alphabetically
      cohort_order = dplyr::if_else(analysis_type == "Main analysis", 0L, 1L)
    ) %>%
    dplyr::filter(!is.na(Domain)) %>%
    dplyr::arrange(Domain, order_within, cohort_order, Cohort) %>%
    dplyr::select(
      Domain,
      Outcome,
      `Cohort excluded` = Cohort,
      Scale,
      `Estimate (95% CI)`,
      `P value`,
      `SNPs (n)`         = nsnp
    ) %>%
    dplyr::mutate(Domain = as.character(Domain))

  # ── CSV ──────────────────────────────────────────────────────────────────────
  readr::write_csv(
    loc_fmt,
    file.path(tables_dir, "Supplementary_Table_S5_leave_one_cohort.csv")
  )
  log_info("Supplementary Table S5 (CSV) saved.")

  # ── Excel ────────────────────────────────────────────────────────────────────
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    wb_s5 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s5, "S5 Leave-one-cohort")

    header_style_s5 <- openxlsx::createStyle(
      fontName       = "Arial", fontSize = 10, textDecoration = "bold",
      fgFill         = "#D9E1F2", border = "Bottom", borderColour = "#4472C4"
    )
    main_style_s5 <- openxlsx::createStyle(
      fontName = "Arial", fontSize = 9,
      textDecoration = "bold", fgFill = "#EEF2F8"
    )
    sens_style_s5 <- openxlsx::createStyle(
      fontName = "Arial", fontSize = 9
    )

    openxlsx::writeDataTable(
      wb_s5, "S5 Leave-one-cohort",
      x           = loc_fmt,
      startRow    = 1, startCol = 1,
      tableStyle  = "TableStyleLight9",
      withFilter  = TRUE
    )

    # Bold + shaded rows for "Main analysis (all cohorts)"
    main_rows_s5 <- which(loc_fmt$`Cohort excluded` == "Main analysis (all cohorts)") + 1
    if (length(main_rows_s5) > 0) {
      openxlsx::addStyle(wb_s5, "S5 Leave-one-cohort",
                         style      = main_style_s5,
                         rows       = main_rows_s5,
                         cols       = seq_len(ncol(loc_fmt)),
                         gridExpand = TRUE)
    }

    openxlsx::setColWidths(wb_s5, "S5 Leave-one-cohort",
                           cols   = seq_len(ncol(loc_fmt)),
                           widths = c(28, 34, 28, 8, 22, 10, 8))

    openxlsx::saveWorkbook(
      wb_s5,
      file.path(tables_dir, "Supplementary_Table_S5_leave_one_cohort.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S5 (Excel) saved.")
  }

  message("Supplementary Table S5 (leave-one-cohort) done.")
}

###############################################################################
# SUPPLEMENTARY TABLE S6 - LOCUS-LEVEL COLOCALISATION (ENDOMETRIOSIS x
#                           PLACENTA PRAEVIA)
#
#   [R2.1] Reviewer #2, comment 1: coloc.abf results at each of the 41
#   independent loci contributing to the endometriosis instrument, testing
#   for a shared causal variant with placenta praevia (FinnGen R12).
#   Source: results/coloc_endo_praevia_by_locus.csv (script 05.3).
#   Ordered by descending PP.H4 (posterior probability of colocalisation).
###############################################################################

coloc_path <- file.path(results_dir, "coloc_endo_praevia_by_locus.csv")

if (!file.exists(coloc_path)) {
  warning("coloc_endo_praevia_by_locus.csv not found - skipping Table S6.")
} else {

  coloc_raw <- readr::read_csv(coloc_path, show_col_types = FALSE)

  # Same window as script 05.3 (coloc.abf regional extraction)
  WINDOW_BP <- 500000L

  hyp_labels <- c(
    PP.H0 = "H0 (no association)",
    PP.H1 = "H1 (endometriosis only)",
    PP.H2 = "H2 (placenta praevia only)",
    PP.H3 = "H3 (distinct causal variants)",
    PP.H4 = "H4 (shared causal variant)"
  )

  coloc_fmt <- coloc_raw %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      `Predominant hypothesis` = hyp_labels[[
        c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")[
          which.max(c(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4))
        ]
      ]]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(PP.H4)) %>%
    dplyr::mutate(
      `Region start` = pos - WINDOW_BP,
      `Region end`   = pos + WINDOW_BP,
      PP.H0 = round(PP.H0, 3),
      PP.H1 = round(PP.H1, 3),
      PP.H2 = round(PP.H2, 3),
      PP.H3 = round(PP.H3, 3),
      PP.H4 = round(PP.H4, 3)
    ) %>%
    dplyr::select(
      `Index SNP`                    = SNP,
      Chromosome                     = chr,
      `Index position`               = pos,
      `Region start`,
      `Region end`,
      `Number of overlapping SNPs`   = n_snps,
      PP.H0, PP.H1, PP.H2, PP.H3, PP.H4,
      `Predominant hypothesis`
    )

  # ── CSV ──────────────────────────────────────────────────────────────────────
  readr::write_csv(
    coloc_fmt,
    file.path(tables_dir, "Supplementary_Table_S6_colocalisation_endo_praevia.csv")
  )
  log_info("Supplementary Table S6 (CSV) saved.")

  # ── Excel ────────────────────────────────────────────────────────────────────
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    wb_s6 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s6, "S6 Colocalisation")

    caption_text_s6 <- paste0(
      "Supplementary Table S6. Locus-level colocalisation analyses of endometriosis and placenta ",
      "praevia. Colocalisation (coloc.abf, Giambartolomei et al. 2014) between genetic liability ",
      "to endometriosis (Rahmioglu et al. 2023) and placenta praevia (FinnGen R12) at each of the ",
      "41 independent loci contributing to the endometriosis instrument. Loci are ordered by ",
      "descending PP.H4 (posterior probability of a shared causal variant)."
    )
    openxlsx::writeData(wb_s6, "S6 Colocalisation",
                        x = caption_text_s6, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb_s6, "S6 Colocalisation",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         wrapText = TRUE
                       ),
                       rows = 1, cols = 1)
    openxlsx::mergeCells(wb_s6, "S6 Colocalisation", cols = 1:ncol(coloc_fmt), rows = 1)

    openxlsx::writeDataTable(
      wb_s6, "S6 Colocalisation",
      x          = coloc_fmt,
      startRow   = 3, startCol = 1,
      tableStyle = "TableStyleLight9",
      withFilter = TRUE
    )

    openxlsx::addStyle(wb_s6, "S6 Colocalisation",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = 3, cols = seq_len(ncol(coloc_fmt)),
                       gridExpand = TRUE)

    # Highlight the top locus (max PP.H4) for readability - informational only,
    # it does not reach the PP.H4 > 0.8 threshold.
    top_row_s6 <- which.max(coloc_fmt$`PP.H4`) + 3
    openxlsx::addStyle(wb_s6, "S6 Colocalisation",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 9, textDecoration = "bold",
                         fgFill = "#FFF2CC"
                       ),
                       rows = top_row_s6, cols = seq_len(ncol(coloc_fmt)),
                       gridExpand = TRUE)

    openxlsx::setColWidths(
      wb_s6, "S6 Colocalisation",
      cols   = seq_len(ncol(coloc_fmt)),
      widths = c(14, 12, 16, 14, 14, 14, 9, 9, 9, 9, 9, 28)
    )

    footnote_row_s6 <- nrow(coloc_fmt) + 5
    footnote_text_s6 <- paste0(
      "Region start/end define the +/-500kb window around the index SNP used for each locus-level ",
      "colocalisation test. PP.H0, no association with either trait at this locus; PP.H1, association ",
      "with endometriosis only; PP.H2, association with placenta praevia only; PP.H3, both traits ",
      "associated but with distinct causal variants; PP.H4, both traits associated with a shared causal ",
      "variant. Predominant hypothesis is the posterior category with the highest probability at each ",
      "locus. PP.H4 > 0.8 is the threshold conventionally used to indicate strong support for ",
      "colocalisation; no locus reached this threshold (maximum PP.H4 = 0.33, at the locus tagged by ",
      "rs66683298, highlighted), consistent with a polygenic mechanism distributed across multiple loci ",
      "of modest individual effect and not contradicting the overall Mendelian randomisation estimate ",
      "for genetic liability to endometriosis on placenta praevia."
    )
    openxlsx::writeData(wb_s6, "S6 Colocalisation",
                        x = footnote_text_s6, startRow = footnote_row_s6, startCol = 1)
    openxlsx::addStyle(wb_s6, "S6 Colocalisation",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 8, fontColour = "#595959",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = footnote_row_s6, cols = 1)
    openxlsx::setRowHeights(wb_s6, "S6 Colocalisation",
                            rows = footnote_row_s6, heights = 90)
    openxlsx::mergeCells(wb_s6, "S6 Colocalisation",
                         cols = 1:ncol(coloc_fmt), rows = footnote_row_s6)

    openxlsx::saveWorkbook(
      wb_s6,
      file.path(tables_dir, "Supplementary_Table_S6_colocalisation_endo_praevia.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S6 (Excel) saved.")
  }

  message("Supplementary Table S6 (colocalisation) done.")
}

###############################################################################
# SUPPLEMENTARY TABLE S7 - KOLLER SENSITIVITY MR ACROSS 30 OUTCOMES
#
#   [R3.3] Reviewer #3 ("similarities or discrepancies... with the one used
#   in their own analysis"): IVW sensitivity analysis using the Koller et al.
#   (2026) endometriosis instrument (86 EUR loci, Supplementary Table 4) across all 30 main
#   outcomes, with direction of effect compared against the primary
#   Rahmioglu et al. (2023) analysis (originally 7/30, extended once the
#   MR-PREG collaboration re-extracted the remaining 23 outcomes at the
#   Koller SNPs).
#   [EDIT] Rahmioglu OR/CI/P/q columns added alongside the Koller ones (not
#   just the categorical Concordant/Discordant column), so the primary vs.
#   sensitivity comparison is directly inspectable in one table rather than
#   requiring a cross-reference to Table 3. Mirrors the ad hoc comparison
#   already prepared as results/Koller_vs_Rahmioglu_publication_table.xlsx,
#   now folded into the main S7 pipeline instead of living as a separate
#   untracked file. ---
#   Sources: results/koller_ivw_results.csv (script 05.1),
#            results/ivw_results.csv (script 04).
###############################################################################

koller_ivw_path_s7 <- file.path(results_dir, "koller_ivw_results.csv")
rahm_ivw_path_s7   <- file.path(results_dir, "ivw_results.csv")

if (!file.exists(koller_ivw_path_s7) || !file.exists(rahm_ivw_path_s7)) {
  warning("koller_ivw_results.csv or ivw_results.csv not found - skipping Table S7.")
} else {

  koller_ivw_s7 <- readr::read_csv(koller_ivw_path_s7, show_col_types = FALSE) %>%
    dplyr::filter(outcome %in% vars_keep, method == "Inverse variance weighted")

  rahm_ivw_s7 <- readr::read_csv(rahm_ivw_path_s7, show_col_types = FALSE) %>%
    dplyr::filter(outcome %in% vars_keep, method == "Inverse variance weighted") %>%
    dplyr::select(
      outcome,
      nsnp_rahmioglu = nsnp,
      b_rahmioglu    = b,
      se_rahmioglu   = se,
      pval_rahmioglu = pval,
      qval_rahmioglu = qval
    )

  Table_S7 <- koller_ivw_s7 %>%
    dplyr::left_join(rahm_ivw_s7, by = "outcome") %>%
    dplyr::left_join(
      domain_meta %>% dplyr::mutate(outcome_id = outcome),
      by = c("outcome" = "outcome_id")
    ) %>%
    dplyr::mutate(
      Outcome = dplyr::recode(outcome, !!!outcome_labels, .default = outcome),
      Domain  = factor(domain, levels = domain_levels_t3),
      Source  = label_source(outcome),
      type    = dplyr::if_else(outcome %in% continuous_outcomes, "continuous", "binary"),
      `Rahmioglu IVW OR or beta` = dplyr::if_else(
        type == "binary",
        sprintf("%.2f", exp(b_rahmioglu)),
        sprintf("%.3f", b_rahmioglu)
      ),
      `Rahmioglu 95% CI` = dplyr::if_else(
        type == "binary",
        sprintf("%.2f-%.2f", exp(b_rahmioglu - 1.96 * se_rahmioglu), exp(b_rahmioglu + 1.96 * se_rahmioglu)),
        sprintf("%.3f-%.3f", b_rahmioglu - 1.96 * se_rahmioglu, b_rahmioglu + 1.96 * se_rahmioglu)
      ),
      `Rahmioglu P-value`     = fmt_p_s3(pval_rahmioglu),
      `Rahmioglu FDR q-value` = sprintf("%.3f", qval_rahmioglu),
      `Koller IVW OR or beta` = dplyr::if_else(
        type == "binary",
        sprintf("%.2f", exp(b)),
        sprintf("%.3f", b)
      ),
      `Koller 95% CI` = dplyr::if_else(
        type == "binary",
        sprintf("%.2f-%.2f", exp(b - 1.96 * se), exp(b + 1.96 * se)),
        sprintf("%.3f-%.3f", b - 1.96 * se, b + 1.96 * se)
      ),
      `Koller P-value`     = fmt_p_s3(pval),
      `Koller FDR q-value` = sprintf("%.3f", qval),
      `Direction compared with the primary Rahmioglu analysis` = dplyr::case_when(
        is.na(b_rahmioglu)          ~ NA_character_,
        sign(b) == sign(b_rahmioglu) ~ "Concordant",
        TRUE                         ~ "Discordant"
      )
    ) %>%
    dplyr::arrange(Domain, order_within) %>%
    dplyr::select(
      Outcome,
      Source,
      `Rahmioglu N SNPs` = nsnp_rahmioglu,
      `Rahmioglu IVW OR or beta`,
      `Rahmioglu 95% CI`,
      `Rahmioglu P-value`,
      `Rahmioglu FDR q-value`,
      `Koller N SNPs` = nsnp,
      `Koller IVW OR or beta`,
      `Koller 95% CI`,
      `Koller P-value`,
      `Koller FDR q-value`,
      `Direction compared with the primary Rahmioglu analysis`
    )

  # ── CSV ──────────────────────────────────────────────────────────────────────
  readr::write_csv(
    Table_S7,
    file.path(tables_dir, "Supplementary_Table_S7_koller_sensitivity_30_outcomes.csv")
  )
  log_info("Supplementary Table S7 (CSV) saved.")

  # ── Excel ────────────────────────────────────────────────────────────────────
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    wb_s7 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s7, "S7 Koller sensitivity")

    caption_text_s7 <- paste0(
      "Supplementary Table S7. Sensitivity Mendelian randomization analyses using the Koller et al. ",
      "(2026) endometriosis instrument across 30 pregnancy and perinatal outcomes, alongside the ",
      "corresponding primary (Rahmioglu et al. 2023) estimates for direct comparison."
    )
    openxlsx::writeData(wb_s7, "S7 Koller sensitivity",
                        x = caption_text_s7, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb_s7, "S7 Koller sensitivity",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         wrapText = TRUE
                       ),
                       rows = 1, cols = 1)
    openxlsx::mergeCells(wb_s7, "S7 Koller sensitivity", cols = 1:ncol(Table_S7), rows = 1)

    openxlsx::writeDataTable(
      wb_s7, "S7 Koller sensitivity",
      x          = Table_S7,
      startRow   = 3, startCol = 1,
      tableStyle = "TableStyleLight9",
      withFilter = TRUE
    )

    openxlsx::addStyle(wb_s7, "S7 Koller sensitivity",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = 3, cols = seq_len(ncol(Table_S7)),
                       gridExpand = TRUE)

    # Highlight rows surviving FDR correction (q < 0.05, Koller instrument)
    q_s7 <- suppressWarnings(as.numeric(Table_S7$`Koller FDR q-value`))
    sig_rows_s7 <- which(!is.na(q_s7) & q_s7 < 0.05) + 3
    if (length(sig_rows_s7) > 0) {
      openxlsx::addStyle(wb_s7, "S7 Koller sensitivity",
                         style = openxlsx::createStyle(
                           fontName = "Arial", fontSize = 9, textDecoration = "bold",
                           fgFill = "#FFF2CC"
                         ),
                         rows = sig_rows_s7, cols = seq_len(ncol(Table_S7)),
                         gridExpand = TRUE)
    }

    openxlsx::setColWidths(
      wb_s7, "S7 Koller sensitivity",
      cols   = seq_len(ncol(Table_S7)),
      widths = c(34, 16, 12, 14, 12, 12, 14, 12, 14, 12, 12, 14, 30)
    )

    footnote_row_s7 <- nrow(Table_S7) + 5
    footnote_text_s7 <- paste0(
      "IVW estimates only. Rahmioglu columns reproduce the primary analysis (Table 3) for direct ",
      "comparison; Rahmioglu instrument: 41 SNPs (Rahmioglu et al. 2023). Koller instrument: 86 EUR ",
      "loci (Koller et al. 2026, Supplementary Table 4); SNP coverage per outcome ranged 68-82/86 ",
      "depending on outcome source and harmonisation. Effect scale is OR for binary outcomes and beta ",
      "for continuous outcomes (gestational age, birthweight z-score), consistently for both ",
      "instruments. Each FDR q-value is the Benjamini-Hochberg-adjusted p-value within that ",
      "instrument's own set of 30 outcomes (Rahmioglu and Koller corrected separately); rows are ",
      "shaded where the Koller q < 0.05. Placenta praevia is the only outcome surviving FDR correction ",
      "under either instrument; several other outcomes are nominally significant (p < 0.05) under the ",
      "Koller instrument but do not survive FDR correction and are reported as hypothesis-generating ",
      "only. Direction compared with the primary Rahmioglu analysis indicates whether the ",
      "Koller-instrument point estimate falls on the same side of the null (Concordant) or the ",
      "opposite side (Discordant) as the corresponding primary IVW estimate; occasional discordance is ",
      "confined to outcomes that are non-significant and imprecise under both instruments and should ",
      "not be interpreted as a true conflict between the two analyses."
    )
    openxlsx::writeData(wb_s7, "S7 Koller sensitivity",
                        x = footnote_text_s7, startRow = footnote_row_s7, startCol = 1)
    openxlsx::addStyle(wb_s7, "S7 Koller sensitivity",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 8, fontColour = "#595959",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = footnote_row_s7, cols = 1)
    openxlsx::setRowHeights(wb_s7, "S7 Koller sensitivity",
                            rows = footnote_row_s7, heights = 130)
    openxlsx::mergeCells(wb_s7, "S7 Koller sensitivity",
                         cols = 1:ncol(Table_S7), rows = footnote_row_s7)

    openxlsx::saveWorkbook(
      wb_s7,
      file.path(tables_dir, "Supplementary_Table_S7_koller_sensitivity_30_outcomes.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S7 (Excel) saved.")
  }

  message("Supplementary Table S7 (Koller sensitivity, 30 outcomes) done.")
}
# SUPPLEMENTARY TABLE S8 - ADENOMYOSIS -> PLACENTA PRAEVIA (EXPLORATORY MR)
#
#   [R2.5] Reviewer #2, comment 5 (adenomyosis misclassification): exploratory
#   MR of genetic liability to adenomyosis (Koller et al. 2026) on placenta
#   praevia (FinnGen R12), to inform (not replace) discussion of possible
#   endometriosis/adenomyosis misclassification.
#   Sources: results/koller_adenomyosis_placenta_praevia_mr.csv,
#            results/koller_adenomyosis_heterogeneity.csv,
#            results/koller_adenomyosis_pleiotropy.csv (script 05.1).
###############################################################################

adeno_mr_path   <- file.path(results_dir, "koller_adenomyosis_placenta_praevia_mr.csv")
adeno_het_path  <- file.path(results_dir, "koller_adenomyosis_heterogeneity.csv")
adeno_plt_path  <- file.path(results_dir, "koller_adenomyosis_pleiotropy.csv")

if (!file.exists(adeno_mr_path)) {
  warning("koller_adenomyosis_placenta_praevia_mr.csv not found - skipping Table S8.")
} else {

  adeno_mr_raw <- readr::read_csv(adeno_mr_path, show_col_types = FALSE)

  method_labels_s8 <- c(
    "Inverse variance weighted" = "IVW",
    "MR Egger"                  = "MR-Egger",
    "Weighted median"           = "Weighted median"
  )

  # Cochran's Q (available for IVW and MR-Egger, not Weighted median)
  adeno_q <- if (file.exists(adeno_het_path)) {
    readr::read_csv(adeno_het_path, show_col_types = FALSE) %>%
      dplyr::mutate(`Cochran's Q p-value` = fmt_p_s3(Q_pval)) %>%
      dplyr::select(method, `Cochran's Q p-value`)
  } else {
    tibble::tibble(method = character(), `Cochran's Q p-value` = character())
  }

  # Egger intercept (a single diagnostic value, attached to the MR-Egger row)
  adeno_intercept_txt <- if (file.exists(adeno_plt_path)) {
    plt_s8 <- readr::read_csv(adeno_plt_path, show_col_types = FALSE)
    if (nrow(plt_s8) > 0) {
      sprintf("beta = %.3f (p = %s)",
              plt_s8$egger_intercept[1], fmt_p_s3(plt_s8$pval[1]))
    } else {
      NA_character_
    }
  } else {
    NA_character_
  }

  Table_S8 <- adeno_mr_raw %>%
    dplyr::left_join(adeno_q, by = "method") %>%
    dplyr::mutate(
      Method   = dplyr::recode(method, !!!method_labels_s8, .default = method),
      OR       = sprintf("%.2f", exp(b)),
      `95% CI` = sprintf("(%.2f-%.2f)", exp(b - 1.96 * se), exp(b + 1.96 * se)),
      `P-value` = fmt_p_s3(pval),
      `Cochran's Q p-value` = dplyr::coalesce(`Cochran's Q p-value`, "-"),
      `Egger intercept (beta), p-value` = dplyr::if_else(
        Method == "MR-Egger" & !is.na(adeno_intercept_txt),
        adeno_intercept_txt,
        "-"
      )
    ) %>%
    dplyr::select(
      `No. SNPs` = nsnp,
      Method,
      OR,
      `95% CI`,
      `P-value`,
      `Cochran's Q p-value`,
      `Egger intercept (beta), p-value`
    )

  # ── CSV ──────────────────────────────────────────────────────────────────────
  readr::write_csv(
    Table_S8,
    file.path(tables_dir, "Supplementary_Table_S8_adenomyosis_praevia_MR.csv")
  )
  log_info("Supplementary Table S8 (CSV) saved.")

  # ── Excel ────────────────────────────────────────────────────────────────────
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    wb_s8 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s8, "S8 Adenomyosis-praevia")

    caption_text_s8 <- paste0(
      "Supplementary Table S8. Exploratory Mendelian randomization analysis of genetic liability ",
      "to adenomyosis and placenta praevia."
    )
    openxlsx::writeData(wb_s8, "S8 Adenomyosis-praevia",
                        x = caption_text_s8, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb_s8, "S8 Adenomyosis-praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         wrapText = TRUE
                       ),
                       rows = 1, cols = 1)
    openxlsx::mergeCells(wb_s8, "S8 Adenomyosis-praevia", cols = 1:ncol(Table_S8), rows = 1)

    openxlsx::writeDataTable(
      wb_s8, "S8 Adenomyosis-praevia",
      x          = Table_S8,
      startRow   = 3, startCol = 1,
      tableStyle = "TableStyleLight9",
      withFilter = TRUE
    )

    openxlsx::addStyle(wb_s8, "S8 Adenomyosis-praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = 3, cols = seq_len(ncol(Table_S8)),
                       gridExpand = TRUE)

    openxlsx::setColWidths(
      wb_s8, "S8 Adenomyosis-praevia",
      cols   = seq_len(ncol(Table_S8)),
      widths = c(9, 16, 8, 14, 12, 16, 26)
    )

    footnote_row_s8 <- nrow(Table_S8) + 5
    footnote_text_s8 <- paste0(
      "Exploratory analysis performed to inform (not replace) discussion of possible endometriosis/",
      "adenomyosis misclassification (Reviewer #2, comment 5); it does not use the manuscript's primary ",
      "endometriosis instrument or outcome data. Genetic instrument: 6 independent, genome-wide ",
      "significant SNPs for adenomyosis (Koller et al. 2026, EUR meta-analysis), clumped and harmonised ",
      "against FinnGen R12 placenta praevia. Given the small number of instrument SNPs, this analysis is ",
      "clearly underpowered and should be interpreted as hypothesis-generating only, not as evidence for ",
      "or against a causal effect of adenomyosis liability on placenta praevia. ",
      "[!] Sample overlap: unlike the primary endometriosis-placenta praevia analysis, for which a ",
      "targeted MRlap/cross-trait LD score regression analysis against FinnGen R12 gave an intercept ",
      "compatible with zero (0.001, SE = 0.005; see response to Reviewer #4) - consistent with, but not ",
      "proof of, absence of sample overlap (an intercept near zero can also reflect low phenotypic ",
      "correlation between the traits even with some overlap) - no equivalent overlap ",
      "assessment was performed for the Koller et al. 2026 adenomyosis GWAS. As a large multi-biobank ",
      "European meta-analysis, this GWAS may include FinnGen among its contributing cohorts; if so, ",
      "FinnGen would then contribute to both the exposure and outcome samples here, which could bias ",
      "this estimate (e.g. toward the confounded/observational association) in a way not captured by ",
      "standard two-sample MR standard errors. This estimate should therefore be interpreted with ",
      "additional caution beyond the sample-size limitation alone. ",
      "Cochran's Q p-value tests heterogeneity of SNP-specific causal estimates (available for IVW and ",
      "MR-Egger only). The Egger intercept tests for directional horizontal pleiotropy and is reported ",
      "once, attached to the MR-Egger row; - indicates not applicable/not computed for that row."
    )
    openxlsx::writeData(wb_s8, "S8 Adenomyosis-praevia",
                        x = footnote_text_s8, startRow = footnote_row_s8, startCol = 1)
    openxlsx::addStyle(wb_s8, "S8 Adenomyosis-praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 8, fontColour = "#595959",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = footnote_row_s8, cols = 1)
    openxlsx::setRowHeights(wb_s8, "S8 Adenomyosis-praevia",
                            rows = footnote_row_s8, heights = 150)
    openxlsx::mergeCells(wb_s8, "S8 Adenomyosis-praevia",
                         cols = 1:ncol(Table_S8), rows = footnote_row_s8)

    openxlsx::saveWorkbook(
      wb_s8,
      file.path(tables_dir, "Supplementary_Table_S8_adenomyosis_praevia_MR.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S8 (Excel) saved.")
  }

  message("Supplementary Table S8 (adenomyosis exploratory MR) done.")
}

###############################################################################
# SUPPLEMENTARY TABLE S9 - ADENOMYOSIS -> ALL 30 OUTCOMES (IVW SCREEN)
#
#   [R2.5] Exploratory screen of the adenomyosis instrument (Koller et al.
#   2026) across all 30 main outcomes. IVW only (see Table S8 for the
#   pleiotropy-robust methods run on placenta praevia specifically).
#   Source: results/koller_adenomyosis_all_outcomes_ivw.csv (script 05.1).
###############################################################################

adeno_ivw30_path <- file.path(results_dir, "koller_adenomyosis_all_outcomes_ivw.csv")

if (!file.exists(adeno_ivw30_path)) {
  warning("koller_adenomyosis_all_outcomes_ivw.csv not found - skipping Table S9.")
} else {

  adeno_ivw30_raw <- readr::read_csv(adeno_ivw30_path, show_col_types = FALSE)

  Table_S9 <- adeno_ivw30_raw %>%
    dplyr::mutate(
      Outcome_label = dplyr::recode(outcome, !!!outcome_labels, .default = outcome),
      type          = dplyr::if_else(outcome %in% continuous_outcomes, "continuous", "binary")
    ) %>%
    dplyr::arrange(pval) %>%
    dplyr::transmute(
      Outcome      = Outcome_label,
      `No. SNPs`   = nsnp,
      `OR or beta` = dplyr::if_else(
        type == "binary",
        sprintf("%.2f", exp(b)),
        sprintf("%.3f", b)
      ),
      `95% CI`     = dplyr::if_else(
        type == "binary",
        sprintf("(%.2f-%.2f)", exp(b - 1.96 * se), exp(b + 1.96 * se)),
        sprintf("(%.3f to %.3f)", b - 1.96 * se, b + 1.96 * se)
      ),
      `P-value`  = fmt_p_s3(pval),
      `q (FDR)`  = sprintf("%.3f", qval)
    )

  # ── CSV ──────────────────────────────────────────────────────────────────────
  readr::write_csv(
    Table_S9,
    file.path(tables_dir, "Supplementary_Table_S9_adenomyosis_30_outcomes.csv")
  )
  log_info("Supplementary Table S9 (CSV) saved.")

  # ── Excel ────────────────────────────────────────────────────────────────────
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    wb_s9 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s9, "S9 Adenomyosis 30 outcomes")

    caption_text_s9 <- paste0(
      "Supplementary Table S9. Exploratory IVW screen of genetic liability to adenomyosis across ",
      "all 30 main outcomes."
    )
    openxlsx::writeData(wb_s9, "S9 Adenomyosis 30 outcomes",
                        x = caption_text_s9, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb_s9, "S9 Adenomyosis 30 outcomes",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         wrapText = TRUE
                       ),
                       rows = 1, cols = 1)
    openxlsx::mergeCells(wb_s9, "S9 Adenomyosis 30 outcomes", cols = 1:ncol(Table_S9), rows = 1)

    openxlsx::writeDataTable(
      wb_s9, "S9 Adenomyosis 30 outcomes",
      x          = Table_S9,
      startRow   = 3, startCol = 1,
      tableStyle = "TableStyleLight9",
      withFilter = TRUE
    )

    openxlsx::addStyle(wb_s9, "S9 Adenomyosis 30 outcomes",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = 3, cols = seq_len(ncol(Table_S9)),
                       gridExpand = TRUE)

    openxlsx::setColWidths(
      wb_s9, "S9 Adenomyosis 30 outcomes",
      cols   = seq_len(ncol(Table_S9)),
      widths = c(30, 9, 8, 14, 12, 10)
    )

    footnote_row_s9 <- nrow(Table_S9) + 5
    footnote_text_s9 <- paste0(
      "Exploratory screen using the same 6-SNP adenomyosis instrument as Table S8, harmonised against ",
      "all 30 main outcomes. Only the inverse-variance weighted method was calculated, given the ",
      "limited number of instruments and the exploratory nature of this analysis; pleiotropy-robust ",
      "methods were reserved for placenta praevia specifically (Table S8), the outcome of primary ",
      "interest for this misclassification check. No outcome survived Benjamini-Hochberg FDR ",
      "correction (q>0.6 throughout). As for Table S8, sample overlap between the Koller et al. 2026 ",
      "adenomyosis GWAS and the outcome GWAS could not be excluded and was not formally assessed."
    )
    openxlsx::writeData(wb_s9, "S9 Adenomyosis 30 outcomes",
                        x = footnote_text_s9, startRow = footnote_row_s9, startCol = 1)
    openxlsx::addStyle(wb_s9, "S9 Adenomyosis 30 outcomes",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 8, fontColour = "#595959",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = footnote_row_s9, cols = 1)
    openxlsx::setRowHeights(wb_s9, "S9 Adenomyosis 30 outcomes",
                            rows = footnote_row_s9, heights = 120)
    openxlsx::mergeCells(wb_s9, "S9 Adenomyosis 30 outcomes",
                         cols = 1:ncol(Table_S9), rows = footnote_row_s9)

    openxlsx::saveWorkbook(
      wb_s9,
      file.path(tables_dir, "Supplementary_Table_S9_adenomyosis_30_outcomes.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S9 (Excel) saved.")
  }

  message("Supplementary Table S9 (adenomyosis 30-outcome IVW screen) done.")
}

###############################################################################
# SUPPLEMENTARY TABLE S10 - MRLAP OVERLAP-AWARE SENSITIVITY ANALYSIS
#                           (ENDOMETRIOSIS x PLACENTA PRAEVIA)
#
#   [R4.2] Reviewer #4 ("Report the overlapping fraction per outcome and
#   overlap-aware estimates (MRlap); the bivariate LDSC intercept gives this
#   as a by-product of point 5."): targeted MRlap sensitivity analysis for
#   the study's single positive finding (endometriosis -> placenta praevia).
#   Source: results/mrlap_praevia_results.csv (script 05.5).
#   NOTE: observed/corrected effects are on MRlap's internal Z/sqrt(N)-
#   standardised scale, NOT log-odds - not directly comparable to the
#   primary IVW OR, which is shown alongside for reference only.
###############################################################################

mrlap_path <- file.path(results_dir, "mrlap_praevia_results.csv")

if (!file.exists(mrlap_path)) {
  warning("mrlap_praevia_results.csv not found - skipping Table S10.")
} else {

  mrlap_raw <- readr::read_csv(mrlap_path, show_col_types = FALSE)
  mv <- setNames(as.character(mrlap_raw$value), mrlap_raw$metric)
  num <- function(x) as.numeric(mv[[x]])

  # Primary IVW estimate for placenta praevia, for reference only (not on
  # the same scale as the MRlap observed/corrected effects above).
  praevia_ivw <- readr::read_csv(file.path(results_dir, "ivw_results.csv"), show_col_types = FALSE) %>%
    dplyr::filter(outcome == "finngen_R12_O15_PLAC_PRAEVIA", method == "Inverse variance weighted")

  Table_S10 <- tibble::tribble(
    ~Metric,                                                  ~Value,
    "MRlap package version",                                  mv[["mrlap_package_version"]],
    "Exposure GWAS (accession, total N)",                     "Rahmioglu et al. 2023, GCST90205183 (475,160)",
    "Outcome GWAS (total N)",                                 "FinnGen R12 placenta praevia (223,001)",
    "Number of independent instruments",                      mv[["n_instruments_used"]],
    "Cross-trait LDSC intercept (SE)",                        sprintf("%.1e (SE %.1e)", num("crosstrait_intercept"), num("crosstrait_intercept_se")),
    "h2 exposure (SE)",                                  sprintf("%.4f (SE %.4f)", num("h2_exposure"), num("h2_exposure_se")),
    "h2 outcome (SE)",                                   sprintf("%.5f (SE %.5f)", num("h2_outcome"), num("h2_outcome_se")),
    "Genetic correlation (rg), exposure-outcome",             sprintf("%.3f", num("rg_exposure_outcome")),
    "Observed effect, standardised scale (SE), P",            sprintf("beta = %.3f (SE %.3f), P = %.1e", num("observed_effect_std"), num("observed_effect_se_std"), num("observed_effect_pval")),
    "Corrected effect, standardised scale (SE), P",           sprintf("beta = %.3f (SE %.3f), P = %.1e", num("corrected_effect_std"), num("corrected_effect_se_std"), num("corrected_effect_pval")),
    "P for difference (observed vs. corrected)",              sprintf("%.1e", num("p_difference")),
    "MR-Egger intercept P (directional pleiotropy)",          sprintf("%.3f", num("egger_intercept_pval")),
    "Warnings raised during MRlap run",                       sprintf("%s (all R package-loading namespace conflicts, e.g. 'replacing previous import ... when loading GenomicSEM'; none relate to model convergence, instrument strength, or data validity - verified individually, see results/mrlap_praevia_warnings.txt and Methods)", mv[["n_warnings"]]),
    "Primary IVW estimate, for reference (not the same scale)", sprintf(
      "OR = %.2f (95%% CI %.2f-%.2f), P = %.1e",
      exp(praevia_ivw$b), exp(praevia_ivw$b - 1.96 * praevia_ivw$se),
      exp(praevia_ivw$b + 1.96 * praevia_ivw$se), praevia_ivw$pval
    )
  )

  # ── CSV ──────────────────────────────────────────────────────────────────────
  readr::write_csv(
    Table_S10,
    file.path(tables_dir, "Supplementary_Table_S10_mrlap_praevia.csv")
  )
  log_info("Supplementary Table S10 (CSV) saved.")

  # ── Excel ────────────────────────────────────────────────────────────────────
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    wb_s10 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s10, "S10 MRlap praevia")

    caption_text_s10 <- paste0(
      "Supplementary Table S10. MRlap overlap-aware sensitivity analysis for genetic liability to ",
      "endometriosis and placenta praevia."
    )
    openxlsx::writeData(wb_s10, "S10 MRlap praevia",
                        x = caption_text_s10, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb_s10, "S10 MRlap praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         wrapText = TRUE
                       ),
                       rows = 1, cols = 1)
    openxlsx::mergeCells(wb_s10, "S10 MRlap praevia", cols = 1:ncol(Table_S10), rows = 1)

    openxlsx::writeDataTable(
      wb_s10, "S10 MRlap praevia",
      x          = Table_S10,
      startRow   = 3, startCol = 1,
      tableStyle = "TableStyleLight9",
      withFilter = FALSE
    )

    openxlsx::addStyle(wb_s10, "S10 MRlap praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = 3, cols = seq_len(ncol(Table_S10)),
                       gridExpand = TRUE)

    # Highlight the observed/corrected effect + P-difference rows (core result)
    core_rows <- which(Table_S10$Metric %in% c(
      "Observed effect, standardised scale (SE), P",
      "Corrected effect, standardised scale (SE), P",
      "P for difference (observed vs. corrected)"
    )) + 3
    openxlsx::addStyle(wb_s10, "S10 MRlap praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 9, textDecoration = "bold",
                         fgFill = "#FFF2CC"
                       ),
                       rows = core_rows, cols = seq_len(ncol(Table_S10)),
                       gridExpand = TRUE)

    openxlsx::setColWidths(
      wb_s10, "S10 MRlap praevia",
      cols   = seq_len(ncol(Table_S10)),
      widths = c(44, 46)
    )

    footnote_row_s10 <- nrow(Table_S10) + 5
    footnote_text_s10 <- paste0(
      "Both GWAS were supplied to MRlap with their total (not per-SNP effective) sample size, as ",
      "required for case-control data. The observed and corrected effects are expressed on MRlap's ",
      "internal Z/sqrt(N)-standardised effect-size scale (see Mounier & Kutalik, 2023), not on the ",
      "raw log-odds scale, and must not be exponentiated or otherwise compared numerically to the ",
      "primary IVW odds ratio shown for reference only. The cross-trait LDSC intercept reflects genetic ",
      "covariance from sample overlap and/or phenotypic correlation between the traits in the ",
      "overlapping samples, and cannot separate the two; the observed-vs-corrected comparison further ",
      "jointly reflects sample overlap, weak-instrument bias, and winner's-curse-related bias and ",
      "cannot isolate overlap alone. A significant P for difference indicates that the ",
      "corrected estimate should be considered more reliable than the observed estimate, per MRlap's ",
      "own recommended interpretation. The 23 warnings are R package-loading namespace conflicts ",
      "(MRlap calls GenomicSEM internally), not convergence or data-validity warnings. Table S11 uses ",
      "the same software and exposure N (475,160); its h2(exposure)=0.0163 and r_g=0.374 corroborate ",
      "the h2=0.0163 and r_g=0.385 reported here."
    )
    openxlsx::writeData(wb_s10, "S10 MRlap praevia",
                        x = footnote_text_s10, startRow = footnote_row_s10, startCol = 1)
    openxlsx::addStyle(wb_s10, "S10 MRlap praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 8, fontColour = "#595959",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = footnote_row_s10, cols = 1)
    openxlsx::setRowHeights(wb_s10, "S10 MRlap praevia",
                            rows = footnote_row_s10, heights = 190)
    openxlsx::mergeCells(wb_s10, "S10 MRlap praevia",
                         cols = 1:ncol(Table_S10), rows = footnote_row_s10)

    openxlsx::saveWorkbook(
      wb_s10,
      file.path(tables_dir, "Supplementary_Table_S10_mrlap_praevia.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S10 (Excel) saved.")
  }

  message("Supplementary Table S10 (MRlap praevia) done.")
}

###############################################################################
# SUPPLEMENTARY TABLE S11 - TARGETED CROSS-TRAIT LDSC (r_g) ANALYSIS
#                            (ENDOMETRIOSIS x PLACENTA PRAEVIA)
#
#   [R2.2]/[R4.5] Genetic correlation between endometriosis and placenta
#   praevia via LDSC. Source: script 05.6's "primary" variant (corrected N,
#   QC-filtered); other variants summarised as a min-max range here, full
#   detail in bivariate_rg.csv. r_g is non-directional and complementary to
#   the MR estimate - not a validation of it - and distinct from the
#   cross-trait intercept (sample-overlap diagnostic, reported in Table S10).
###############################################################################

# [SYNC] This dated folder name is also hardcoded in
# scripts/05.6_ldsc_targeted_endo_praevia_endoMR-PREG.R (variable `out_dir`).
# If script 05.6 is ever re-run on a different date, update BOTH places -
# this block will otherwise keep silently reading the stale 2026-07-28
# folder without any error, since it only checks file.exists(), not recency.
ldsc_dir <- file.path(results_dir, "ldsc_endometriosis_placenta_praevia_20260728")
univar_path <- file.path(ldsc_dir, "univariate_heritability.csv")
bivar_path  <- file.path(ldsc_dir, "bivariate_rg.csv")

if (!file.exists(univar_path) || !file.exists(bivar_path)) {
  warning("LDSC endo-praevia results not found - skipping Table S11.")
} else {

  univar_raw <- readr::read_csv(univar_path, show_col_types = FALSE)
  bivar_raw  <- readr::read_csv(bivar_path, show_col_types = FALSE)

  primary_label <- "primary (N=475,160 + ambig/MHC/indel/dup filter)"
  u_endo   <- dplyr::filter(univar_raw, variant == primary_label, trait == "endometriosis")
  u_praev  <- dplyr::filter(univar_raw, variant == primary_label, trait == "placenta_praevia")
  b_primary <- dplyr::filter(bivar_raw, variant == primary_label)

  # Sensitivity rg range across all three variants (original / n_corrected / primary)
  rg_range <- range(bivar_raw$rg)

  Table_S11 <- tibble::tribble(
    ~Metric,                                              ~Value,
    "LDSC implementation",                                "GenomicSEM v0.0.5 (R reimplementation of Bulik-Sullivan et al. 2015)",
    "Exposure GWAS (accession, total N)",                 "Rahmioglu et al. 2023 endometriosis, GCST90205183 (N=475,160)",
    "Outcome GWAS (accession, total N)",                  "FinnGen R12 O15_PLAC_PRAEVIA (N=223,001)",
    "LD reference",                                       "1000 Genomes European (eur_w_ld_chr), HapMap3 SNPs",
    "N SNPs, endometriosis (post-QC)",                    format(u_endo$n_snps, big.mark = ","),
    "N SNPs, placenta praevia (post-QC)",                 format(u_praev$n_snps, big.mark = ","),
    "N SNPs, cross-trait overlap",                        format(b_primary$n_snps_overlap, big.mark = ","),
    "SNP-h2 endometriosis, observed scale (SE), Z, P",    sprintf("%.4f (SE %.4f), Z=%.1f, P=%.1e", u_endo$h2_obs, u_endo$h2_obs_se, u_endo$h2_z, u_endo$h2_p),
    "SNP-h2 placenta praevia, observed scale (SE), Z, P", sprintf("%.4f (SE %.4f), Z=%.2f, P=%.2f", u_praev$h2_obs, u_praev$h2_obs_se, u_praev$h2_z, u_praev$h2_p),
    "Cross-trait LDSC intercept (SE)",                    sprintf("%.4f (SE %.3f)", b_primary$crosstrait_intercept, b_primary$crosstrait_intercept_se),
    "Genetic correlation (r_g), endo-praevia (SE), P",    sprintf("%.3f (SE %.3f), P=%.3f", b_primary$rg, b_primary$rg_se, b_primary$rg_p),
    "Sensitivity: r_g range across N/QC specifications",  sprintf("%.3f - %.3f across 3 analytic specifications (see Methods)", rg_range[1], rg_range[2])
  )

  # ── CSV ──────────────────────────────────────────────────────────────────────
  readr::write_csv(
    Table_S11,
    file.path(tables_dir, "Supplementary_Table_S11_ldsc_endo_praevia.csv")
  )
  log_info("Supplementary Table S11 (CSV) saved.")

  # ── Excel ────────────────────────────────────────────────────────────────────
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    wb_s11 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s11, "S11 LDSC endo praevia")

    caption_text_s11 <- paste0(
      "Supplementary Table S11. Targeted cross-trait LD Score Regression analysis of the genome-wide ",
      "genetic correlation between endometriosis and placenta praevia."
    )
    openxlsx::writeData(wb_s11, "S11 LDSC endo praevia",
                        x = caption_text_s11, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb_s11, "S11 LDSC endo praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         wrapText = TRUE
                       ),
                       rows = 1, cols = 1)
    openxlsx::mergeCells(wb_s11, "S11 LDSC endo praevia", cols = 1:ncol(Table_S11), rows = 1)

    openxlsx::writeDataTable(
      wb_s11, "S11 LDSC endo praevia",
      x          = Table_S11,
      startRow   = 3, startCol = 1,
      tableStyle = "TableStyleLight9",
      withFilter = FALSE
    )

    openxlsx::addStyle(wb_s11, "S11 LDSC endo praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = 3, cols = seq_len(ncol(Table_S11)),
                       gridExpand = TRUE)

    # Highlight the core rg/intercept result rows
    core_rows_s11 <- which(Table_S11$Metric %in% c(
      "Cross-trait LDSC intercept (SE)",
      "Genetic correlation (r_g), endo-praevia (SE), P"
    )) + 3
    openxlsx::addStyle(wb_s11, "S11 LDSC endo praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 9, textDecoration = "bold",
                         fgFill = "#FFF2CC", wrapText = TRUE, valign = "top"
                       ),
                       rows = core_rows_s11, cols = seq_len(ncol(Table_S11)),
                       gridExpand = TRUE)

    openxlsx::setColWidths(
      wb_s11, "S11 LDSC endo praevia",
      cols   = seq_len(ncol(Table_S11)),
      widths = c(44, 58)
    )

    footnote_row_s11 <- nrow(Table_S11) + 5
    footnote_text_s11 <- paste0(
      "HapMap3 SNPs, European 1000 Genomes LD scores. Exposure N=475,160 is GCST90205183 (the file ",
      "used here), not the 762,600 figure from GCST90258638 used in the primary MR analysis. rg is ",
      "stable across N/QC specifications (0.372-0.374); placenta praevia h2 is imprecisely estimated ",
      "given its case count, which limits precision but not validity. The cross-trait intercept is a ",
      "sample-overlap diagnostic, not rg. Table S10 (MRlap) uses the same N and software and gives a ",
      "corroborating h2=0.0163, rg=0.385."
    )
    openxlsx::writeData(wb_s11, "S11 LDSC endo praevia",
                        x = footnote_text_s11, startRow = footnote_row_s11, startCol = 1)
    openxlsx::addStyle(wb_s11, "S11 LDSC endo praevia",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 8, fontColour = "#595959",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = footnote_row_s11, cols = 1)
    openxlsx::setRowHeights(wb_s11, "S11 LDSC endo praevia",
                            rows = footnote_row_s11, heights = 150)
    openxlsx::mergeCells(wb_s11, "S11 LDSC endo praevia",
                         cols = 1:ncol(Table_S11), rows = footnote_row_s11)

    openxlsx::saveWorkbook(
      wb_s11,
      file.path(tables_dir, "Supplementary_Table_S11_ldsc_endo_praevia.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S11 (Excel) saved.")
  }

  message("Supplementary Table S11 (LDSC endo-praevia) done.")
}

###############################################################################
# SUPPLEMENTARY TABLE S12 - CHARACTERISTICS OF THE SENSITIVITY INSTRUMENTS
#                            (KOLLER ET AL. 2026, ENDOMETRIOSIS + ADENOMYOSIS)
#
#   Parallel to Supplementary Table 2 (primary Rahmioglu instrument), for the
#   two Koller-derived instruments used in the revision follow-up analyses
#   (Sections 05.1 / Supplementary Tables S7-S9). F-statistic computed as
#   (beta/se)^2, consistent with the summary values reported in the Results
#   ("Genetic instruments" paragraph: Koller endometriosis mean F=51.8,
#   Koller adenomyosis mean F=42.6) - no per-SNP effective N is available for
#   the endometriosis instrument (published Beta/SE used directly, see
#   05.1_koller_sensitivity_endoMR-PREG.R REVISION LOG), so R2/F are not
#   computed via the N-dependent formula used in Table 2.
#   Sources: results/Koller_Endometriosis_Table4_EUR_snps.tsv,
#            results/Koller_Adenomyosis_clumped_snps.tsv (script 05.1).
###############################################################################

message("\n=== Creating Supplementary Table S12: Koller instrument characteristics ===")

koller_endo_snp_file  <- file.path(results_dir, "Koller_Endometriosis_Table4_EUR_snps.tsv")
koller_adeno_snp_file <- file.path(results_dir, "Koller_Adenomyosis_clumped_snps.tsv")

if (!file.exists(koller_endo_snp_file) || !file.exists(koller_adeno_snp_file)) {
  warning("Koller instrument SNP files not found - skipping Table S12.")
} else {

  fmt_koller_instrument <- function(path, instrument_label) {
    data.table::fread(path) %>%
      dplyr::transmute(
        Instrument      = instrument_label,
        SNP             = SNP,
        Beta            = round(beta.exposure, 4),
        SE              = round(se.exposure, 4),
        `P-value`       = formatC(pval.exposure, format = "e", digits = 2),
        `Effect allele` = effect_allele.exposure,
        `Other allele`  = other_allele.exposure,
        EAF             = round(eaf.exposure, 4),
        F               = round((beta.exposure / se.exposure)^2, 1)
      )
  }

  Table_S12 <- dplyr::bind_rows(
    fmt_koller_instrument(koller_endo_snp_file,  "Endometriosis (Koller et al. 2026, Supplementary Table 4, 86 EUR loci)"),
    fmt_koller_instrument(koller_adeno_snp_file, "Adenomyosis (Koller et al. 2026)")
  )

  readr::write_csv(
    Table_S12,
    file.path(tables_dir, "Supplementary_Table_S12_koller_instrument_characteristics.csv")
  )
  log_info("Supplementary Table S12 (CSV) saved.")

  if (requireNamespace("openxlsx", quietly = TRUE)) {

    wb_s12 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s12, "S12 Koller instruments")

    caption_text_s12 <- paste0(
      "Supplementary Table S12. Characteristics of the sensitivity endometriosis and adenomyosis ",
      "instruments (Koller et al. 2026) used in the revision follow-up analyses."
    )
    openxlsx::writeData(wb_s12, "S12 Koller instruments",
                        x = caption_text_s12, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb_s12, "S12 Koller instruments",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         wrapText = TRUE
                       ),
                       rows = 1, cols = 1)
    openxlsx::mergeCells(wb_s12, "S12 Koller instruments", cols = 1:ncol(Table_S12), rows = 1)

    openxlsx::writeDataTable(
      wb_s12, "S12 Koller instruments",
      x          = Table_S12,
      startRow   = 3, startCol = 1,
      tableStyle = "TableStyleLight9",
      withFilter = TRUE
    )

    openxlsx::addStyle(wb_s12, "S12 Koller instruments",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = 3, cols = seq_len(ncol(Table_S12)),
                       gridExpand = TRUE)

    openxlsx::setColWidths(
      wb_s12, "S12 Koller instruments",
      cols   = seq_len(ncol(Table_S12)),
      widths = c(46, 14, 10, 10, 12, 12, 12, 10, 10)
    )

    footnote_row_s12 <- nrow(Table_S12) + 5
    footnote_text_s12 <- paste0(
      "SNPs for the endometriosis instrument (86 EUR loci) use the published Beta/SE/EAF from Koller ",
      "et al. (2026) Supplementary Table 4 directly, not re-derived from Z-scores; no per-SNP effective ",
      "sample size is available for this instrument, so F is computed as (Beta/SE)² rather than via ",
      "the R²-based formula used in Supplementary Table 2. The adenomyosis instrument (6 SNPs) uses ",
      "the same (Beta/SE)² definition for consistency between the two instruments in this table, ",
      "although a per-SNP effective sample size is available for it (see Methods, Section 1.5, and ",
      "05.1_koller_sensitivity_endoMR-PREG.R). Neither instrument was re-clumped for this analysis: the ",
      "endometriosis loci are already LD-independent in the source publication, and the adenomyosis ",
      "instrument was clumped upstream (r²<0.001, 10,000kb window, 1000G EUR) exactly as for the ",
      "primary Rahmioglu instrument (Supplementary Table 2)."
    )
    openxlsx::writeData(wb_s12, "S12 Koller instruments",
                        x = footnote_text_s12, startRow = footnote_row_s12, startCol = 1)
    openxlsx::addStyle(wb_s12, "S12 Koller instruments",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 8, fontColour = "#595959",
                         wrapText = TRUE, valign = "top"
                       ),
                       rows = footnote_row_s12, cols = 1)
    openxlsx::setRowHeights(wb_s12, "S12 Koller instruments",
                            rows = footnote_row_s12, heights = 110)
    openxlsx::mergeCells(wb_s12, "S12 Koller instruments",
                         cols = 1:ncol(Table_S12), rows = footnote_row_s12)

    openxlsx::saveWorkbook(
      wb_s12,
      file.path(tables_dir, "Supplementary_Table_S12_koller_instrument_characteristics.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S12 (Excel) saved.")
  }

  message("Supplementary Table S12 (Koller instrument characteristics) done.")
}

