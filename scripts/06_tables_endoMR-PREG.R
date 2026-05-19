#!/usr/bin/env Rscript

###############################################################################
# 06_tables_endoMR-PREG.R
#
# Generate main and supplementary tables for:
#   Genetic liability to endometriosis → pregnancy outcomes (endoMR-PREG)
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

outcome_labels <- c(
    # Bleeding (4)
    Antepartum_bleeding                       = "Antepartum bleeding",
    Postpartum_hemorrhage                     = "Postpartum haemorrhage (any)",
    Postpartum_hemorrhage_due_to_atony        = "PPH due to uterine atony",
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
    anaemia_preg_all                          = "Pregnancy anaemia",
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
# 3) MAIN TABLE 1 — DESCRIPTION OF GWAS DATASETS
###############################################################################

log_info("Creating Table 1: description of GWAS datasets...")

table1_gwas <- tibble::tribble(
  ~Study,                      ~Year, ~Dataset_or_consortium,                                            ~Phenotype_group,                                                      ~Ancestry,                                             ~Sample_size,                                             ~Sex,           ~Main_adjustments,                                                   ~Notes_or_accession,
  "Rahmioglu et al.",          2023,  "International Endometriosis Genetics Consortium + UK Biobank",    "Endometriosis (overall and subtypes)",                              "Predominantly European (~98%) + Japanese (~2%)",     "60,674 cases; 701,926 controls",                       "Female/mixed", "Age, batch, principal components, study-specific covariates", "Primary exposure GWAS; clumped instruments used in endoMR-PREG",
  "FinnGen R12",               2024,  "FinnGen (release 12)",                                            "Pregnancy and fertility ICD-10 phenotypes",                          "Finnish (European)",                                 "Varies by phenotype",                                   "Female",       "Age, batch, principal components",                       "Used for placental phenotypes (O43–O45) and related outcomes",
  "McBride et al. (MR-PREG)",  2026,  "MR-PREG collaboration (ALSPAC, BiB, MoBa, UKB, FinnGen + GWAS)",  "Adverse pregnancy and perinatal outcomes (binary and continuous)",   "Predominantly European",                              "Up to 678,001 women (outcome-specific)",               "Female",       "Age, study/centre, principal components",              "Core source for hypertensive disorders, GDM, PTB, SGA/LGA, CS, PROM, stillbirth, Apgar, NICU",
  "Westergaard et al.",        2024,  "Nordic registry-based GWAS (6 cohorts)",                           "Bleeding in pregnancy and postpartum haemorrhage (overall, subtypes)", "Northern European (registry-based)",               "Up to 331,792 women; outcome-specific case counts",     "Female",       "Age, parity, calendar year, cohort",                   "Provides GWAS for antepartum bleeding and PPH subtypes used in endoMR-PREG"
)

table1_file <- file.path(tables_dir, "Table1_GWAS_sources.csv")
write.csv(table1_gwas, table1_file, row.names = FALSE)
log_info("Table 1 (GWAS sources) saved: ", table1_file)

table1_gwas_gt <- table1_gwas %>%
  gt::gt() %>%
  gt::tab_header(
    title = gt::md("**Table 1. Description of exposure and outcome GWAS datasets used in endoMR-PREG**")
  ) %>%
  gt::cols_label(
    Study                = "Study",
    Year                 = "Year",
    Dataset_or_consortium = "Dataset / Consortium",
    Phenotype_group      = "Phenotype group",
    Ancestry             = "Ancestry",
    Sample_size          = "Sample size",
    Sex                  = "Sex",
    Main_adjustments     = "Main covariate adjustments",
    Notes_or_accession   = "Notes / accession"
  ) %>%
  gt::tab_options(
    table.font.size = 11,
    heading.align   = "left"
  )

table1_gwas_gt

###############################################################################
# 4) MAIN TABLE 2 — SUMMARY OF 30 OUTCOMES & SAMPLE SIZES
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
  # FinnGen R12 — sample sizes from FinnGen R12 phenotype manifest
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
# 5) MAIN TABLE 3 — PRIMARY IVW MR ESTIMATES (FDR-CORRECTED)
#    Ordered by domain (same order as Figure 2)
###############################################################################

log_info("Creating Table 3: primary IVW MR estimates with FDR correction...")

# Domain ordering — mirrors outcome_meta in script 07
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
    Effect_scale = if_else(Type == "Binary", "OR", "β"),
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
      sprintf("(%.2f–%.2f)", CI_low_val, CI_high_val),
      sprintf("(%.3f–%.3f)", CI_low_val, CI_high_val)
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
# 6) SUPPLEMENTARY TABLES 1A–C — DEFINITIONS OF PERINATAL OUTCOMES
###############################################################################

log_info("Creating Supplementary Tables 1A–C (definitions of perinatal outcomes)...")

Supp_1A <- tibble::tribble(
  ~`Binary outcomes`,                ~`Case definition`,                                                ~`Control definition`,                                        ~`Exclusion criteria`,                                                                                         ~`Contributing studies`,
  "Pregnancy loss outcomes",         NA_character_,                                                    NA_character_,                                                NA_character_,                                                                                                 NA_character_,
  "Miscarriage",                     "≥1 pregnancy loss before 20 gestational weeks",                  "No pregnancy loss",                                          "Multiple births",                                                                                             "ALSPAC, MoBa, UKB, FinnGen",
  "Stillbirth",                      "≥1 pregnancy loss at or after 20 gestational weeks",             "No pregnancy loss",                                          "Multiple births",                                                                                             "ALSPAC, BiB, MoBa, UKB",
  "Maternal morbidity outcomes",     NA_character_,                                                    NA_character_,                                                NA_character_,                                                                                                 NA_character_,
  "GDM",                             "Diabetes mellitus diagnosed in pregnancy",                       "No GDM or pre-existing diabetes mellitus",                   "Pre-existing diabetes, multiple births, non-live births",                                                     "BiB, MoBa, UKB, GenDIP",
  "Perinatal depression",            "Maternal depression during pregnancy and up to one year after birth", "No maternal depression",                               "Pre-existing depression, multiple births, non-live births",                                                   "ALSPAC, MoBa, UKB, PGC",
  "Labour outcomes",                 NA_character_,                                                    NA_character_,                                                NA_character_,                                                                                                 NA_character_,
  "Induction of labour",             "Artificial stimulation of uterine contractions",                 "No artificial stimulation of uterine contractions",         "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "Prelabour rupture of membranes",  "Rupture of the amniotic sac before 37 gestational weeks",       "No rupture of the amniotic sac before 37 gestational weeks","Multiple births, non-live births",                                                                           "ALSPAC, MoBa, FinnGen",
  "Caesarean section",               "Delivery by caesarean section",                                 "No delivery by caesarean section",                          "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB, FinnGen",
  "Offspring birth outcomes",        NA_character_,                                                    NA_character_,                                                NA_character_,                                                                                                 NA_character_,
  "LBW",                             "<2,500 g",                                                       "≥2,500 to 4,500 g",                                          "Multiple births, non-live births, PTBs (GA <37 weeks)",                                                      "ALSPAC, MoBa, UKB",
  "HBW",                             ">4,500 g",                                                       "≥2,500 to 4,500 g",                                          "Multiple births, non-live births, PTBs (GA <37 weeks)",                                                      "ALSPAC, MoBa, UKB",
  "Pre-term birth",                  "GA <37 weeks",                                                   "GA ≥37 to <42 weeks",                                       "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB, FinnGen, EGG",
  "Post-term birth",                 "GA ≥42 weeks",                                                   "GA ≥37 to <42 weeks",                                       "Multiple births, non-live births, elective caesarean section, physician-induced labour",                      "ALSPAC, BiB, MoBa, UKB, FinnGen, EGG",
  "SGA",                             "Birth weight <10th percentile for GA",                          "Birth weight ≥10th percentile for GA",                      "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "LGA",                             "Birth weight >90th percentile for GA",                          "Birth weight ≤90th percentile for GA",                      "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "Low Apgar score at 1 minute",     "<7 points",                                                      "≥7 points",                                                  "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa",
  "Low Apgar score at 5 minutes",    "<7 points",                                                      "≥7 points",                                                  "Multiple births, non-live births",                                                                           "ALSPAC, MoBa",
  "NICU admission",                  "Offspring admitted to the NICU",                                "Offspring not admitted to the NICU",                        "Multiple births, non-live births",                                                                           "ALSPAC, MoBa"
)

s1A_file <- file.path(tables_dir, "Supp_Table_1A_primary_perinatal_definitions.csv")
write.csv(Supp_1A, s1A_file, row.names = FALSE)
log_info("Supplementary Table 1A saved: ", s1A_file)

Supp_1B <- tibble::tribble(
  ~`Binary outcomes`,          ~`Case definition`,                                      ~`Control definition`,                   ~`Exclusion criteria`,                                                                                         ~`Contributing studies`,
  "Pregnancy loss outcomes",   NA_character_,                                          NA_character_,                          NA_character_,                                                                                                 NA_character_,
  "Sporadic miscarriage",      "1–2 pregnancy loss before 20 gestational weeks",       "No pregnancy loss",                    "Age at menarche <9 or >17, underlying conditions leading to miscarriage (see supplement), multiple births",   "ALSPAC, MoBa, UKB",
  "Recurrent miscarriage",     "≥3 pregnancy loss before 20 gestational weeks",        "No pregnancy loss",                    "Age at menarche <9 or >17, underlying conditions leading to miscarriage (see supplement), multiple births",   "ALSPAC, MoBa, UKB, FinnGen",
  "Labour outcomes",           NA_character_,                                          NA_character_,                          NA_character_,                                                                                                 NA_character_,
  "Emergency caesarean section","Delivery by emergency caesarean section",            "No delivery by caesarean section",     "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "Elective caesarean section","Delivery by elective caesarean section",               "No delivery by caesarean section",     "Multiple births, non-live births",                                                                           "ALSPAC, BiB, MoBa, UKB",
  "Offspring birth outcomes",  NA_character_,                                          NA_character_,                          NA_character_,                                                                                                 NA_character_,
  "Very PTB",                  "GA <34 weeks",                                         "GA ≥37 to <42 weeks",                  "Multiple births, non-live births, elective caesarean section, physician-induced labour",                      "ALSPAC, BiB, MoBa, UKB",
  "Spontaneous PTB",           "GA <37 weeks",                                         "GA ≥37 to <42 weeks",                  "Multiple births, non-live births, elective caesarean section, physician-induced labour",                      "ALSPAC, BiB, MoBa, UKB"
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
# 7) TABLE S1 — HARMONISED DATASET SUMMARY (SNPs PER OUTCOME)
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
# 8) SUPPLEMENTARY TABLE 2 — GENETIC INSTRUMENTS FOR ENDOMETRIOSIS
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
# 9) SUPPLEMENTARY TABLE S3 — SENSITIVITY ANALYSES (WIDE FORMAT, eBioMedicine)
#
#   One row per outcome, ordered by the 9 clinical domains (same as Figure 2 /
#   Table 3). Columns: Domain, Scale (OR / β), IVW, MR-Egger, Weighted Median,
#   Weighted Mode effect estimates, Egger intercept, Cochran's Q p-value, and
#   MR-PRESSO as a diagnostic block.
#   Effect estimates are ORs for binary outcomes and β for continuous outcomes;
#   all column headers use "effect (95% CI)" to avoid misleading OR labelling
#   for continuous outcomes.
###############################################################################

log_info("Creating Table S3: wide-format sensitivity analyses (eBioMedicine)...")

# ── Helper functions ──────────────────────────────────────────────────────────

# Scalar: returns "X.XX (X.XX–X.XX)" for OR (binary) or "X.XXX (X.XXX–X.XXX)" for β
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
    # "effect (95% CI)" headers — generic label covers both OR and β
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
# Format: "β = X.XXX (p = X.XXX)" — explicit β symbol clarifies it is always
# a regression intercept on the log-OR / β scale, regardless of outcome type.
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
  log_warn("mr_pleiotropy.csv not found — Egger intercept column will be NA.")
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
  log_warn("mr_heterogeneity.csv not found — computed on the fly.")
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
  log_warn("mr_presso_results.csv not found — PRESSO columns will be NA.")
}

# ── Join all components + domain ordering ─────────────────────────────────────
Table_S3_wide <- s3_wide %>%
  dplyr::left_join(pleio_raw,     by = "outcome") %>%
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
    Scale,                                  # OR or β — clarifies metric per row
    `IVW effect (95% CI)`,
    `IVW p-value`,
    `MR-Egger effect (95% CI)`,
    `Egger intercept (beta), p-value`,
    `Weighted median effect (95% CI)`,
    `Weighted mode effect (95% CI)`,
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
  #                Egger intercept, WM eff, WMode eff, Q p, PRESSO n,
  #                PRESSO corr, PRESSO distort
  openxlsx::setColWidths(
    wb_s3, "S3 Sensitivity analyses",
    cols   = seq_len(ncol(tbl_s3_export)),
    widths = c(26, 32, 7, 22, 12, 22, 30, 22, 22, 14, 16, 26, 22)
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
# 10) SUPPLEMENTARY TABLE S4 — TRIO-BASED MR ESTIMATES
#
#   Maternal, fetal, and paternal genetic effects of endometriosis liability
#   on pregnancy outcomes (DONUTS mutually adjusted estimates from MR-PREG).
#   Source: results/trios_adj_mr_results_by_outcome_long.csv (script 04.2).
#   Ordered by domain (same as Table 3 / Figure 2).
###############################################################################

log_info("Creating Supplementary Table S4: trio-based MR estimates...")

trios_long_path <- file.path(results_dir, "trios_adj_mr_results_by_outcome_long.csv")

if (!file.exists(trios_long_path)) {
  warning("trios_adj_mr_results_by_outcome_long.csv not found — skipping Table S4.")
} else {

  trios_long <- readr::read_csv(trios_long_path, show_col_types = FALSE)

  # Origin display labels (all four components for the supplementary table)
  origin_levels <- c(
    "Maternal (unadj.)",
    "Maternal (adj. fetal)",
    "Fetal (adj. maternal)",
    "Paternal (adj.)"
  )

  fmt_trios_est <- function(or, lcl, ucl, b, b_lci, b_uci, sc) {
    dplyr::if_else(
      sc == "binary",
      sprintf("%.2f (%.2f–%.2f)", or, lcl, ucl),
      sprintf("%.3f (%.3f–%.3f)", b, b_lci, b_uci)
    )
  }

  Table_S4_trios <- trios_long %>%
    dplyr::filter(
      origin %in% origin_levels,
      outcome_id %in% vars_keep
    ) %>%
    dplyr::left_join(
      domain_meta %>% dplyr::select(outcome = outcome, domain, order_within),
      by = c("outcome_id" = "outcome")
    ) %>%
    dplyr::mutate(
      Outcome  = dplyr::recode(outcome_id, !!!outcome_labels, .default = outcome_id),
      Domain   = factor(domain, levels = domain_levels_t3),
      Origin   = factor(origin, levels = origin_levels),
      Scale    = dplyr::if_else(scale == "binary", "OR", "β"),
      b_est    = log(OR),
      b_lci_v  = log(LCL),
      b_uci_v  = log(UCL),
      `Effect (95% CI)` = dplyr::if_else(
        scale == "binary",
        sprintf("%.2f (%.2f–%.2f)", OR, LCL, UCL),
        sprintf("%.3f (%.3f–%.3f)", OR, LCL, UCL)
      ),
      `P value` = dplyr::case_when(
        pval < 0.001 ~ formatC(pval, format = "e", digits = 1),
        TRUE         ~ sprintf("%.3f", pval)
      ),
      `q (FDR)` = dplyr::case_when(
        !is.na(qval) & qval < 0.001 ~ formatC(qval, format = "e", digits = 1),
        !is.na(qval)                ~ sprintf("%.3f", qval),
        TRUE                        ~ NA_character_
      )
    ) %>%
    dplyr::arrange(Domain, order_within, Origin) %>%
    dplyr::select(
      Domain,
      Outcome,
      `Genetic effect` = Origin,
      Scale,
      `Effect (95% CI)`,
      `P value`,
      `q (FDR)`
    )

  # CSV
  readr::write_csv(
    Table_S4_trios,
    file.path(tables_dir, "Supplementary_Table_S4_trio_based_MR.csv")
  )
  log_info("Supplementary Table S4 (trio-based MR, CSV) saved.")

  # Excel
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_s4 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s4, "S4 Trio-based MR")

    caption_s4 <- paste0(
      "Supplementary Table S4. Trio-based Mendelian randomization estimates of ",
      "maternal, fetal, and paternal genetic effects of endometriosis liability on ",
      "pregnancy outcomes. Estimates are derived from mutually adjusted (DONUTS) ",
      "models using the MR-PREG trio-based GWAS. OR = odds ratio for binary outcomes; ",
      "β = beta coefficient for continuous outcomes (gestational age, birthweight z-score). ",
      "FDR q-values are provided where available."
    )
    openxlsx::writeData(wb_s4, "S4 Trio-based MR",
                        x = caption_s4, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb_s4, "S4 Trio-based MR",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10,
                         textDecoration = "bold", wrapText = TRUE
                       ),
                       rows = 1, cols = 1)

    openxlsx::writeDataTable(
      wb_s4, "S4 Trio-based MR",
      x          = Table_S4_trios,
      startRow   = 3, startCol = 1,
      tableStyle = "TableStyleLight9",
      withFilter = TRUE
    )

    openxlsx::setColWidths(wb_s4, "S4 Trio-based MR",
                           cols   = seq_len(ncol(Table_S4_trios)),
                           widths = c(28, 32, 24, 8, 22, 12, 12))

    openxlsx::addStyle(wb_s4, "S4 Trio-based MR",
                       style = openxlsx::createStyle(
                         fontName = "Arial", fontSize = 10, textDecoration = "bold",
                         fgFill = "#D9E1F2", border = "Bottom", borderColour = "#4472C4",
                         wrapText = TRUE
                       ),
                       rows = 3, cols = seq_len(ncol(Table_S4_trios)),
                       gridExpand = TRUE)

    openxlsx::saveWorkbook(
      wb_s4,
      file.path(tables_dir, "Supplementary_Table_S4_trio_based_MR.xlsx"),
      overwrite = TRUE
    )
    log_info("Supplementary Table S4 (trio-based MR, Excel) saved.")
  }

  message("Supplementary Table S4 (trio-based MR) done.")
}

log_info("=== Table generation for endoMR-PREG completed successfully ===")

###############################################################################
# 12) EXPORT ALL TABLES TO HTML
###############################################################################

log_info("Exporting all tables to HTML...")

suppressPackageStartupMessages(library(gt))

html_dir <- file.path(tables_dir, "html")
dir.create(html_dir, showWarnings = FALSE, recursive = TRUE)

tables_to_export <- list(
  Table1_GWAS_sources                      = table1_gwas,
  Table2_pregnancy_sample_sizes            = table1_preg,
  Table3_IVW_FDR_main_results              = table3_main,
  Supp_Table_S1_phenotype_definitions_A    = Supp_1A,
  Supp_Table_S1_phenotype_definitions_B    = Supp_1B,
  Supp_Table_S1_phenotype_definitions_C    = Supp_1C,
  Supp_Table_S2_instruments                = supp_table2,
  Supp_Table_S3_sensitivity_analyses       = Table_S3_csv,
  Supp_Table_S4_trio_based_MR             = if (exists("Table_S4_trios")) Table_S4_trios else NULL
)

table_titles <- list(
  Table1_GWAS_sources =
    "Table 1. Overview of GWAS datasets for endometriosis and pregnancy outcomes in the endoMR-PREG study",

  Table2_pregnancy_sample_sizes =
    "Table 2. Outcome definitions, sample sizes, case/control counts, and SNP availability per outcome",

  Table3_IVW_FDR_main_results =
    "Table 3. Primary IVW MR results across 30 outcomes (with FDR q-values)",

  Supp_Table_S1_phenotype_definitions_A =
    "Supplementary Table S1A. Detailed phenotype definitions and contributing cohorts — primary binary outcomes",

  Supp_Table_S1_phenotype_definitions_B =
    "Supplementary Table S1B. Detailed phenotype definitions and contributing cohorts — secondary binary outcomes",

  Supp_Table_S1_phenotype_definitions_C =
    "Supplementary Table S1C. Detailed phenotype definitions and contributing cohorts — continuous outcomes",

  Supp_Table_S2_instruments =
    "Supplementary Table S2. Endometriosis instrument SNP list and characteristics (post-clumping)",

  Supp_Table_S3_sensitivity_analyses =
    "Supplementary Table S3. Sensitivity analyses including MR-Egger, weighted median, and weighted mode estimates, heterogeneity (Cochran's Q), pleiotropy (Egger intercept), and MR-PRESSO results",

  Supp_Table_S4_trio_based_MR =
    "Supplementary Table S4. Trio-based MR estimates by maternal, fetal, and paternal genetic effects"
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
# SUPPLEMENTARY TABLE S5 — LEAVE-ONE-COHORT SENSITIVITY ANALYSIS
#
#   IVW estimates across all outcomes with each contributing study
#   sequentially excluded (leave-one-cohort-out design).
#   Source: results/leave_one_cohort_results.csv (produced by script 05).
#   Outcomes and domains ordered as in Table 3 / Figure 2.
###############################################################################

loc_path <- file.path(results_dir, "leave_one_cohort_results.csv")

if (!file.exists(loc_path)) {
  warning("leave_one_cohort_results.csv not found — skipping Table S5.")
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
