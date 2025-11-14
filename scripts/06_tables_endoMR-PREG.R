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
# 1) LOAD HARMONISED DATA & DEFINE 29 PRIMARY OUTCOMES
###############################################################################

dat <- data.table::fread(harm_file)
log_info("Loaded harmonised data with ", nrow(dat), " rows")

stopifnot(all(c("id.exposure", "beta.exposure", "beta.outcome") %in% colnames(dat)))
stopifnot("outcome" %in% colnames(dat))

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
  
  # Growth / GA / weight (7)
  "ga_all",
  "ga_subsamp",
  "sga",
  "lbw_all",
  "hbw_all",
  "lga",
  "zbw_all",
  
  # Neonatal / Apgar (3)
  "lowapgar1",
  "lowapgar5",
  "nicu",
  
  # Maternal complications (5)
  "anaemia_preg_all",
  "gdm_subsamp",
  "gh_subsamp",
  "hdp_subsamp",
  "pe_subsamp",
  
  # Other obstetric ≥20 SA (2)
  "induction",
  "posttb_all"
)

outcome_labels <- c(
  Antepartum_bleeding                       = "Antepartum bleeding",
  Postpartum_hemorrhage                     = "Postpartum hemorrhage",
  Postpartum_hemorrhage_due_to_atony        = "PPH due to atony",
  Postpartum_hemorrhage_due_to_retained_placenta = "PPH due to retained placenta",
  finngen_R12_O15_PLAC_PRAEVIA              = "Placenta praevia",
  finngen_R12_O15_PLAC_DISORD               = "Placental disorders",
  finngen_R12_O15_PLAC_PREMAT_SEPAR         = "Premature placental separation",
  rup_memb                                  = "Premature rupture of membranes",
  pretb_all                                 = "Preterm birth (all)",
  vpretb_all                                = "Very preterm birth",
  ga_all                                    = "Gestational age (all)",
  ga_subsamp                                = "Gestational age (subsample)",
  sga                                       = "Small for gestational age",
  lbw_all                                   = "Low birthweight",
  hbw_all                                   = "High birthweight",
  lga                                       = "Large for gestational age",
  zbw_all                                   = "Z-score birthweight",
  lowapgar1                                 = "Low Apgar score at 1 min",
  lowapgar5                                 = "Low Apgar score at 5 min",
  nicu                                      = "NICU admission",
  anaemia_preg_all                          = "Pregnancy anemia",
  gdm_subsamp                               = "Gestational diabetes",
  gh_subsamp                                = "Gestational hypertension",
  hdp_subsamp                               = "Hypertensive disorders of pregnancy",
  pe_subsamp                                = "Preeclampsia",
  induction                                 = "Labour induction",
  posttb_all                                = "Post-term birth"
)

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

continuous_outcomes <- c("zbw_all")

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

snps_file <- file.path(results_dir, "endometriosis_clumped_snps.tsv")
if (file.exists(snps_file)) {
  clumped_snps  <- data.table::fread(snps_file)
  snp_whitelist <- unique(clumped_snps$SNP)
} else {
  snp_whitelist <- unique(dat$SNP)
  log_warn("Clumped SNP file not found. Using SNPs from harmonised data (n = ",
           length(snp_whitelist), ").")
}
log_info("Using ", length(snp_whitelist), " SNP instruments")

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
  "McBride et al. (MR-PREG)",  2025,  "MR-PREG collaboration (ALSPAC, BiB, MoBa, UKB, FinnGen + GWAS)",  "Adverse pregnancy and perinatal outcomes (binary and continuous)",   "Predominantly European",                              "Up to 678,001 women (outcome-specific)",               "Female",       "Age, study/centre, principal components",              "Core source for hypertensive disorders, GDM, PTB, SGA/LGA, CS, PROM, stillbirth, Apgar, NICU",
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
# 4) MAIN TABLE 2 — SUMMARY OF 29 OUTCOMES & SAMPLE SIZES
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

manual_sizes <- tibble::tribble(
  ~Source,             ~Outcome,                             ~`N total`, ~Cases,  ~Controls,
  "Westergaard (PPH)", "Antepartum bleeding",                  331792L,     3236L,   328556L,
  "Westergaard (PPH)", "Early bleeding (live birth)",          331792L,     6356L,   325436L,
  "Westergaard (PPH)", "Early bleeding (any outcome)",         331792L,    28898L,   302894L,
  "Westergaard (PPH)", "PPH due to atony",                     274857L,    13048L,   261809L,
  "Westergaard (PPH)", "PPH due to retained placenta",         272683L,     6256L,   266427L,
  "Westergaard (PPH)", "Postpartum hemorrhage",                331792L,    21521L,   310271L
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
###############################################################################

log_info("Creating Table 3: primary IVW MR estimates with FDR correction...")

ivw_res <- load_mr_results("ivw_results") %>%
  dplyr::filter(
    outcome %in% vars_keep,
    method == "Inverse variance weighted"
  ) %>%
  dplyr::left_join(type_df, by = "outcome") %>%
  dplyr::mutate(
    Outcome_label = dplyr::recode(outcome, !!!outcome_labels, .default = outcome),
    `Data source` = label_source(outcome),
    q_fdr         = p.adjust(pval, method = "fdr")
  )

table3_main <- ivw_res %>%
  dplyr::mutate(
    Effect_scale = if_else(Type == "Binary", "OR", "Beta"),
    Estimate_val = if_else(Type == "Binary", exp(b), b),
    CI_low_val   = if_else(
      Type == "Binary",
      exp(b - 1.96 * se),
      b - 1.96 * se
    ),
    CI_high_val  = if_else(
      Type == "Binary",
      exp(b + 1.96 * se),
      b + 1.96 * se
    ),
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
  dplyr::select(
    `Data source`,
    Outcome      = outcome,
    Outcome_label,
    Type,
    Effect_scale,
    `Estimate`,
    `95% CI`,
    `P-value`,
    `q (FDR)`,
    `No. SNPs` = nsnp
  ) %>%
  dplyr::arrange(`Data source`, Outcome_label)

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
  dplyr::select(
    `Data source`,
    Outcome       = outcome,
    Outcome_label,
    Type,
    `No. SNPs`
  ) %>%
  dplyr::arrange(`Data source`, Outcome_label)

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

if (!exists("clumped2")) {
  stop("Object 'clumped2' is not available. Please load the clumped instrument dataset (clumped2) before running this script.")
}

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
# 9) TABLE S3 — FULL MR RESULTS FOR ALL METHODS
###############################################################################

log_info("Creating Table S3: all MR methods results...")

all_res <- load_mr_results("all_mr_methods")

label_group <- function(outcome, outcome_label) {
  placental_disorders <- c(
    "finngen_R12_O15_PLAC_PRAEVIA",
    "finngen_R12_O15_PLAC_DISORD",
    "finngen_R12_O15_PLAC_PREMAT_SEPAR"
  )
  
  hdp_outcomes <- c(
    "hdp_subsamp",
    "gh_subsamp",
    "pe_subsamp"
  )
  
  pregnancy_timing <- c(
    "ga_all",
    "ga_subsamp",
    "pretb_all",
    "vpretb_all",
    "posttb_all"
  )
  
  fetal_growth_bw <- c(
    "hbw_all", "lbw_all", "sga",
    "lga", "lga_all",
    "zbw_all",
    "bw_z", "bw_zscore", "bw_z_score",
    "birthweight_z", "birth_weight_z"
  )
  
  labour_delivery <- c(
    "induction",
    "rup_memb",
    "cs",
    "el_cs",
    "em_cs"
  )
  
  bleeding_haem <- c(
    "Antepartum_bleeding",
    "Early_bleeding_with_any_outcome",
    "Early_bleeding_ending_in_live_birth",
    "Postpartum_hemorrhage",
    "Postpartum_hemorrhage_due_to_atony",
    "Postpartum_hemorrhage_due_to_retained_placenta",
    "antepartum_bleeding",
    "ap_bleeding",
    "pph_all",
    "pph_atony",
    "pph_retained_placenta",
    "pph_other"
  )
  
  maternal_metab_haem <- c(
    "gdm_subsamp",
    "anaemia_preg_all"
  )
  
  maternal_mental <- c(
    "depr_subsamp"
  )
  
  neonatal_condition <- c(
    "lowapgar1",
    "lowapgar5",
    "nicu",
    "sb_subsamp"
  )
  
  dplyr::case_when(
    outcome %in% placental_disorders ~ "Placental disorders",
    outcome %in% hdp_outcomes ~ "Hypertensive disorders of pregnancy",
    outcome %in% pregnancy_timing |
      grepl("^ga", outcome, ignore.case = TRUE) ~ "Pregnancy timing",
    outcome %in% fetal_growth_bw |
      outcome_label %in% c("Z-score birthweight") ~ "Fetal growth and birthweight",
    outcome %in% labour_delivery ~ "Labour and delivery complications",
    outcome %in% bleeding_haem |
      outcome_label %in% c(
        "Antepartum bleeding",
        "PPH due to atony",
        "PPH due to retained placenta",
        "Postpartum hemorrhage"
      ) ~ "Bleeding and haemorrhage",
    outcome %in% maternal_metab_haem ~ "Maternal metabolic/haematologic complications",
    outcome %in% maternal_mental |
      outcome_label %in% c("Postpartum depression", "Perinatal depression") ~
      "Maternal mental health",
    outcome %in% neonatal_condition ~ "Neonatal condition at birth",
    TRUE ~ "Other"
  )
}

domain_outcomes <- c(
  "finngen_R12_O15_PLAC_PRAEVIA",
  "finngen_R12_O15_PLAC_DISORD",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",
  "hdp_subsamp", "gh_subsamp", "pe_subsamp",
  "ga_all", "ga_subsamp", "pretb_all", "vpretb_all", "posttb_all",
  "hbw_all", "lbw_all", "sga",
  "lga", "lga_all",
  "zbw_all",
  "bw_z", "bw_zscore", "bw_z_score",
  "birthweight_z", "birth_weight_z",
  "induction", "rup_memb", "cs", "el_cs", "em_cs",
  "Antepartum_bleeding", "Early_bleeding_with_any_outcome", "Early_bleeding_ending_in_live_birth",
  "Postpartum_hemorrhage", "Postpartum_hemorrhage_due_to_atony", "Postpartum_hemorrhage_due_to_retained_placenta",
  "antepartum_bleeding", "ap_bleeding",
  "pph_all", "pph_atony", "pph_retained_placenta", "pph_other",
  "gdm_subsamp", "anaemia_preg_all",
  "depr_subsamp",
  "lowapgar1", "lowapgar5", "nicu", "sb_subsamp"
)

vars_domains <- union(vars_keep, domain_outcomes)

has_exposure_source <- "exposure_source" %in% names(all_res)

Table_S3 <- all_res %>%
  dplyr::filter(outcome %in% vars_domains) %>%
  dplyr::mutate(
    Outcome    = dplyr::recode(outcome, !!!outcome_labels, .default = outcome),
    group      = label_group(outcome, Outcome),
    instrument = if (has_exposure_source) exposure_source else exposure,
    exposure   = "Genetic liability to endometriosis",
    type       = dplyr::if_else(outcome %in% continuous_outcomes, "continuous", "binary"),
    priority   = "primary"
  ) %>%
  dplyr::select(
    group,
    instrument,
    exposure,
    Outcome,
    method,
    nsnp,
    b,
    se,
    pval,
    type,
    priority
  ) %>%
  dplyr::arrange(group, Outcome, method)

s3_file <- file.path(tables_dir, "Table_S3_all_MR_methods.csv")
write.csv(Table_S3, s3_file, row.names = FALSE)
log_info("Table S3 (all MR methods) saved: ", s3_file)

###############################################################################
# 10) TABLE S4 — HETEROGENEITY STATISTICS (COCHRAN'S Q, IVW ONLY)
###############################################################################

log_info("Creating Table S4: heterogeneity statistics (IVW only)...")

dat_het <- data.table::fread(harm_file)
dat_het <- dat_het[dat_het$outcome %in% vars_keep, , drop = FALSE]

het_res <- TwoSampleMR::mr_heterogeneity(dat_het)

Table_S4 <- het_res %>%
  dplyr::filter(
    outcome %in% vars_keep,
    method %in% c("Inverse variance weighted", "IVW")
  ) %>%
  dplyr::mutate(
    Outcome  = dplyr::recode(outcome, !!!outcome_labels, .default = outcome),
    group    = label_group(outcome, Outcome),
    exposure = "Genetic liability to endometriosis",
    type     = dplyr::if_else(outcome %in% continuous_outcomes,
                              "continuous", "binary"),
    priority = "primary"
  ) %>%
  dplyr::select(
    group,
    exposure,
    Outcome,
    Q,
    Q_df,
    Q_pval,
    type,
    priority
  ) %>%
  dplyr::arrange(group, Outcome)

s4_file <- file.path(tables_dir, "Table_S4_heterogeneity_IVW.csv")
write.csv(Table_S4, s4_file, row.names = FALSE)
log_info("Table S4 (heterogeneity, IVW) saved: ", s4_file)

###############################################################################
# 11) TABLE S5 — MR-EGGER INTERCEPT (HORIZONTAL PLEIOTROPY)
###############################################################################

log_info("Creating Table S5: MR-Egger intercept (pleiotropy)...")

pleio_dat <- data.table::fread(harm_file)
pleio_dat <- pleio_dat[pleio_dat$outcome %in% vars_keep, , drop = FALSE]

egger_res <- TwoSampleMR::mr_pleiotropy_test(pleio_dat)

Table_S5 <- egger_res %>%
  dplyr::filter(outcome %in% vars_keep) %>%
  dplyr::mutate(
    outcome_code = outcome,
    Outcome      = dplyr::recode(outcome_code, !!!outcome_labels, .default = outcome_code),
    group        = label_group(outcome_code, Outcome),
    exposure     = "Genetic liability to endometriosis",
    type         = dplyr::if_else(outcome_code %in% continuous_outcomes,
                                  "continuous", "binary"),
    priority     = "primary"
  ) %>%
  dplyr::select(
    group,
    exposure,
    Outcome,
    egger_intercept,
    se,
    pval,
    type,
    priority
  ) %>%
  dplyr::arrange(group, Outcome)

s5_file <- file.path(tables_dir, "Table_S5_egger_intercept.csv")
write.csv(Table_S5, s5_file, row.names = FALSE)
log_info("Table S5 (MR-Egger intercept) saved: ", s5_file)

log_info("=== Table generation for endoMR-PREG completed successfully ===")
