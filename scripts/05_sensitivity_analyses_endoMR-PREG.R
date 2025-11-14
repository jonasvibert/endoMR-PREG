#!/usr/bin/env Rscript
###############################################################################
# Script: 05_sensitivity_analyses_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   1) Run standard TwoSampleMR sensitivity analyses from the harmonised dataset:
#        - Cochran’s Q (mr_heterogeneity)
#        - MR-Egger intercept (mr_pleiotropy_test)
#        - Single-SNP MR (mr_singlesnp)
#        - Leave-one-out by SNP (mr_leaveoneout)
#   2) Run leave-one-cohort-out MR analyses using per-study MR-PREG outcome GWAS.
#
# Inputs:
#   - results/harmonised_rahmioglu_bpo.csv
#   - data/OUTCOME_MR-PREG/stu_out_dat.txt
#
# Outputs (saved under results/):
#   - mr_heterogeneity_*.csv
#   - mr_pleiotropy_*.csv
#   - mr_single_snp_*.csv
#   - mr_leaveoneout_snp_*.csv
#   - leave_one_cohort_results.csv
#   - cohort_composition.csv
###############################################################################

### 1) SETUP ####################################################################

required_pkgs <- c(
  "TwoSampleMR", "dplyr", "ggplot2", "here", "readr",
  "data.table", "stringr", "tidyr", "purrr"
)

safe_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

invisible(lapply(required_pkgs, safe_install))

results_dir  <- here::here("results")
plots_dir    <- here::here("plots")
data_dir     <- here::here("data")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)

harm_file    <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")
outcome_file <- file.path(data_dir, "OUTCOME_MR-PREG", "stu_out_dat.txt")

if (!file.exists(harm_file)) {
  stop("Harmonised file not found: ", harm_file)
}
if (!file.exists(outcome_file)) {
  stop("Outcome file not found: ", outcome_file)
}

dat       <- data.table::fread(harm_file) |> as.data.frame()
harm_data <- dat  # reuse for exposure reconstruction

stopifnot(all(c("SNP", "beta.exposure", "se.exposure",
                "effect_allele.exposure", "other_allele.exposure",
                "eaf.exposure", "id.exposure", "exposure") %in% names(harm_data)))
stopifnot(all(c("beta.outcome", "se.outcome", "pval.outcome",
                "effect_allele.outcome", "other_allele.outcome",
                "eaf.outcome", "id.outcome", "outcome") %in% names(dat)))

message("Loaded harmonised dataset: ", nrow(dat), " rows.")


### 2) HELPERS ##################################################################

if (!exists("export_csv")) {
  export_csv <- function(df, stem, dir_out = NULL) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
      message("Nothing to export for ", stem)
      return(invisible(NULL))
    }
    
    if (is.null(dir_out)) {
      if (exists("results_dir", inherits = TRUE)) {
        dir_out <- get("results_dir", inherits = TRUE)
      } else {
        dir_out <- "."
      }
    }
    
    if (!dir.exists(dir_out)) {
      dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    }
    
    fn <- file.path(
      dir_out,
      paste0(stem, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
    )
    readr::write_csv(df, fn)
    message("Written: ", fn)
    invisible(fn)
  }
}

if (!exists("labels_df")) {
  labels_df <- data.frame(
    outcome      = character(),
    outcome_full = character(),
    stringsAsFactors = FALSE
  )
}


### 3) STANDARD SENSITIVITY ANALYSES ###########################################

message("=== Standard TwoSampleMR sensitivity analyses ===")

het_res     <- mr_heterogeneity(dat)      # Cochran’s Q
plt_res     <- mr_pleiotropy_test(dat)    # MR-Egger intercept
single_res  <- mr_singlesnp(dat)          # Single-SNP MR
loo_snp_res <- mr_leaveoneout(dat)        # Leave-one-out by SNP

export_csv(het_res,     "mr_heterogeneity")
export_csv(plt_res,     "mr_pleiotropy")
export_csv(single_res,  "mr_single_snp")
export_csv(loo_snp_res, "mr_leaveoneout_snp")


### 4) LEAVE-ONE-COHORT-OUT MR (MR-PREG) #######################################

message("=== Leave-one-cohort-out MR analysis ===")
set.seed(1)

# 4.1 Load raw outcome data ----------------------------------------------------

mr_data <- readr::read_table(outcome_file, col_types = cols()) |> as.data.frame()

required_cols_out <- c(
  "SNP", "beta", "se", "pval", "samplesize", "study",
  "Phenotype", "effect_allele", "other_allele", "eaf",
  "ncase", "ncontrol"
)
missing_out <- setdiff(required_cols_out, names(mr_data))
if (length(missing_out)) {
  stop("Missing columns in outcome data: ", paste(missing_out, collapse = ", "))
}

message("Outcome rows: ", nrow(mr_data))

# 4.2 Cohort composition -------------------------------------------------------

cohort_summary <- mr_data |>
  dplyr::group_by(study) |>
  dplyr::summarise(
    n_snps         = dplyr::n_distinct(SNP),
    n_phenotypes   = dplyr::n_distinct(Phenotype),
    mean_samplesize = mean(samplesize, na.rm = TRUE),
    .groups        = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_snps))

message("Available cohorts: ", paste(cohort_summary$study, collapse = ", "))
if (nrow(cohort_summary) < 2) {
  stop("Need ≥ 2 cohorts for leave-one-cohort-out analysis.")
}

# 4.3 Outcome labels and types (restricted to the 29 kept outcomes) -----------

outcome_labels <- c(
  # Placental outcomes
  finngen_R12_O15_PLAC_PRAEVIA = "Placenta praevia",
  finngen_R12_O15_PLAC_DISORD  = "Placental disorders",
  finngen_R12_O15_PLAC_PREMAT_SEPAR = "Premature placental separation",
  
  # Caesarean delivery
  el_cs   = "Elective caesarean section",
  em_cs   = "Emergency caesarean section",
  cs      = "Caesarean section",
  
  # Apgar scores
  lowapgar1 = "Low Apgar score at 1 min",
  lowapgar5 = "Low Apgar score at 5 min",
  
  # Labour and delivery complications
  rup_memb  = "Premature rupture of membranes",
  induction = "Labour induction",
  
  # Gestational age and timing
  ga_all    = "Gestational age",
  pretb_all = "Preterm birth (any)",
  vpretb_all = "Very preterm birth",
  posttb_all = "Post-term birth",
  
  # Birth weight outcomes
  hbw_all   = "High birthweight (>4000g)",
  lbw_all   = "Low birthweight (<2500g)",
  sga       = "Small for gestational age",
  
  # Maternal health
  depr_subsamp     = "Postpartum depression",
  anaemia_preg_all = "Pregnancy anaemia",
  
  # Pregnancy complications
  gdm_subsamp = "Gestational diabetes mellitus",
  hdp_subsamp = "Hypertensive disorders of pregnancy",
  gh_subsamp  = "Gestational hypertension",
  pe_subsamp  = "Preeclampsia",
  
  # Neonatal outcomes
  nicu       = "NICU admission",
  sb_subsamp = "Stillbirth",
  
  # Haemorrhage and bleeding
  Antepartum_bleeding_filtered                        = "Antepartum bleeding",
  Postpartum_hemorrhage_filtered                      = "Postpartum haemorrhage",
  Postpartum_hemorrhage_due_to_atony_filtered         = "PPH due to atony",
  Postpartum_hemorrhage_due_to_retained_placenta_filtered = "PPH due to retained placenta"
)

outcome_types <- c(
  # Placental outcomes
  finngen_R12_O15_PLAC_PRAEVIA = "binary",
  finngen_R12_O15_PLAC_DISORD  = "binary",
  finngen_R12_O15_PLAC_PREMAT_SEPAR = "binary",
  
  # Caesarean delivery
  el_cs   = "binary",
  em_cs   = "binary",
  cs      = "binary",
  
  # Apgar scores
  lowapgar1 = "binary",
  lowapgar5 = "binary",
  
  # Labour and delivery complications
  rup_memb  = "binary",
  induction = "binary",
  
  # Gestational age and timing
  ga_all    = "continuous",
  pretb_all = "binary",
  vpretb_all = "binary",
  posttb_all = "binary",
  
  # Birth weight outcomes
  hbw_all   = "binary",
  lbw_all   = "binary",
  sga       = "binary",
  
  # Maternal health
  depr_subsamp     = "binary",
  anaemia_preg_all = "binary",
  
  # Pregnancy complications
  gdm_subsamp = "binary",
  hdp_subsamp = "binary",
  gh_subsamp  = "binary",
  pe_subsamp  = "binary",
  
  # Neonatal outcomes
  nicu       = "binary",
  sb_subsamp = "binary",
  
  # Haemorrhage and bleeding
  Antepartum_bleeding_filtered                        = "binary",
  Postpartum_hemorrhage_filtered                      = "binary",
  Postpartum_hemorrhage_due_to_atony_filtered         = "binary",
  Postpartum_hemorrhage_due_to_retained_placenta_filtered = "binary"
)

mr_data <- mr_data |>
  dplyr::mutate(
    id.outcome    = Phenotype,
    outcome_label = dplyr::coalesce(outcome_labels[Phenotype], Phenotype),
    outcome_type  = dplyr::coalesce(outcome_types[Phenotype], "binary")
  )


# 4.4 Exposure reconstruction from harmonised data -----------------------------

message("Preparing exposure dataset from harmonised file...")

exposure_df_raw <- harm_data |>
  dplyr::select(
    SNP, beta.exposure, se.exposure, pval.exposure,
    effect_allele.exposure, other_allele.exposure, eaf.exposure,
    samplesize.exposure, id.exposure, exposure
  ) |>
  dplyr::distinct(SNP, .keep_all = TRUE)

exposure_for_format <- exposure_df_raw |>
  dplyr::transmute(
    SNP,
    beta          = beta.exposure,
    se            = se.exposure,
    pval          = pval.exposure,
    effect_allele = effect_allele.exposure,
    other_allele  = other_allele.exposure,
    eaf           = eaf.exposure,
    samplesize    = samplesize.exposure,
    ncase         = NA_real_,
    ncontrol      = NA_real_,
    id            = id.exposure,
    exposure      = exposure,
    phenotype     = exposure
  )

exposure_formatted <- TwoSampleMR::format_data(
  exposure_for_format,
  type          = "exposure",
  phenotype_col = "phenotype"
)

message("Exposure SNPs (after deduplication): ",
        dplyr::n_distinct(exposure_formatted$SNP))


# 4.5 Helper functions: meta-analysis and formatting ---------------------------

meta_outcome_after_drop <- function(dat_one_outcome) {
  dat_one_outcome |>
    dplyr::filter(
      !is.na(SNP),
      is.finite(beta.outcome),
      is.finite(se.outcome)
    ) |>
    dplyr::group_by(SNP) |>
    dplyr::summarise(
      w    = sum(1 / (se.outcome^2), na.rm = TRUE),
      beta = ifelse(
        w > 0,
        sum(beta.outcome / (se.outcome^2), na.rm = TRUE) / w,
        NA_real_
      ),
      se   = ifelse(w > 0, sqrt(1 / w), NA_real_),
      pval = ifelse(
        w > 0,
        2 * pnorm(-abs(beta / se)),
        NA_real_
      ),
      effect_allele = {
        tmp <- stats::na.omit(effect_allele.outcome)
        if (length(tmp)) tmp[1] else NA_character_
      },
      other_allele  = {
        tmp <- stats::na.omit(other_allele.outcome)
        if (length(tmp)) tmp[1] else NA_character_
      },
      ss_sum     = sum(samplesize.outcome, na.rm = TRUE),
      eaf        = ifelse(
        ss_sum > 0 & any(!is.na(eaf.outcome)),
        sum(eaf.outcome * samplesize.outcome, na.rm = TRUE) / ss_sum,
        mean(eaf.outcome, na.rm = TRUE)
      ),
      samplesize = ss_sum,
      ncase      = sum(ncase.outcome, na.rm = TRUE),
      ncontrol   = sum(ncontrol.outcome, na.rm = TRUE),
      .groups    = "drop"
    ) |>
    dplyr::select(-w, -ss_sum) |>
    dplyr::filter(
      is.finite(beta),
      is.finite(se),
      !is.na(effect_allele),
      !is.na(other_allele)
    )
}

build_outcome_formatted <- function(meta_df, outcome_id, outcome_label) {
  if (is.null(meta_df) || nrow(meta_df) == 0) return(NULL)
  
  tmp <- meta_df |>
    dplyr::transmute(
      SNP, beta, se, pval,
      effect_allele, other_allele, eaf,
      samplesize, ncase, ncontrol,
      id        = outcome_id,
      outcome   = outcome_label,
      phenotype = outcome_label
    ) |>
    dplyr::filter(
      !is.na(SNP),
      is.finite(beta),
      is.finite(se),
      !is.na(effect_allele),
      !is.na(other_allele)
    )
  
  if (nrow(tmp) == 0) return(NULL)
  
  TwoSampleMR::format_data(tmp, type = "outcome", phenotype_col = "phenotype")
}


# 4.6 Run main + leave-one-cohort-out MR per outcome ---------------------------

message("Preparing outcome data frame for leave-one-cohort-out...")

mr_outcome_raw <- mr_data |>
  dplyr::transmute(
    SNP,
    beta.outcome         = beta,
    se.outcome           = se,
    pval.outcome         = pval,
    effect_allele.outcome = effect_allele,
    other_allele.outcome  = other_allele,
    eaf.outcome          = eaf,
    samplesize.outcome   = samplesize,
    ncase.outcome        = ncase,
    ncontrol.outcome     = ncontrol,
    id.outcome,
    outcome              = outcome_label,
    outcome_type,
    study
  )

outcomes_list <- split(mr_outcome_raw, mr_outcome_raw$id.outcome)

run_loo_for_outcome <- function(one_outcome_df) {
  outcome_id   <- unique(one_outcome_df$id.outcome)[1]
  outcome_name <- unique(one_outcome_df$outcome)[1]
  cohorts      <- sort(unique(one_outcome_df$study))
  
  message(">> Outcome: ", outcome_id, " | ", outcome_name,
          " | Cohorts: ", length(cohorts))
  
  # Main meta-analysis across all cohorts
  main_meta <- meta_outcome_after_drop(one_outcome_df)
  main_outcome_fmt <- build_outcome_formatted(main_meta, outcome_id, outcome_name)
  if (is.null(main_outcome_fmt) || nrow(main_outcome_fmt) < 3) return(NULL)
  
  main_harm <- tryCatch(
    TwoSampleMR::harmonise_data(exposure_formatted, main_outcome_fmt),
    error = function(e) NULL
  )
  if (is.null(main_harm) || dplyr::n_distinct(main_harm$SNP) < 3) return(NULL)
  
  main_ivw <- TwoSampleMR::mr(main_harm, method_list = "mr_ivw")
  if (!nrow(main_ivw)) return(NULL)
  
  main_ivw$analysis_type       <- "Main analysis"
  main_ivw$nsnp                <- dplyr::n_distinct(main_harm$SNP)
  main_ivw$left_out_cohort     <- "None"
  main_ivw$n_snps_remaining    <- main_ivw$nsnp
  main_ivw$n_cohorts_remaining <- length(cohorts)
  main_ivw$outcome             <- outcome_name
  main_ivw$id.outcome          <- outcome_id
  
  # Leave-one-cohort-out analyses
  loo_list <- lapply(cohorts, function(left_out) {
    sub_df <- dplyr::filter(one_outcome_df, study != left_out)
    
    meta   <- meta_outcome_after_drop(sub_df)
    out_fmt <- build_outcome_formatted(meta, outcome_id, outcome_name)
    if (is.null(out_fmt) || nrow(out_fmt) < 3) return(NULL)
    
    harm <- tryCatch(
      TwoSampleMR::harmonise_data(exposure_formatted, out_fmt),
      error = function(e) NULL
    )
    if (is.null(harm) || dplyr::n_distinct(harm$SNP) < 3) return(NULL)
    
    res <- TwoSampleMR::mr(harm, method_list = "mr_ivw")
    if (!nrow(res)) return(NULL)
    
    res$analysis_type       <- paste0("Without ", left_out)
    res$left_out_cohort     <- left_out
    res$n_snps_remaining    <- dplyr::n_distinct(harm$SNP)
    res$n_cohorts_remaining <- dplyr::n_distinct(sub_df$study)
    res$outcome             <- outcome_name
    res$id.outcome          <- outcome_id
    
    res
  })
  
  dplyr::bind_rows(main_ivw, dplyr::bind_rows(loo_list))
}

message("Running main + leave-one-cohort-out analyses...")
results_list    <- lapply(outcomes_list, run_loo_for_outcome)
all_results_raw <- dplyr::bind_rows(results_list)

if (!nrow(all_results_raw)) {
  stop("No MR results produced in leave-one-cohort-out analysis.")
}


# 4.7 Summary and export -------------------------------------------------------

all_results <- all_results_raw |>
  dplyr::mutate(
    outcome      = dplyr::coalesce(outcome, outcome_labels[id.outcome], id.outcome),
    outcome_type = dplyr::coalesce(outcome_types[id.outcome], "binary"),
    effect_label = ifelse(outcome_type == "binary", "OR", "Beta"),
    est_display  = ifelse(outcome_type == "binary", exp(b), b),
    ci_low       = ifelse(outcome_type == "binary", exp(b - 1.96 * se), b - 1.96 * se),
    ci_high      = ifelse(outcome_type == "binary", exp(b + 1.96 * se), b + 1.96 * se),
    est_ci       = sprintf("%.3f (%.3f–%.3f)", est_display, ci_low, ci_high),
    p_display    = dplyr::case_when(
      is.na(pval)    ~ NA_character_,
      pval < 0.001   ~ "<0.001",
      pval < 0.01    ~ sprintf("%.3f", pval),
      TRUE           ~ sprintf("%.3f", pval)
    )
  ) |>
  dplyr::select(
    id.outcome, outcome, method, b, se, pval, nsnp, analysis_type,
    left_out_cohort, n_snps_remaining, n_cohorts_remaining,
    outcome_type, effect_label, est_display, ci_low, ci_high,
    est_ci, p_display
  )

readr::write_csv(all_results,
                 file.path(results_dir, "leave_one_cohort_results.csv"))
readr::write_csv(cohort_summary,
                 file.path(results_dir, "cohort_composition.csv"))

message("Sensitivity analyses complete. Results saved to: ", results_dir)
