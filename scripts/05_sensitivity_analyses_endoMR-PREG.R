#!/usr/bin/env Rscript
################################################################################
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
################################################################################

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

set.seed(42)

dat       <- data.table::fread(harm_file) |> as.data.frame()
harm_data <- dat  # reuse for exposure reconstruction

stopifnot(all(c("SNP", "beta.exposure", "se.exposure",
                "effect_allele.exposure", "other_allele.exposure",
                "eaf.exposure", "id.exposure", "exposure") %in% names(harm_data)))
stopifnot(all(c("beta.outcome", "se.outcome", "pval.outcome",
                "effect_allele.outcome", "other_allele.outcome",
                "eaf.outcome", "id.outcome", "outcome") %in% names(dat)))

message("Loaded harmonised dataset: ", nrow(dat), " rows.")

### 1b) RESTRICT TO THE 30 MAIN OUTCOMES #######################################

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

# Restrict harmonised dataset to the 30 outcomes
dat <- dat[dat$outcome %in% vars_keep, , drop = FALSE]

# Label mapping (if needed later)
labels_df <- data.frame(
  outcome      = names(outcome_labels),
  outcome_full = unname(outcome_labels),
  stringsAsFactors = FALSE
)

message("Outcomes included in sensitivity analyses (n = ",
        length(unique(dat$outcome)), "):")
print(sort(unique(dat$outcome)))

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

### 3) STANDARD SENSITIVITY ANALYSES ###########################################

message("=== Standard TwoSampleMR sensitivity analyses (30 outcomes) ===")

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
    n_snps          = dplyr::n_distinct(SNP),
    n_phenotypes    = dplyr::n_distinct(Phenotype),
    mean_samplesize = mean(samplesize, na.rm = TRUE),
    .groups         = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_snps))

message("Available cohorts: ", paste(cohort_summary$study, collapse = ", "))
if (nrow(cohort_summary) < 2) {
  stop("Need ≥ 2 cohorts for leave-one-cohort-out analysis.")
}

# 4.3 Restrict to the same 30 outcomes and define types ------------------------

# Re-use the same 30 outcomes
vars_keep_loo <- vars_keep

# Outcome types: ga_all and zbw_all are continuous, all others binary
outcome_types <- c(
  Antepartum_bleeding                       = "binary",
  Postpartum_hemorrhage                     = "binary",
  Postpartum_hemorrhage_due_to_atony        = "binary",
  Postpartum_hemorrhage_due_to_retained_placenta = "binary",
  finngen_R12_O15_PLAC_PRAEVIA              = "binary",
  finngen_R12_O15_PLAC_DISORD               = "binary",
  finngen_R12_O15_PLAC_PREMAT_SEPAR         = "binary",
  
  rup_memb                                  = "binary",
  
  pretb_all                                 = "binary",
  vpretb_all                                = "binary",
  
  ga_all                                    = "continuous",
  sga                                       = "binary",
  lbw_all                                   = "binary",
  hbw_all                                   = "binary",
  lga                                       = "binary",
  zbw_all                                   = "continuous",
  
  lowapgar1                                 = "binary",
  lowapgar5                                 = "binary",
  nicu                                      = "binary",
  sb_subsamp                                = "binary",

  anaemia_preg_all                          = "binary",
  gdm_subsamp                               = "binary",
  gh_subsamp                                = "binary",
  hdp_subsamp                               = "binary",
  pe_subsamp                                = "binary",
  depr_subsamp                              = "binary",
  
  induction                                 = "binary",
  posttb_all                                = "binary",
  el_cs                                     = "binary",
  em_cs                                     = "binary"
)

# Restrict raw MR-PREG outcome data to the 30 outcomes
mr_data <- mr_data |>
  dplyr::filter(Phenotype %in% vars_keep_loo) |>
  dplyr::mutate(
    id.outcome    = Phenotype,
    outcome_label = dplyr::coalesce(outcome_labels[Phenotype], Phenotype),
    outcome_type  = dplyr::coalesce(outcome_types[Phenotype], "binary")
  )

message("Phenotypes in leave-one-cohort-out analysis: ",
        paste(sort(unique(mr_data$Phenotype)), collapse = ", "))

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
    beta.outcome          = beta,
    se.outcome            = se,
    pval.outcome          = pval,
    effect_allele.outcome = effect_allele,
    other_allele.outcome  = other_allele,
    eaf.outcome           = eaf,
    samplesize.outcome    = samplesize,
    ncase.outcome         = ncase,
    ncontrol.outcome      = ncontrol,
    id.outcome,
    outcome               = outcome_label,
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

readr::write_csv(
  all_results,
  file.path(results_dir, "leave_one_cohort_results.csv")
)
readr::write_csv(
  cohort_summary,
  file.path(results_dir, "cohort_composition.csv")
)

### 5) MR-PRESSO GLOBAL AND OUTLIER TESTS ######################################

message("=== MR-PRESSO global and outlier tests (30 outcomes) ===")

# Ensure MRPRESSO is available -------------------------------------------------
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::install_github("rondolab/MR-PRESSO")
}
suppressPackageStartupMessages(library(MRPRESSO))

# We assume `dat` is the harmonised dataset restricted to the 30 outcomes
# with at least: SNP, outcome, beta.exposure, se.exposure, beta.outcome, se.outcome

# Split by outcome
dat_list <- split(dat, dat$outcome)

run_mr_presso_for_outcome <- function(df, outcome_id) {
  n_snps <- length(unique(df$SNP))
  message("[MR-PRESSO] Outcome: ", outcome_id, " | SNPs: ", n_snps)
  
  # Need at least 3 SNPs
  if (n_snps < 3) {
    message("[MR-PRESSO] Skipping ", outcome_id, ": <3 SNPs.")
    return(NULL)
  }
  
  # Prepare data frame in the format expected by mr_presso()
  presso_dat <- data.frame(
    SNP          = df$SNP,
    BetaExposure = df$beta.exposure,
    BetaOutcome  = df$beta.outcome,
    SdExposure   = df$se.exposure,
    SdOutcome    = df$se.outcome,
    stringsAsFactors = FALSE
  )
  
  # Drop rows with missing / non-finite values
  presso_dat <- presso_dat[
    is.finite(presso_dat$BetaExposure) &
      is.finite(presso_dat$BetaOutcome) &
      is.finite(presso_dat$SdExposure) &
      is.finite(presso_dat$SdOutcome),
  ]
  
  if (nrow(presso_dat) < 3) {
    message("[MR-PRESSO] Skipping ", outcome_id, ": <3 valid SNPs after filtering.")
    return(NULL)
  }
  
  # Run MR-PRESSO (NbDistribution can be increased later if needed)
  message("[MR-PRESSO] Running MR-PRESSO for ", outcome_id, " ...")
  res <- tryCatch(
    {
      MRPRESSO::mr_presso(
        BetaOutcome    = "BetaOutcome",
        BetaExposure   = "BetaExposure",
        SdOutcome      = "SdOutcome",
        SdExposure     = "SdExposure",
        OUTLIERtest    = TRUE,
        DISTORTIONtest = TRUE,
        data           = presso_dat,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )
    },
    error = function(e) {
      message("[MR-PRESSO] Failed for ", outcome_id, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(res)) return(NULL)
  
  main_res   <- res[["Main MR results"]]
  presso_res <- res[["MR-PRESSO results"]]
  distortion <- res[["Distortion Test"]]
  
  global_p     <- NA_real_
  n_outliers   <- 0L
  outlier_snps <- NA_character_
  causal_orig  <- NA_real_
  se_orig      <- NA_real_
  pval_orig    <- NA_real_
  causal_corr  <- NA_real_
  se_corr      <- NA_real_
  pval_corr    <- NA_real_

  # --- Global test p-value ---------------------------------------------------
  # MRPRESSO returns Global Test as a named list ($RSSobs, $Pvalue),
  # not a data frame, so colnames() fails. Handle both formats.
  if (!is.null(presso_res) && !is.null(presso_res[["Global Test"]])) {
    gt       <- presso_res[["Global Test"]]
    pval_val <- if (is.data.frame(gt)) gt[["Pvalue"]][1] else gt[["Pvalue"]]
    if (!is.null(pval_val) && length(pval_val) == 1) {
      global_p <- suppressWarnings(as.numeric(pval_val))
    }
  }

  # --- Outlier test: identify outlier SNPs -----------------------------------
  # MRPRESSO Outlier Test returns a data frame with one row per SNP and a
  # "Pvalue" column. Outliers are those where Pvalue < SignifThreshold.
  if (!is.null(presso_res) && !is.null(presso_res[["Outlier Test"]])) {
    ot <- presso_res[["Outlier Test"]]
    if (is.data.frame(ot) && "Pvalue" %in% colnames(ot)) {
      # Use significance threshold matching what was passed to mr_presso()
      sig_thresh   <- 0.05
      pvals_num    <- suppressWarnings(as.numeric(
        gsub("<", "", as.character(ot[["Pvalue"]]))
      ))
      outlier_rows <- which(!is.na(pvals_num) & pvals_num < sig_thresh)
      if (length(outlier_rows) > 0) {
        outlier_vec  <- presso_dat$SNP[outlier_rows]
        outlier_vec  <- unique(outlier_vec[!is.na(outlier_vec)])
        n_outliers   <- length(outlier_vec)
        outlier_snps <- paste(outlier_vec, collapse = ";")
      }
    }
  }

  # --- Causal estimates from Main MR results ---------------------------------
  # main_res has rows "Raw" and (if outliers) "Outlier-corrected".
  # Columns: "MR Analysis", "Causal Estimate", "Sd", "T-stat", "P-value"
  if (!is.null(main_res) && is.data.frame(main_res) &&
      "Causal Estimate" %in% colnames(main_res)) {
    raw_row  <- main_res[grepl("Raw",              main_res[["MR Analysis"]], ignore.case = TRUE), ]
    corr_row <- main_res[grepl("Outlier-corrected", main_res[["MR Analysis"]], ignore.case = TRUE), ]

    if (nrow(raw_row) > 0) {
      causal_orig <- suppressWarnings(as.numeric(raw_row[1, "Causal Estimate"]))
      if ("Sd"      %in% colnames(main_res)) se_orig   <- suppressWarnings(as.numeric(raw_row[1, "Sd"]))
      if ("P-value" %in% colnames(main_res)) pval_orig <- suppressWarnings(as.numeric(raw_row[1, "P-value"]))
    }
    if (nrow(corr_row) > 0) {
      causal_corr <- suppressWarnings(as.numeric(corr_row[1, "Causal Estimate"]))
      if ("Sd"      %in% colnames(main_res)) se_corr   <- suppressWarnings(as.numeric(corr_row[1, "Sd"]))
      if ("P-value" %in% colnames(main_res)) pval_corr <- suppressWarnings(as.numeric(corr_row[1, "P-value"]))
    }
  }

  data.frame(
    outcome_id         = outcome_id,
    outcome_label      = if (!is.null(outcome_labels[outcome_id]))
      unname(outcome_labels[outcome_id]) else outcome_id,
    n_snps_input       = nrow(presso_dat),
    global_p           = global_p,
    n_outliers         = n_outliers,
    outlier_snps       = outlier_snps,
    causal_estimate    = causal_orig,
    se_raw             = se_orig,
    pval_raw           = pval_orig,
    corrected_estimate = causal_corr,
    se_corrected       = se_corr,
    pval_corrected     = pval_corr,
    stringsAsFactors   = FALSE
  )
}

# Run MR-PRESSO for each outcome
mrpresso_summary_list <- mapply(
  FUN        = run_mr_presso_for_outcome,
  df         = dat_list,
  outcome_id = names(dat_list),
  SIMPLIFY   = FALSE
)

mrpresso_summary <- dplyr::bind_rows(mrpresso_summary_list)

# Export results
readr::write_csv(
  mrpresso_summary,
  file.path(results_dir, "mr_presso_results.csv")
)

message("MR-PRESSO analyses complete. Results saved to: ", results_dir)
message("Sensitivity analyses complete. Results saved to: ", results_dir)
