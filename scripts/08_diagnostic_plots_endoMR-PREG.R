#!/usr/bin/env Rscript
###############################################################################
# Script: 08_diagnostic_plots_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   Generate diagnostic plots for all outcomes with IVW p-value < 0.05:
#     - Scatter plot (mr_scatter_plot)
#     - Single-SNP forest plot (mr_forest_plot)
#     - Leave-one-out plot (mr_leaveoneout_plot)
#     - Funnel plot (mr_funnel_plot)
#
# Inputs:
#   - results/harmonised_rahmioglu_bpo.csv   (harmonised SNP-outcome data)
#   - results/ivw_results.csv                (IVW MR results with p-values)
#
# Outputs (saved to results/plots/diagnostics/):
#   - 08_<outcome>_scatter.(png|pdf)
#   - 08_<outcome>_forest.(png|pdf)
#   - 08_<outcome>_leaveoneout.(png|pdf)
#   - 08_<outcome>_funnel.(png|pdf)
#
# Author: Jonas Vibert
###############################################################################

### 1) SETUP ###################################################################

required_pkgs <- c(
  "TwoSampleMR",
  "dplyr",
  "ggplot2",
  "here",
  "readr",
  "data.table",
  "stringr",
  "purrr"
)

safe_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_pkgs, safe_install))

# Paths
results_dir <- here::here("results")
plots_dir   <- file.path(results_dir, "plots")
diag_dir    <- file.path(plots_dir, "diagnostics")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(diag_dir,    showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

read_csv_safe <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

message("=== 08_diagnostic_plots_endoMR-PREG.R ===")

### Outcome labels (same mapping as in 06/07 scripts) ##########################

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


### 2) LOAD HARMONISED DATA + IVW RESULTS ######################################

harm_file        <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")
ivw_results_path <- file.path(results_dir, "ivw_results.csv")

harm_df     <- data.table::fread(harm_file)
ivw_results <- read_csv_safe(ivw_results_path)

# Basic checks
required_cols_harm <- c(
  "outcome", "id.exposure", "id.outcome",
  "beta.exposure", "beta.outcome",
  "se.exposure", "se.outcome", "pval.outcome"
)
stopifnot(all(required_cols_harm %in% names(harm_df)))

required_cols_ivw <- c("outcome", "method", "pval")
stopifnot(all(required_cols_ivw %in% names(ivw_results)))

if (!"outcome_full" %in% names(ivw_results)) {
  ivw_results$outcome_full <- ivw_results$outcome
}

message("Loaded harmonised data: ", nrow(harm_df), " rows")
message("Loaded IVW results: ", nrow(ivw_results), " rows")

### 3) SELECT OUTCOMES WITH IVW p < 0.05 #######################################

ivw_results <- ivw_results %>%
  mutate(p_num = suppressWarnings(as.numeric(pval)))

ivw_sig <- ivw_results %>%
  filter(
    method == "Inverse variance weighted",
    !is.na(p_num),
    p_num < 0.05
  ) %>%
  arrange(p_num)

if (!nrow(ivw_sig)) {
  message("No outcomes with IVW p < 0.05. No diagnostic plots produced.")
  quit(save = "no")
}

sig_outcomes <- ivw_sig %>%
  distinct(outcome, outcome_full, p_num)

message("Significant outcomes (IVW p < 0.05):")
print(sig_outcomes)

### 4) FUNCTION TO GENERATE ALL PLOTS FOR ONE OUTCOME ##########################

run_plots_for_outcome <- function(outcome_id,
                                  outcome_label,
                                  dat_all,
                                  out_dir) {
  message("\n--- Outcome: ", outcome_id, " (", outcome_label, ") ---")
  
  # Extract relevant rows
  dat_out <- dat_all[dat_all$outcome == outcome_id, , drop = FALSE]
  
  if (nrow(dat_out) == 0L) {
    message("No harmonised data for outcome: ", outcome_id)
    return(invisible(NULL))
  }
  
  # Clean invalid SNPs
  dat_out <- dat_out %>%
    filter(
      !is.na(beta.exposure),
      !is.na(beta.outcome),
      !is.na(se.exposure),
      !is.na(se.outcome),
      se.exposure > 0,
      se.outcome  > 0
    )
  
  if (nrow(dat_out) == 0L) {
    message("No valid SNPs after filtering for ", outcome_id)
    return(invisible(NULL))
  }
  
  # Require ≥ 3 SNPs
  n_snps <- dplyr::n_distinct(dat_out$SNP)
  if (n_snps < 3L) {
    message("Too few SNPs (n = ", n_snps, "). Skipping.")
    return(invisible(NULL))
  }
  
  # MR methods for scatter plot
  mr_methods_main <- c(
    "mr_egger_regression",
    "mr_ivw",
    "mr_weighted_median",
    "mr_weighted_mode"
  )
  
  # MR results
  res <- tryCatch(
    TwoSampleMR::mr(dat_out, method_list = mr_methods_main),
    error = function(e) {
      message("mr() failed for ", outcome_id, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(res)) return(invisible(NULL))
  
  # Single SNP MR
  res_single <- tryCatch(
    TwoSampleMR::mr_singlesnp(dat_out),
    error = function(e) {
      message("mr_singlesnp() failed for ", outcome_id)
      return(NULL)
    }
  )
  if (is.null(res_single)) return(invisible(NULL))
  
  # Leave-one-out MR
  res_loo <- tryCatch(
    TwoSampleMR::mr_leaveoneout(dat_out, method = TwoSampleMR::mr_ivw),
    error = function(e) {
      message("mr_leaveoneout() failed for ", outcome_id)
      return(NULL)
    }
  )
  if (is.null(res_loo)) return(invisible(NULL))
  
  # Filename prefix
  label_short <- ifelse(
    is.na(outcome_label) | outcome_label == "",
    outcome_id,
    outcome_label
  )
  prefix <- paste0("08_", safe_name(label_short), "_")
  
  ### 4.1 Scatter plot ----------------------------------------------------------
  p_scatter_list <- tryCatch(
    TwoSampleMR::mr_scatter_plot(res, dat_out),
    error = function(e) NULL
  )
  
  if (!is.null(p_scatter_list) && length(p_scatter_list) > 0) {
    p <- p_scatter_list[[1]] + ggtitle(outcome_label)
    ggsave(file.path(out_dir, paste0(prefix, "scatter.png")), p, width = 7, height = 7, dpi = 300)
    ggsave(file.path(out_dir, paste0(prefix, "scatter.pdf")), p, width = 7, height = 7)
  }
  
  ### 4.2 Single SNP forest plot ------------------------------------------------
  p_forest_list <- tryCatch(
    TwoSampleMR::mr_forest_plot(res_single),
    error = function(e) NULL
  )
  
  if (!is.null(p_forest_list) && length(p_forest_list) > 0) {
    p <- p_forest_list[[1]] + ggtitle(outcome_label)
    ggsave(file.path(out_dir, paste0(prefix, "forest.png")), p, width = 7, height = 7, dpi = 300)
    ggsave(file.path(out_dir, paste0(prefix, "forest.pdf")), p, width = 7, height = 7)
  }
  
  ### 4.3 Leave-one-out plot ----------------------------------------------------
  p_loo_list <- tryCatch(
    TwoSampleMR::mr_leaveoneout_plot(res_loo),
    error = function(e) NULL
  )
  
  if (!is.null(p_loo_list) && length(p_loo_list) > 0) {
    p <- p_loo_list[[1]] + ggtitle(outcome_label)
    ggsave(file.path(out_dir, paste0(prefix, "leaveoneout.png")), p, width = 7, height = 7, dpi = 300)
    ggsave(file.path(out_dir, paste0(prefix, "leaveoneout.pdf")), p, width = 7, height = 7)
  }
  
  ### 4.4 Funnel plot -----------------------------------------------------------
  p_funnel_list <- tryCatch(
    TwoSampleMR::mr_funnel_plot(res_single),
    error = function(e) NULL
  )
  
  if (!is.null(p_funnel_list) && length(p_funnel_list) > 0) {
    p <- p_funnel_list[[1]] + ggtitle(outcome_label)
    ggsave(file.path(out_dir, paste0(prefix, "funnel.png")), p, width = 7, height = 7, dpi = 300)
    ggsave(file.path(out_dir, paste0(prefix, "funnel.pdf")), p, width = 7, height = 7)
  }
  
  invisible(NULL)
}

### 5) LOOP ACROSS SIGNIFICANT OUTCOMES ########################################

for (i in seq_len(nrow(sig_outcomes))) {
  run_plots_for_outcome(
    outcome_id    = sig_outcomes$outcome[i],
    outcome_label = sig_outcomes$outcome_full[i],
    dat_all       = harm_df,
    out_dir       = diag_dir
  )
}

### 6) SUMMARY ##################################################################

message("\n=== DIAGNOSTIC PLOTS COMPLETED ===")
message("All plots saved to: ", diag_dir)
