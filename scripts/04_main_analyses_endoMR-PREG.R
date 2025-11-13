#!/usr/bin/env Rscript
###############################################################################
# Two-Sample MR analysis
# Input  : results/harmonised_rahmioglu_bpo.csv
# Output : CSVs and plots for IVW, Egger, WM, heterogeneity, pleiotropy, LOO,
#          and cohort leave-one-out (results saved under results/)
###############################################################################

### 1) SETUP ####################################################################
# Packages
required_pkgs <- c("TwoSampleMR", "MRPRESSO", "dplyr", "ggplot2",
                   "here", "readr", "data.table", "devtools")

safe_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "MRPRESSO") {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools", repos = "https://cloud.r-project.org")
      }
      devtools::install_github("rondolab/MR-PRESSO")
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_pkgs, safe_install))

# Paths
results_dir <- here::here("results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

harm_file    <- here::here("results", "harmonised_rahmioglu_bpo.csv")
# Kept for completeness (not used below but may be needed later)
path_outcome <- here::here("data", "OUTCOME_MR-PREG", "stu_out_dat.txt")

### 2) LOAD + RESTRICT OUTCOMES #################################################
# Read harmonised dataset
dat <- data.table::fread(harm_file)

# Minimal sanity checks
stopifnot(all(c("id.exposure", "beta.exposure", "beta.outcome") %in% colnames(dat)))
stopifnot("outcome" %in% colnames(dat))

# Dictionary of human-readable outcome names (exactly 29)
outcome_labels <- c(
  # Placental outcomes
  finngen_R12_O15_PLAC_PRAEVIA                      = "Placenta praevia",
  finngen_R12_O15_PLAC_DISORD                       = "Placental disorders",
  finngen_R12_O15_PLAC_PREMAT_SEPAR                 = "Premature placental separation",
  
  # Caesarean delivery
  el_cs                                            = "Elective caesarean section",
  em_cs                                            = "Emergency caesarean section",
  cs                                               = "Caesarean section",
  
  # Apgar scores
  lowapgar1                                        = "Low Apgar score at 1 min",
  lowapgar5                                        = "Low Apgar score at 5 min",
  
  # Labor and delivery complications
  rup_memb                                         = "Premature rupture of membranes",
  induction                                        = "Labour induction",
  
  # Gestational age and timing
  ga_all                                           = "Gestational age",
  pretb_all                                        = "Preterm birth (any)",
  vpretb_all                                       = "Very preterm birth",
  posttb_all                                       = "Post-term birth",
  
  # Birth weight outcomes
  hbw_all                                          = "High birthweight (>4000g)",
  lbw_all                                          = "Low birthweight (<2500g)",
  sga                                              = "Small for gestational age",
  
  # Maternal health
  depr_subsamp                                     = "Postpartum Depression",
  anaemia_preg_all                                 = "Pregnancy anemia",
  
  # Pregnancy complications
  gdm_subsamp                                      = "Gestational diabetes mellitus",
  hdp_subsamp                                      = "Hypertensive disorders of pregnancy",
  gh_subsamp                                       = "Gestational hypertension",
  pe_subsamp                                       = "Preeclampsia",
  
  # Neonatal outcomes
  nicu                                             = "NICU admission",
  sb_subsamp                                       = "Stillbirth",
  
  # Hemorrhage and bleeding
  Antepartum_bleeding_filtered                     = "Antepartum bleeding",
  Postpartum_hemorrhage_filtered                   = "Postpartum hemorrhage",
  Postpartum_hemorrhage_due_to_atony_filtered      = "PPH due to atony",
  Postpartum_hemorrhage_due_to_retained_placenta_filtered = "PPH due to retained placenta"
)

# Restrict to the 29 target outcomes only
dat <- dat[dat$outcome %in% names(outcome_labels), , drop = FALSE]

# Attach readable labels
labels_df <- data.frame(
  outcome      = names(outcome_labels),
  outcome_full = unname(outcome_labels),
  stringsAsFactors = FALSE
)

dat <- merge(dat, labels_df, by = "outcome", all.x = TRUE, sort = FALSE)
dat$outcome_full[is.na(dat$outcome_full)] <- dat$outcome[is.na(dat$outcome_full)]
dat$id.outcome <- dat$outcome

# Basic cleaning
dat <- dat %>%
  dplyr::filter(
    !is.na(beta.outcome),
    !is.na(se.outcome),
    !is.infinite(se.outcome),
    !is.na(pval.outcome),
    se.outcome > 0
  )

# Show the final list of outcomes (for a quick visual check)
message("Outcomes included (n = ", length(unique(dat$outcome_full)), "):")
print(sort(unique(dat$outcome_full)))

### 3) HELPERS ##################################################################
run_mr_methods <- function(data, methods) {
  # data must be a harmonised TwoSampleMR format tibble/data.frame
  out <- mr(data, method_list = methods)
  as.data.frame(out)
}

export_csv <- function(df, stem) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    warning("Nothing to export for: ", stem)
    return(invisible(NULL))
  }
  
  # Ensure readable labels are present when an 'outcome' column exists
  if ("outcome" %in% colnames(df) && !"outcome_full" %in% colnames(df)) {
    df <- merge(df, labels_df, by = "outcome", all.x = TRUE, sort = FALSE)
    df$outcome_full[is.na(df$outcome_full)] <- df$outcome[is.na(df$outcome_full)]
  }
  
  # Format numeric columns: round normally; scientific for very small values
  num_cols <- vapply(df, is.numeric, logical(1L))
  if (any(num_cols)) {
    df[num_cols] <- lapply(df[num_cols], function(x) {
      sapply(x, function(val) {
        if (is.na(val)) return(NA)
        if (abs(val) < 1e-3) format(val, scientific = TRUE, digits = 3)
        else round(val, 3)
      })
    })
  }
  
  out_path <- file.path(results_dir, paste0(stem, ".csv"))
  write.csv(df, out_path, row.names = FALSE)
  message("Wrote: ", out_path)
}

### 4) MAIN MR ESTIMATES ########################################################
# IVW
ivw_res <- run_mr_methods(dat, "mr_ivw")
ivw_res_raw <- ivw_res

# MR-Egger
egger_res <- run_mr_methods(dat, "mr_egger_regression")

# Weighted median
wm_res <- run_mr_methods(dat, "mr_weighted_median")

# All default MR methods (IVW, Egger, WM, simple/weighted modes, etc.)
all_res <- mr(dat)

### 4b) MULTIPLE TESTING CORRECTION (FDR) #######################################
ivw_res$qval   <- p.adjust(ivw_res$pval,   method = "fdr")
egger_res$qval <- p.adjust(egger_res$pval, method = "fdr")
wm_res$qval    <- p.adjust(wm_res$pval,    method = "fdr")
all_res$qval   <- p.adjust(all_res$pval,   method = "fdr")

# Export
export_csv(ivw_res,   "ivw_results")
export_csv(egger_res, "egger_results")
export_csv(wm_res,    "weighted_median_results")
export_csv(all_res,   "all_mr_methods")
