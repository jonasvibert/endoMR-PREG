#!/usr/bin/env Rscript
################################################################################
# Script: 03_harmonise_data_endoMR-PREG.R
# Project: endoMR-PREG
# Purpose:
#   Harmonise Rahmioglu endometriosis instruments with pregnancy outcome GWAS
#   - Exposure: Rahmioglu endometriosis instruments (clumped SNPs)
#   - Outcomes: MR-PREG, Westergaard PPH, FinnGen R12
#
# Output:
#   results/harmonised_rahmioglu_bpo.csv
################################################################################

### 1) Setup ###################################################################

# Packages ---------------------------------------------------------------------
required_pkgs <- c(
  "TwoSampleMR",
  "dplyr",
  "data.table",
  "here"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

# Directories ------------------------------------------------------------------
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# 2) Exposure data: Rahmioglu endometriosis instruments                        #
################################################################################

message("\n=== Loading exposure instruments (Rahmioglu endometriosis) ===")

# Preferred file: full clumped object exported in script 01
exp_file_candidates <- c(
  file.path(results_dir, "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv"),
  file.path(results_dir, "endometriosis_clumped_snps.tsv")  # legacy name
)

exp_file <- exp_file_candidates[file.exists(exp_file_candidates)][1]

if (is.na(exp_file)) {
  stop(
    "No exposure instrument file found. Expected one of:\n  - ",
    paste(exp_file_candidates, collapse = "\n  - "),
    "\nRun 01_select_instruments_endoMR-PREG.R first."
  )
}

message("Exposure file: ", basename(exp_file))

exposure_dat <- data.table::fread(exp_file, data.table = FALSE)

# Sanity checks ---------------------------------------------------------------
if (!"id.exposure" %in% colnames(exposure_dat)) {
  warning("Column 'id.exposure' not found in exposure dataset. ",
          "Harmonisation will still work but some metadata may be missing.")
}

message("Exposure dataset: ", nrow(exposure_dat), " rows x ",
        ncol(exposure_dat), " columns")
message("  - Unique SNPs: ", dplyr::n_distinct(exposure_dat$SNP))

################################################################################
# 3) Outcome data: pregnancy outcomes                                          #
################################################################################

message("\n=== Loading pregnancy outcome GWAS ===")

out_file <- file.path(results_dir, "formatted_all_pregnancy_outcomes.tsv")
if (!file.exists(out_file)) {
  stop(
    "Outcome file 'formatted_all_pregnancy_outcomes.tsv' not found in ",
    results_dir, ".\nRun 02_prepare_outcomes_endoMR-PREG.R first."
  )
}

message("Outcome file: ", basename(out_file))

outcome_dat <- data.table::fread(out_file, data.table = FALSE)

# Sanity checks ---------------------------------------------------------------
required_cols <- c("SNP", "id.outcome", "outcome")
missing_cols  <- setdiff(required_cols, colnames(outcome_dat))

if (length(missing_cols)) {
  stop("Outcome dataset is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

message("Outcome dataset: ", nrow(outcome_dat), " rows x ",
        ncol(outcome_dat), " columns")
message("  - Unique outcomes: ", dplyr::n_distinct(outcome_dat$outcome))
message("  - Unique SNPs:     ", dplyr::n_distinct(outcome_dat$SNP))

################################################################################
# 4) Harmonisation                                                             #
################################################################################

message("\n=== Harmonising exposure and outcome datasets ===")
message("Harmonisation settings: action = 2 (drop ambiguous palindromic SNPs).")

dat <- TwoSampleMR::harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat  = outcome_dat,
  action       = 2
)

# Quick summary of harmonisation ----------------------------------------------
message("\nHarmonisation summary:")

if ("action" %in% colnames(dat)) {
  message("  - action codes (alignment status):")
  print(table(dat$action, useNA = "ifany"))
}

if ("remove" %in% colnames(dat)) {
  message("  - remove flag (SNPs dropped after harmonisation):")
  print(table(dat$remove, useNA = "ifany"))
}

if ("palindromic" %in% colnames(dat)) {
  message("  - palindromic SNPs:")
  print(table(dat$palindromic, useNA = "ifany"))
}

if ("outcome" %in% colnames(dat)) {
  message("  - SNP counts per outcome (after harmonisation):")
  print(table(dat$outcome))
}

if ("mr_keep" %in% colnames(dat)) {
  message("  - palindromic SNPs suggested for keep  per outcome (after harmonisation):")
  print(table(dat$mr_keep))
}

dat <- dat %>% filter(mr_keep==TRUE)

message("Final harmonised dataset: ", nrow(dat), " rows.")

# Optional: inspect first rows (can be commented out in production) ----------
# print(head(dat))
# str(dat)

################################################################################
# 5) Save harmonised dataset                                                   #
################################################################################

out_harmonised <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")

write.csv(dat, out_harmonised, row.names = FALSE)
message("\nSaved harmonised dataset to: ", out_harmonised)
