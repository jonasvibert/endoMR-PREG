#!/usr/bin/env Rscript
################################################################################
# Script: 05.2_effective_tests_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   [R4.4] Reviewer #4 ("the effective number of independent tests is closer
#   to 12 to 15... estimate the effective test count from the eigenvalues of
#   the outcome z-statistic correlation matrix (Li and Ji, or Nyholt)"):
#   estimate the effective number of independent tests (Meff) among the 30
#   main outcomes, using the eigenvalues of the outcome z-statistic
#   correlation matrix (computed across the shared instrument SNPs), and
#   re-assess significance of the main IVW results against a Bonferroni
#   threshold based on Meff instead of the nominal 30.
#
# Method:
#   Li & Ji (2005, Heredity) eigenvalue-based effective number of tests:
#     Meff = sum_i f(lambda_i),  f(lambda) = 1                if lambda >= 1
#                                 f(lambda) = lambda            if lambda <  1
#   (equivalent to floor(lambda)=0 for lambda<1, so f(lambda)=lambda-0=lambda)
#   Nyholt (2004) is reported alongside for comparison (as suggested by the
#   reviewer): Meff = M - sum_i [I(lambda_i > 1) * (lambda_i - 1)]
#   Both operate on the same correlation matrix of outcome z-statistics.
#
# Input:
#   - results/harmonised_rahmioglu_bpo.csv (script 03)
#   - results/ivw_results.csv (script 04)
#
# Outputs (saved under results/):
#   - outcome_zscore_correlation_matrix.csv
#   - effective_number_of_tests.csv (Meff estimates + resulting threshold)
#   - ivw_results_with_effective_test_correction.csv
################################################################################

### 1) SETUP ####################################################################

required_pkgs <- c("dplyr", "tidyr", "data.table", "here")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

results_dir <- here::here("results")

harm_file <- file.path(results_dir, "harmonised_rahmioglu_bpo.csv")
ivw_file  <- file.path(results_dir, "ivw_results.csv")

if (!file.exists(harm_file)) stop("Missing file: ", harm_file, ". Run script 03 first.")
if (!file.exists(ivw_file))  stop("Missing file: ", ivw_file,  ". Run script 04 first.")

vars_keep <- c(
  "Antepartum_bleeding", "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta",
  "finngen_R12_O15_PLAC_PRAEVIA", "finngen_R12_O15_PLAC_DISORD",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR", "rup_memb",
  "pretb_all", "vpretb_all", "ga_all", "sga", "lbw_all", "hbw_all", "lga",
  "zbw_all", "lowapgar1", "lowapgar5", "nicu", "sb_subsamp",
  "anaemia_preg_all", "gdm_subsamp", "gh_subsamp", "hdp_subsamp",
  "pe_subsamp", "depr_subsamp", "induction", "posttb_all", "el_cs", "em_cs"
)

### 2) BUILD THE SNP x OUTCOME Z-STATISTIC MATRIX ##############################

dat <- data.table::fread(harm_file, data.table = FALSE)
dat <- dat[dat$outcome %in% vars_keep, ]

dat$z <- dat$beta.outcome / dat$se.outcome

z_wide <- dat |>
  dplyr::select(SNP, outcome, z) |>
  dplyr::distinct(SNP, outcome, .keep_all = TRUE) |>
  tidyr::pivot_wider(names_from = outcome, values_from = z) |>
  tibble::column_to_rownames("SNP")

message("Z-statistic matrix: ", nrow(z_wide), " SNPs x ", ncol(z_wide), " outcomes.")
message("Missing outcomes (expected 30): ",
        paste(setdiff(vars_keep, colnames(z_wide)), collapse = ", "))

### 3) CORRELATION MATRIX + EIGENVALUES #########################################

cor_mat <- cor(as.matrix(z_wide), use = "pairwise.complete.obs")

write.csv(
  cor_mat,
  file.path(results_dir, "outcome_zscore_correlation_matrix.csv")
)

eig <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values
eig <- pmax(eig, 0)  # clip tiny negative numerical artefacts to 0
M   <- length(eig)

message("Eigenvalues of the ", M, "x", M, " outcome correlation matrix computed.")
message("Sum of eigenvalues (should equal M = ", M, "): ", round(sum(eig), 3))

### 4) EFFECTIVE NUMBER OF TESTS (LI & JI; NYHOLT FOR COMPARISON) ##############

meff_liji <- sum(ifelse(eig >= 1, 1, eig))

meff_nyholt <- M - sum(ifelse(eig > 1, eig - 1, 0))

message(sprintf("Effective number of tests (Li & Ji, 2005):  %.2f (of M = %d nominal)", meff_liji, M))
message(sprintf("Effective number of tests (Nyholt, 2004):   %.2f (of M = %d nominal)", meff_nyholt, M))

alpha_nominal   <- 0.05 / M
alpha_liji      <- 0.05 / meff_liji
alpha_nyholt    <- 0.05 / meff_nyholt

meff_summary <- data.frame(
  method              = c("Nominal (M = 30)", "Li & Ji (2005)", "Nyholt (2004)"),
  M_effective         = c(M, meff_liji, meff_nyholt),
  bonferroni_alpha    = c(alpha_nominal, alpha_liji, alpha_nyholt)
)
write.csv(
  meff_summary,
  file.path(results_dir, "effective_number_of_tests.csv"),
  row.names = FALSE
)
print(meff_summary)

### 5) RE-ASSESS SIGNIFICANCE OF MAIN IVW RESULTS ###############################

ivw_res <- read.csv(ivw_file, stringsAsFactors = FALSE)

ivw_res <- ivw_res |>
  dplyr::mutate(
    sig_nominal_bonferroni = pval < alpha_nominal,
    sig_liji_bonferroni    = pval < alpha_liji,
    sig_nyholt_bonferroni  = pval < alpha_nyholt,
    sig_fdr                = qval < 0.05
  ) |>
  dplyr::arrange(pval)

write.csv(
  ivw_res,
  file.path(results_dir, "ivw_results_with_effective_test_correction.csv"),
  row.names = FALSE
)

message("\nOutcomes surviving each correction:")
message("  Nominal Bonferroni (M=30):        ",
        paste(ivw_res$outcome[ivw_res$sig_nominal_bonferroni], collapse = ", "))
message("  Li & Ji effective-test Bonferroni: ",
        paste(ivw_res$outcome[ivw_res$sig_liji_bonferroni], collapse = ", "))
message("  Nyholt effective-test Bonferroni:  ",
        paste(ivw_res$outcome[ivw_res$sig_nyholt_bonferroni], collapse = ", "))
message("  FDR (q < 0.05, as in main manuscript): ",
        paste(ivw_res$outcome[ivw_res$sig_fdr], collapse = ", "))

message("\n=== 05.2_effective_tests_endoMR-PREG.R completed ===")
