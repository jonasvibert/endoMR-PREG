#!/usr/bin/env Rscript
################################################################################
# Script: 05.5_mrlap_praevia_endoMR-PREG.R
# Project: endoMR-PREG
#
# REVISION LOG (revision-round1-reviewers)
#   [R4.2] Reviewer #4 ("Report the overlapping fraction per outcome and
#   overlap-aware estimates (MRlap); the bivariate LDSC intercept gives this
#   as a by-product of point 5."): overlap-aware, winner's-curse- and
#   weak-instrument-robust MR estimate for the endometriosis -> placenta
#   praevia association (the study's single positive finding, and therefore
#   the one for which residual sample overlap would most affect
#   interpretation), using the MRlap method (Mounier & Kutalik, 2023).
#
# Method: MRlap performs its own internal cross-trait LDSC (using the same
#   eur_w_ld_chr / w_hm3.snplist reference as script 05.4) to estimate the
#   genetic covariance intercept (a proxy for sample overlap), then reports
#   both the standard/uncorrected IVW estimate and an overlap + weak
#   instrument + winner's-curse corrected estimate, using its own
#   distance-pruned (not LD-clumped; no genotype reference panel available)
#   instrument selection at the default 5e-8 threshold, 500kb pruning
#   window. This instrument set is selected internally by MRlap from the
#   full exposure GWAS and will not be identical to the manuscript's primary
#   40/41-SNP clumped instrument (that is expected and standard for how
#   MRlap is applied/reported elsewhere).
#
# Data sources: same genome-wide files as script 05.3/05.4 (Rahmioglu 2023
#   full GWAS GCST90205183, GRCh37; FinnGen R12 O15_PLAC_PRAEVIA, GRCh38).
#   MRlap merges exposure/outcome on rsID internally so the build mismatch
#   is not an issue (same rationale as script 05.3). N is the fixed TOTAL
#   sample size for each trait (Rahmioglu: 475,160; FinnGen praevia:
#   223,001), per MRlap's own documentation ("If (at least) one of the
#   datasets is coming from a case-control GWAS, the sample size column
#   should correspond to the total sample size") - NOT per-SNP effective N.
#
# REVISION LOG (2026-07-27, verification following external review)
#   [FIX] N: replaced per-SNP effective N (1/(2*eaf*(1-eaf)*se^2)) for the
#   FinnGen outcome, and an undocumented RAHM_N=205904 for the exposure,
#   with a fixed total N per trait. This directly contradicted MRlap's
#   documented N requirement for case-control GWAS.
#
# REVISION LOG (2026-07-29)
#   [FIX] Rahmioglu N corrected 762,600 -> 475,160 (wrong accession,
#   GCST90258638 vs GCST90205183 actually read here - see script 05.6).
#   Affects h2_exposure, rg, and std_beta = Z/sqrt(N) effect sizes.
#   [FIX] Effect scale: MRlap's internal MR (see package source,
#   tidy_inputGWAS()/run_MR()) is run on Z/sqrt(N)-STANDARDISED effect
#   sizes (std_beta = Z/sqrt(N)), not on the raw log-odds scale. The
#   previous version of this script computed exp(corrected_effect) and
#   reported it as an odds ratio - this is not a valid transformation:
#   corrected_effect is a coefficient between two standardised effect-size
#   scales, not a log-odds ratio, and exponentiating it does not produce
#   an interpretable OR. corrected_OR/CI columns have been removed; the
#   observed/corrected effects are now reported on their native
#   (standardised) scale only, per MRlap's own recommended reporting
#   (compare observed_effect vs. corrected_effect on the same scale, use
#   p_difference to decide which to emphasise).
#   [ADD] h2_outcome_se (available in the package output but not
#   previously extracted), package version, and any warnings raised
#   during the run are now captured explicitly.
#
# Output:
#   - results/mrlap_praevia_results.csv
################################################################################

### 1) SETUP ####################################################################

required_pkgs <- c("data.table", "MRlap", "here")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

source(here::here("config", "config.R"))

results_dir <- here::here("results")
data_dir    <- here::here("data")

rahmioglu_gw_file    <- RAHMIOGLU_GW_FILE
finngen_praevia_file <- file.path(data_dir, "OUTCOME_FINNGEN", "finngen_R12_O15_PLAC_PRAEVIA")
ldsc_ref_dir         <- file.path(data_dir, "LDSC_REF", "eur_w_ld_chr")
hm3_file             <- file.path(data_dir, "LDSC_REF", "w_hm3.snplist")

if (!file.exists(rahmioglu_gw_file)) stop("Rahmioglu genome-wide file not found: ", rahmioglu_gw_file, " - set RAHMIOGLU_GW_FILE env var or edit config/config.R")
if (!file.exists(finngen_praevia_file)) stop("Missing file: ", finngen_praevia_file)

# N corrected to match the accession actually read above (GCST90205183 =
# 475,160), not GCST90258638 (762,600, a different meta-analysis stage from
# the same publication - see script 05.6 header for the full explanation).
RAHM_N   <- 475160L  # total N, GCST90205183: 23,779 cases + 450,668 controls
PRAEV_N  <- 223001L  # total N, FinnGen R12 placenta praevia: 1,815 cases + 221,186 controls

### 2) PREPARE EXPOSURE (RAHMIOGLU) - NEEDS SNP/CHR/POS/A1/A2/BETA/SE/N ########

message("Loading + preparing exposure (Rahmioglu) for MRlap...")
exposure <- data.table::fread(rahmioglu_gw_file)
data.table::setnames(exposure,
  c("chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error"),
  c("chr", "pos", "a1", "a2", "beta", "se"))
exposure$N <- RAHM_N
# MRlap needs an rsID column too (used to merge with outcome + hm3 aligning file)
message("Building rsID <-> CHR:BP map from eur_w_ld_chr...")
ldscore_files <- list.files(ldsc_ref_dir, pattern = "^[0-9]+\\.l2\\.ldscore\\.gz$", full.names = TRUE)
rsid_map <- data.table::rbindlist(lapply(ldscore_files, function(f) {
  data.table::fread(f, select = c("CHR", "SNP", "BP"))
}))
data.table::setnames(rsid_map, c("CHR", "SNP", "BP"), c("chr", "snp", "pos"))
exposure$chr <- as.character(exposure$chr)
rsid_map$chr <- as.character(rsid_map$chr)
exposure <- merge(exposure, rsid_map, by = c("chr", "pos"), all = FALSE)
exposure <- exposure[, c("snp", "chr", "pos", "a1", "a2", "beta", "se", "N")]
message("  Exposure prepared: ", nrow(exposure), " SNPs.")

### 3) PREPARE OUTCOME (FINNGEN PLACENTA PRAEVIA) ###############################

message("Loading + preparing outcome (FinnGen placenta praevia) for MRlap...")
outcome <- data.table::fread(finngen_praevia_file,
                              select = c("#chrom", "pos", "rsids", "ref", "alt",
                                         "beta", "sebeta", "af_alt"))
data.table::setnames(outcome, c("#chrom", "pos", "rsids", "ref", "alt", "beta", "sebeta", "af_alt"),
                      c("chr", "pos", "snp", "a2", "a1", "beta", "se", "eaf"))
outcome$snp <- sub(";.*$", "", outcome$snp)
outcome <- outcome[outcome$snp != "" & !is.na(outcome$snp) & !duplicated(outcome$snp), ]
outcome$N <- PRAEV_N  # fixed total N (see REVISION LOG) - not per-SNP effective N
outcome$chr <- as.character(outcome$chr)
outcome <- outcome[, c("snp", "chr", "pos", "a1", "a2", "beta", "se", "N")]
message("  Outcome prepared: ", nrow(outcome), " SNPs. Total N: ", PRAEV_N)

### 4) RUN MRLAP #################################################################

message("\nRunning MRlap (endometriosis -> placenta praevia)...")
message("MRlap package version: ", as.character(utils::packageVersion("MRlap")))

# Capture any warnings raised during the call (e.g. negative-h2 bootstrap
# warning from get_correction()) without interrupting execution.
mrlap_warnings <- character(0)
mrlap_res <- withCallingHandlers(
  MRlap::MRlap(
    exposure       = as.data.frame(exposure),
    exposure_name  = "Endometriosis_Rahmioglu2023",
    outcome        = as.data.frame(outcome),
    outcome_name   = "PlacentaPraevia_FinnGenR12",
    ld             = ldsc_ref_dir,
    hm3            = hm3_file,
    save_logfiles  = FALSE,
    verbose        = TRUE
  ),
  warning = function(w) {
    mrlap_warnings <<- c(mrlap_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (length(mrlap_warnings) == 0) {
  message("MRlap run completed with no warnings.")
} else {
  message("MRlap run raised ", length(mrlap_warnings), " warning(s):")
  for (w in mrlap_warnings) message("  - ", w)
}

mr_res <- mrlap_res$MRcorrection
lambda <- mrlap_res$LDSC

# --- [FIX] Effect scale: observed_effect/corrected_effect are on MRlap's
#     internal Z/sqrt(N)-standardised scale (see package source: std_beta =
#     Z/sqrt(N) in tidy_inputGWAS(); the IVW regression in run_MR() is run
#     on std_beta.exp/std_beta.out). This is NOT the raw log-odds scale, so
#     exp(corrected_effect) is not a valid odds ratio and is no longer
#     computed. Report observed vs. corrected on their native (shared,
#     comparable) standardised scale instead, per MRlap's own recommended
#     use of p_difference to decide which estimate to emphasise. ---
summary_tbl <- data.frame(
  metric = c(
    "mrlap_package_version",
    "n_instruments_used",
    "observed_effect_std", "observed_effect_se_std", "observed_effect_pval",
    "corrected_effect_std", "corrected_effect_se_std", "corrected_effect_pval",
    "p_difference",
    "egger_intercept_pval",
    "crosstrait_intercept", "crosstrait_intercept_se",
    "h2_exposure", "h2_exposure_se",
    "h2_outcome", "h2_outcome_se",
    "rg_exposure_outcome",
    "n_warnings"
  ),
  value = c(
    as.character(utils::packageVersion("MRlap")),
    mr_res$m_IVs,
    mr_res$observed_effect, mr_res$observed_effect_se, mr_res$observed_effect_p,
    mr_res$corrected_effect, mr_res$corrected_effect_se, mr_res$corrected_effect_p,
    mr_res$p_difference,
    mr_res$egger_intercept_p,
    lambda$int_crosstrait, lambda$int_crosstrait_se,
    lambda$h2_exp, lambda$h2_exp_se,
    lambda$h2_out, lambda$h2_out_se,
    lambda$rg,
    length(mrlap_warnings)
  )
)

write.csv(summary_tbl, file.path(results_dir, "mrlap_praevia_results.csv"), row.names = FALSE)
if (length(mrlap_warnings) > 0) {
  writeLines(mrlap_warnings, file.path(results_dir, "mrlap_praevia_warnings.txt"))
}
print(summary_tbl)

message("\nNOTE: 'std' fields are on MRlap's internal Z/sqrt(N)-standardised ",
        "scale, comparable to each other but NOT interpretable as log-odds ",
        "or a raw beta - do not exponentiate to an OR (see REVISION LOG).")
message("\n=== 05.5_mrlap_praevia_endoMR-PREG.R completed ===")
