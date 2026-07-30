#!/usr/bin/env Rscript
################################################################################
# Script: 05.6_ldsc_targeted_endo_praevia_endoMR-PREG.R
# Project: endoMR-PREG
#
# [R2.2]/[R4.5] Dedicated endometriosis <-> placenta praevia cross-trait LDSC,
# superseding the relevant pair from script 05.4. Two fixes vs. 05.4/05.5:
#
#   1) N: 05.4/05.5 used N=762,600 for the Rahmioglu exposure. That's GWAS
#      Catalog accession GCST90258638, not GCST90205183 (the file actually
#      read here, N=475,160). Fixed to 475,160. rg is invariant to a uniform
#      per-trait N (cancels in the h2/covariance ratio), so this changes
#      h2(endometriosis) but not rg materially. The 762,600 figure is still
#      correct for the primary MR instrument (different data source, Supp
#      Table 33 Top10K-SNPs, includes 23andMe) - not a project-wide error.
#   2) QC: GenomicSEM::munge() doesn't drop strand-ambiguous SNPs or exclude
#      the MHC region by default (checked against source). Added here:
#      ambiguous-SNP removal, MHC exclusion (chr6:25-35Mb, both builds),
#      indel/multiallelic removal, duplicate-rsID removal.
#
# Three variants, run for comparison: "original" (05.4 as-is, wrong N, no
# QC filter, reused not rerun), "n_corrected" (N fixed, no QC filter),
# "primary" (N fixed + QC filter - use this one).
#
# Method: GenomicSEM::munge()/ldsc() (Bulik-Sullivan et al. 2015 LDSC,
# R implementation) - no Python ldsc.py available on this machine.
# European 1000G HapMap3 LD scores (data/LDSC_REF/eur_w_ld_chr).
# sample.prev/population.prev left NA (rg unaffected; no placenta praevia
# population prevalence assumed). Observed-scale h2 reported throughout.
#
# Data: exposure = Rahmioglu 2023 GCST90205183, GRCh37, no native rsID
# (assigned via CHR:BP join to HapMap3). Outcome = FinnGen R12
# O15_PLAC_PRAEVIA, GRCh38, native rsID, N=223,001 (matches manifest).
# Merged on rsID, sidesteps the build mismatch. See config/config.R for
# the exposure file path.
#
# Output: results/ldsc_endometriosis_placenta_praevia_20260728/
################################################################################

### 1) SETUP ####################################################################

required_pkgs <- c("data.table", "dplyr", "GenomicSEM", "here")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

source(here::here("config", "config.R"))

data_dir    <- here::here("data")
munge_dir   <- here::here("data", "LDSC_MUNGE")
# [SYNC] This dated folder name is also hardcoded in scripts/06_tables_endoMR-PREG.R
# (variable `ldsc_dir`, Supplementary Table S10 section). If this script is
# ever re-run on a different date (new dated folder), update BOTH places -
# script 06 will otherwise keep silently reading this stale 2026-07-28 folder
# without any error, since it only checks file.exists(), not recency.
out_dir     <- here::here("results", "ldsc_endometriosis_placenta_praevia_20260728")
munge_log_dir <- file.path(out_dir, "munge_logs")
ldsc_log_dir  <- file.path(out_dir, "ldsc_logs")
dir.create(munge_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(munge_log_dir, showWarnings = FALSE)
dir.create(ldsc_log_dir, showWarnings = FALSE)

rahmioglu_gw_file    <- RAHMIOGLU_GW_FILE
finngen_praevia_file <- file.path(data_dir, "OUTCOME_FINNGEN", "finngen_R12_O15_PLAC_PRAEVIA")
ldsc_ref_dir         <- file.path(data_dir, "LDSC_REF", "eur_w_ld_chr")
hm3_file             <- file.path(data_dir, "LDSC_REF", "w_hm3.snplist")

if (!file.exists(rahmioglu_gw_file)) stop("Rahmioglu genome-wide file not found: ", rahmioglu_gw_file, " - set RAHMIOGLU_GW_FILE env var or edit config/config.R")
if (!file.exists(finngen_praevia_file)) stop("Missing file: ", finngen_praevia_file)
if (!dir.exists(ldsc_ref_dir)) stop("Missing LDSC reference dir: ", ldsc_ref_dir)
if (!file.exists(hm3_file)) stop("Missing file: ", hm3_file)

# NOTE: the "original" variant (N=762,600, the wrong figure - see FIX 1 above)
# is NOT rebuilt by this script; it is reused as-is from the existing
# endometriosis.sumstats.gz / placenta_praevia.sumstats.gz written by script
# 05.4 (see munged_orig below). No N constant for it is needed here.
RAHM_N_CORR    <- 475160L   # CORRECT total N for GCST90205183 (the file actually used), confirmed via GWAS Catalog API 2026-07-28
PRAEV_N        <- 223001L   # confirmed against finngen_R12_manifest.csv (1,815 cases / 221,186 controls)
MHC_CHR        <- "6"
MHC_START      <- 25000000L
MHC_END        <- 35000000L
AMBIGUOUS_PAIRS <- c("AT", "TA", "CG", "GC")

qc_counts <- list()  # collects filtering counts for qc_filtering_counts.csv

### 2) RSID <-> CHR:BP (GRCh37, HapMap3) MAP FROM LDSC REFERENCE ###############

message("Building rsID <-> CHR:BP map from eur_w_ld_chr (HapMap3 SNPs, GRCh37)...")
ldscore_files <- list.files(ldsc_ref_dir, pattern = "^[0-9]+\\.l2\\.ldscore\\.gz$", full.names = TRUE)
stopifnot(length(ldscore_files) == 22)
rsid_map <- data.table::rbindlist(lapply(ldscore_files, function(f) {
  data.table::fread(f, select = c("CHR", "SNP", "BP"))
}))
data.table::setnames(rsid_map, c("CHR", "SNP", "BP"), c("chr", "rsid", "bp"))
rsid_map$chr <- as.character(rsid_map$chr)
rsid_map <- rsid_map[!duplicated(rsid_map, by = c("chr", "bp")), ]  # multi-mapping guard

### 3) PREPARE EXPOSURE (RAHMIOGLU) - BASE + QC FLAGS ##########################

message("Loading Rahmioglu genome-wide sumstats...")
rahm <- data.table::fread(rahmioglu_gw_file)
data.table::setnames(rahm,
  c("chromosome", "base_pair_location", "effect_allele", "other_allele",
    "effect_allele_frequency", "beta", "standard_error", "p_value"),
  c("chr", "bp", "A1", "A2", "MAF", "effect", "se", "P"))
rahm$chr <- as.character(rahm$chr)
n0 <- nrow(rahm)
rahm <- merge(rahm, rsid_map, by = c("chr", "bp"), all = FALSE)
qc_counts[["rahm_no_rsid_in_hm3"]] <- n0 - nrow(rahm)

rahm$A1 <- toupper(rahm$A1); rahm$A2 <- toupper(rahm$A2)
n1 <- nrow(rahm)
is_acgt <- function(x) x %in% c("A", "C", "G", "T")
rahm <- rahm[is_acgt(A1) & is_acgt(A2), ]
qc_counts[["rahm_indel_multiallelic"]] <- n1 - nrow(rahm)

n2 <- nrow(rahm)
rahm <- rahm[!duplicated(rahm$rsid), ]
qc_counts[["rahm_duplicate_rsid"]] <- n2 - nrow(rahm)

rahm$pair <- paste0(rahm$A1, rahm$A2)
rahm$is_ambiguous <- rahm$pair %in% AMBIGUOUS_PAIRS
rahm$is_mhc <- (rahm$chr == MHC_CHR & rahm$bp >= MHC_START & rahm$bp <= MHC_END)
qc_counts[["rahm_ambiguous_flagged"]] <- sum(rahm$is_ambiguous)
qc_counts[["rahm_mhc_flagged"]] <- sum(rahm$is_mhc)

message("  Rahmioglu after rsID/indel/dup QC: ", nrow(rahm), " SNPs (",
        sum(rahm$is_ambiguous), " ambiguous, ", sum(rahm$is_mhc), " in extended MHC)")

### 4) PREPARE OUTCOME (FINNGEN PRAEVIA) - BASE + QC FLAGS #####################

message("Loading FinnGen placenta praevia genome-wide sumstats...")
fg <- data.table::fread(finngen_praevia_file,
                         select = c("#chrom", "pos", "rsids", "ref", "alt", "beta", "sebeta", "af_alt", "pval"))
data.table::setnames(fg, c("#chrom", "pos", "rsids", "ref", "alt", "beta", "sebeta", "af_alt", "pval"),
                      c("chr", "bp", "rsid", "A2", "A1", "effect", "se", "MAF", "P"))
fg$chr <- as.character(fg$chr)
fg$rsid <- sub(";.*$", "", fg$rsid)
n0 <- nrow(fg)
fg <- fg[fg$rsid != "" & !is.na(fg$rsid), ]
qc_counts[["finngen_missing_rsid"]] <- n0 - nrow(fg)

fg$A1 <- toupper(fg$A1); fg$A2 <- toupper(fg$A2)
n1 <- nrow(fg)
fg <- fg[is_acgt(A1) & is_acgt(A2), ]
qc_counts[["finngen_indel_multiallelic"]] <- n1 - nrow(fg)

n2 <- nrow(fg)
fg <- fg[!duplicated(fg$rsid), ]
qc_counts[["finngen_duplicate_rsid"]] <- n2 - nrow(fg)

fg$pair <- paste0(fg$A1, fg$A2)
fg$is_ambiguous <- fg$pair %in% AMBIGUOUS_PAIRS
fg$is_mhc <- (fg$chr == MHC_CHR & fg$bp >= MHC_START & fg$bp <= MHC_END)
qc_counts[["finngen_ambiguous_flagged"]] <- sum(fg$is_ambiguous)
qc_counts[["finngen_mhc_flagged"]] <- sum(fg$is_mhc)

message("  FinnGen praevia after rsID/indel/dup QC: ", nrow(fg), " SNPs (",
        sum(fg$is_ambiguous), " ambiguous, ", sum(fg$is_mhc), " in extended MHC)")

write.csv(data.frame(step = names(qc_counts), n_removed_or_flagged = unlist(qc_counts)),
          file.path(out_dir, "qc_filtering_counts.csv"), row.names = FALSE)

### 5) BUILD MUNGE-READY INPUT FILES FOR THE THREE NEW VARIANTS #################
# ("original" reuses the existing endometriosis.sumstats.gz / placenta_praevia.sumstats.gz
#  from script 05.4 as-is - not rebuilt here.)

to_munge_input <- function(df, N) {
  out <- df[, c("rsid", "A1", "A2", "MAF", "effect", "P")]
  data.table::setnames(out, "rsid", "SNP")
  out$N <- N
  out
}

rahm_ncorr        <- to_munge_input(rahm, RAHM_N_CORR)                                   # N fix only
rahm_primary      <- to_munge_input(rahm[!is_ambiguous & !is_mhc, ], RAHM_N_CORR)         # N fix + QC filter
fg_primary        <- to_munge_input(fg[!is_ambiguous & !is_mhc, ], PRAEV_N)               # QC filter only (N already correct)

f_rahm_ncorr   <- file.path(munge_dir, "rahmioglu_ncorr_for_munge.tsv")
f_rahm_primary <- file.path(munge_dir, "rahmioglu_primary_for_munge.tsv")
f_fg_primary   <- file.path(munge_dir, "finngen_praevia_primary_for_munge.tsv")
data.table::fwrite(rahm_ncorr, f_rahm_ncorr, sep = "\t")
data.table::fwrite(rahm_primary, f_rahm_primary, sep = "\t")
data.table::fwrite(fg_primary, f_fg_primary, sep = "\t")

message("Prepared: ", nrow(rahm_ncorr), " SNPs (endo, N-corrected) / ",
        nrow(rahm_primary), " SNPs (endo, N-corrected + QC) / ",
        nrow(fg_primary), " SNPs (praevia, QC-filtered)")

### 6) MUNGE THE THREE NEW FILES #################################################

old_wd <- getwd()
setwd(munge_log_dir)

trait_names_new <- c("endometriosis_ncorr", "endometriosis_primary", "placenta_praevia_primary")
input_files_new <- c(f_rahm_ncorr, f_rahm_primary, f_fg_primary)

message("\nMunging 3 new variant files against HapMap3...")
GenomicSEM::munge(
  files       = input_files_new,
  hm3         = hm3_file,
  trait.names = trait_names_new,
  info.filter = 0,
  maf.filter  = 0.01
)
munged_new <- paste0(trait_names_new, ".sumstats.gz")
stopifnot(all(file.exists(munged_new)))
file.copy(munged_new, munge_dir, overwrite = TRUE)
file.remove(munged_new)
setwd(old_wd)

munged_paths_new <- setNames(file.path(munge_dir, munged_new), trait_names_new)
# reuse existing (unfiltered, N=762,600 for endo) munged files from script 05.4
munged_orig <- c(
  endometriosis_orig    = file.path(munge_dir, "endometriosis.sumstats.gz"),
  placenta_praevia_orig = file.path(munge_dir, "placenta_praevia.sumstats.gz")
)
stopifnot(all(file.exists(munged_orig)))

### 7) BIVARIATE LDSC FOR EACH VARIANT PAIRING ###################################

run_pairwise_rg <- function(exp_path, out_path, exp_label, out_label) {
  setwd(ldsc_log_dir)
  res <- tryCatch(
    GenomicSEM::ldsc(
      traits          = c(exp_path, out_path),
      sample.prev     = c(NA, NA),
      population.prev = c(NA, NA),
      ld              = ldsc_ref_dir,
      wld             = ldsc_ref_dir,
      trait.names     = c(exp_label, out_label),
      stand           = TRUE
    ),
    error = function(e) e
  )
  setwd(old_wd)
  if (inherits(res, "error") || is.null(res$S_Stand)) {
    return(list(rg = NA, rg_se = NA, rg_p = NA,
                h2_exp = NA, h2_out = NA, h2_exp_se = NA, h2_out_se = NA,
                crosstrait_intercept = NA, crosstrait_intercept_se = NA,
                status = "not_estimable",
                note = if (inherits(res, "error")) conditionMessage(res) else "LDSC could not standardise (negative h2)"))
  }
  SE_mat <- matrix(0, 2, 2)
  SE_mat[lower.tri(SE_mat, diag = TRUE)] <- sqrt(diag(res$V_Stand))
  rg    <- res$S_Stand[2, 1]
  rg_se <- SE_mat[2, 1]
  list(rg = rg, rg_se = rg_se, rg_p = 2 * pnorm(-abs(rg / rg_se)),
       h2_exp = res$S[1, 1], h2_out = res$S[2, 2],
       h2_exp_se = sqrt(res$V[1, 1]), h2_out_se = sqrt(res$V[3, 3]),
       crosstrait_intercept = NA, crosstrait_intercept_se = NA,  # parsed from log text below
       status = "estimated", note = "")
}

message("\n--- Variant b: n_corrected (endo N=475,160; no ambiguous/MHC filter) ---")
res_ncorr <- run_pairwise_rg(munged_paths_new["endometriosis_ncorr"], munged_orig["placenta_praevia_orig"],
                              "endometriosis_ncorr", "placenta_praevia_orig")

message("\n--- Variant c: primary (endo N=475,160 + ambiguous/MHC/indel/dup filter; praevia same filter) ---")
res_primary <- run_pairwise_rg(munged_paths_new["endometriosis_primary"], munged_paths_new["placenta_praevia_primary"],
                                "endometriosis_primary", "placenta_praevia_primary")

# Catch standardisation failure here (from the R object) rather than let it
# surface as a silent NA later when parsing the log text.
if (res_ncorr$status != "estimated") {
  warning("Variant 'n_corrected': LDSC could not standardise (", res_ncorr$note,
          "). Downstream log-parsed values for this variant will likely be NA.")
}
if (res_primary$status != "estimated") {
  warning("Variant 'primary': LDSC could not standardise (", res_primary$note,
          "). Downstream log-parsed values for this variant will likely be NA.")
}

### 8) PARSE LDSC LOGS FOR UNIVARIATE + CROSS-TRAIT STATS #######################

# "Chi^2" and "h2" contain a literal digit, so a naive "grab every number in
# the line" regex picks up the wrong one - extract only after the colon.
num_after_colon <- function(line) {
  if (is.na(line) || length(line) == 0) return(NA_real_)
  after <- sub("^.*:", "", line)
  as.numeric(regmatches(after, gregexpr("-?[0-9]+\\.?[0-9]*[eE]?[+-]?[0-9]*", after))[[1]])
}
num_before_word <- function(line, word) {
  if (is.na(line) || length(line) == 0) return(NA_real_)
  m <- regmatches(line, regexpr(paste0("[0-9]+(?=\\s*", word, ")"), line, perl = TRUE))
  if (length(m) == 0 || m == "") return(NA_real_)
  as.numeric(m)
}

parse_ldsc_log <- function(log_path) {
  if (!file.exists(log_path)) return(NULL)
  txt <- readLines(log_path, warn = FALSE)

  herit_idx <- grep("^Heritability Results for trait:", txt)
  stopifnot(length(herit_idx) == 2)

  parse_herit_block <- function(start_idx) {
    block <- txt[start_idx:(start_idx + 6)]
    list(
      mean_chi2 = num_after_colon(grep("Mean Chi\\^2", block, value = TRUE)[1])[1],
      lambda_gc = num_after_colon(grep("Lambda GC", block, value = TRUE)[1])[1],
      intercept = num_after_colon(grep("^Intercept:", block, value = TRUE)[1])[1],
      intercept_se = num_after_colon(grep("^Intercept:", block, value = TRUE)[1])[2],
      ratio = num_after_colon(grep("^Ratio:", block, value = TRUE)[1])[1],
      ratio_se = num_after_colon(grep("^Ratio:", block, value = TRUE)[1])[2],
      h2_obs = num_after_colon(grep("Total Observed Scale h2:", block, value = TRUE)[1])[1],
      h2_obs_se = num_after_colon(grep("Total Observed Scale h2:", block, value = TRUE)[1])[2],
      h2_z = num_after_colon(grep("^h2 Z:", block, value = TRUE)[1])[1]
    )
  }
  t1 <- parse_herit_block(herit_idx[1])
  t2 <- parse_herit_block(herit_idx[2])

  chi2_removal_lines <- grep("^Removing .* SNPs with Chi\\^2", txt, value = TRUE)
  merge_lines <- grep("remain after merging with LD-score files", txt, value = TRUE)
  n_snps_t1_final <- if (length(chi2_removal_lines) >= 1) num_before_word(chi2_removal_lines[1], "remain") else num_before_word(merge_lines[1], "remain")
  n_snps_t2_final <- if (length(chi2_removal_lines) >= 2) num_before_word(chi2_removal_lines[2], "remain") else num_before_word(merge_lines[2], "remain")

  covline  <- grep("^Cross trait Intercept:", txt, value = TRUE)[1]
  gcovline <- grep("^Total Observed Scale Genetic Covariance", txt, value = TRUE)[1]
  gcovZ    <- grep("^g_cov Z:", txt, value = TRUE)[1]
  gcovP    <- grep("^g_cov P-value:", txt, value = TRUE)[1]
  rgline   <- grep("^Genetic Correlation between", txt, value = TRUE)[1]
  n_overlap_line <- grep("SNPs remain after merging .* summary statistics$", txt, value = TRUE)[1]

  list(
    n_snps_trait1 = n_snps_t1_final, n_snps_trait2 = n_snps_t2_final,
    n_snps_overlap = num_before_word(n_overlap_line, "SNPs"),
    trait1 = t1, trait2 = t2,
    crosstrait_intercept = num_after_colon(covline)[1], crosstrait_intercept_se = num_after_colon(covline)[2],
    g_cov = num_after_colon(gcovline)[1], g_cov_se = num_after_colon(gcovline)[2],
    g_cov_z = num_after_colon(gcovZ)[1], g_cov_p = num_after_colon(gcovP)[1],
    rg = num_after_colon(rgline)[1], rg_se = num_after_colon(rgline)[2],
    warnings = grep("[Ww]arning|WARNING", txt, value = TRUE)
  )
}

# Locate the auto-generated ldsc log files (named from input basenames by GenomicSEM)
log_ncorr_path   <- file.path(ldsc_log_dir, paste0(basename(munged_paths_new["endometriosis_ncorr"]), "_",
                                                    basename(munged_orig["placenta_praevia_orig"]), "_ldsc.log"))
log_primary_path <- file.path(ldsc_log_dir, paste0(basename(munged_paths_new["endometriosis_primary"]), "_",
                                                    basename(munged_paths_new["placenta_praevia_primary"]), "_ldsc.log"))
log_orig_path    <- file.path(ldsc_log_dir, "original_endo_vs_praevia_ldsc.log")  # copied in section 9 below, parsed after copy

message("ldsc log (n_corrected variant): ", log_ncorr_path, " exists=", file.exists(log_ncorr_path))
message("ldsc log (primary variant): ", log_primary_path, " exists=", file.exists(log_primary_path))

### 9) ASSEMBLE FINAL SUMMARY TABLES #############################################

# --- copy the "original" (2026-07-27, script 05.4) log into this package's
#     log folder so all three variants are parsed identically from log text ---
orig_log <- here::here("endometriosis.sumstats.gz_placenta_praevia.sumstats.gz_ldsc.log")
if (file.exists(orig_log)) file.copy(orig_log, log_orig_path, overwrite = TRUE)
orig_munge_log <- here::here("endometriosis_placenta_praevia_PLAC_DISORD_PLAC_PREMAT_SEPAR_PPH_munge.log")
if (file.exists(orig_munge_log)) file.copy(orig_munge_log, file.path(munge_log_dir, "original_5trait_munge.log"), overwrite = TRUE)

p_orig    <- parse_ldsc_log(log_orig_path)
p_ncorr   <- parse_ldsc_log(log_ncorr_path)
p_primary <- parse_ldsc_log(log_primary_path)

variant_labels <- c("original (05.4, N=762,600 - superseded, no ambig/MHC filter)",
                     "n_corrected (N=475,160, no ambig/MHC filter)",
                     "primary (N=475,160 + ambig/MHC/indel/dup filter)")
parsed_list <- list(p_orig, p_ncorr, p_primary)

univariate_heritability <- do.call(rbind, lapply(seq_along(parsed_list), function(i) {
  p <- parsed_list[[i]]
  if (is.null(p)) {
    return(data.frame(variant = variant_labels[i], trait = c("endometriosis", "placenta_praevia"),
                       n_snps = NA, mean_chi2 = NA, lambda_gc = NA, intercept = NA, intercept_se = NA,
                       ratio = NA, ratio_se = NA, h2_obs = NA, h2_obs_se = NA, h2_z = NA,
                       h2_p = NA, warnings = "log not found/parseable"))
  }
  rbind(
    data.frame(variant = variant_labels[i], trait = "endometriosis", n_snps = p$n_snps_trait1,
               mean_chi2 = p$trait1$mean_chi2, lambda_gc = p$trait1$lambda_gc,
               intercept = p$trait1$intercept, intercept_se = p$trait1$intercept_se,
               ratio = p$trait1$ratio, ratio_se = p$trait1$ratio_se,
               h2_obs = p$trait1$h2_obs, h2_obs_se = p$trait1$h2_obs_se, h2_z = p$trait1$h2_z,
               h2_p = 2 * pnorm(-abs(p$trait1$h2_z)),
               warnings = paste(p$warnings, collapse = "; ")),
    data.frame(variant = variant_labels[i], trait = "placenta_praevia", n_snps = p$n_snps_trait2,
               mean_chi2 = p$trait2$mean_chi2, lambda_gc = p$trait2$lambda_gc,
               intercept = p$trait2$intercept, intercept_se = p$trait2$intercept_se,
               ratio = p$trait2$ratio, ratio_se = p$trait2$ratio_se,
               h2_obs = p$trait2$h2_obs, h2_obs_se = p$trait2$h2_obs_se, h2_z = p$trait2$h2_z,
               h2_p = 2 * pnorm(-abs(p$trait2$h2_z)),
               warnings = paste(p$warnings, collapse = "; "))
  )
}))
write.csv(univariate_heritability, file.path(out_dir, "univariate_heritability.csv"), row.names = FALSE)
print(univariate_heritability)

bivariate_rg <- do.call(rbind, lapply(seq_along(parsed_list), function(i) {
  p <- parsed_list[[i]]
  if (is.null(p)) {
    return(data.frame(variant = variant_labels[i], n_snps_overlap = NA, rg = NA, rg_se = NA, rg_p = NA,
                       crosstrait_intercept = NA, crosstrait_intercept_se = NA,
                       g_cov = NA, g_cov_se = NA, g_cov_z = NA, g_cov_p = NA,
                       status = "log not found/parseable"))
  }
  data.frame(variant = variant_labels[i], n_snps_overlap = p$n_snps_overlap,
             rg = p$rg, rg_se = p$rg_se, rg_p = 2 * pnorm(-abs(p$rg / p$rg_se)),
             crosstrait_intercept = p$crosstrait_intercept, crosstrait_intercept_se = p$crosstrait_intercept_se,
             g_cov = p$g_cov, g_cov_se = p$g_cov_se, g_cov_z = p$g_cov_z, g_cov_p = p$g_cov_p,
             status = "estimated")
}))
write.csv(bivariate_rg, file.path(out_dir, "bivariate_rg.csv"), row.names = FALSE)
print(bivariate_rg)

message("\n=== 05.6_ldsc_targeted_endo_praevia_endoMR-PREG.R completed ===")
message("Outputs written to: ", out_dir)
