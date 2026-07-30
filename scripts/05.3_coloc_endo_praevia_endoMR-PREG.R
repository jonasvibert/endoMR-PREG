#!/usr/bin/env Rscript
################################################################################
# Script: 05.3_coloc_endo_praevia_endoMR-PREG.R
# Project: endoMR-PREG
#
# REVISION LOG (revision-round1-reviewers)
#   [R2.1] Reviewer #2 ("I recommend performing colocalization analysis...
#   for the positive finding (endometriosis and placenta previa)"):
#   locus-level Bayesian colocalisation (coloc.abf, Giambartolomei et al. 2014)
#   at each of the loci contributing to the main endometriosis instrument, to
#   assess whether the endometriosis and placenta praevia GWAS signals at
#   each region are consistent with a shared causal variant vs. two distinct
#   variants in LD (confounding by local LD) vs. no signal in one/both traits.
#
# Data sources:
#   - Exposure: Rahmioglu et al. 2023 endometriosis GWAS, full genome-wide
#     summary statistics (GWAS Catalog GCST90205183, build GRCh37), obtained
#     from external drive (My Passport/GWAS/endometriosis_Rahmioglu_2023/).
#     This file has NO rsID column (chromosome/base_pair_location only), so
#     rsIDs are assigned via a CHR+BP join against the HapMap3 SNP positions
#     bundled in the LDSC eur_w_ld_chr reference (see below) - this also
#     conveniently restricts the exposure side to a well-characterised,
#     LD-representative SNP set (~1.2M HapMap3 SNPs), which is standard
#     practice for regional coloc when full imputed dosages aren't available.
#   - Outcome: FinnGen R12 placenta praevia (O15_PLAC_PRAEVIA), full
#     genome-wide, build GRCh38, already has native rsIDs
#     (data/OUTCOME_FINNGEN/finngen_R12_O15_PLAC_PRAEVIA).
#   - rsID<->CHR:BP(GRCh37) map: data/LDSC_REF/eur_w_ld_chr/*.l2.ldscore.gz
#     (CHR, SNP, BP columns). Because both exposure and outcome are merged on
#     rsID (a build-independent identifier), no GRCh37->GRCh38 liftover is
#     needed despite the build mismatch between the two source files.
#   - Loci: the 41 independent genome-wide-significant SNPs of the main
#     Rahmioglu endometriosis instrument used throughout the manuscript
#     (results/Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv), each
#     defining a +/-500kb region.
#
# Method:
#   For each locus, coloc.abf() is run on the harmonised beta/varbeta of both
#   traits within the +/-500kb window (Wakefield approximate Bayes factors;
#   default priors p1=p2=1e-4, p12=1e-5). PP.H4 = posterior probability that
#   both traits share a single causal variant in the region; PP.H3 = two
#   distinct causal variants; PP.H0-H2 = no/one-trait-only signal.
#   No case/control N or sample overlap correction is used here (coloc.abf's
#   ABF only requires beta/varbeta, not N) - this analysis addresses LD
#   confounding, not the same overlap/winner's-curse concerns that MRlap
#   (script 05.5) targets.
#
# Output:
#   - results/coloc_endo_praevia_by_locus.csv (PP.H0-H4 per locus)
#   - results/coloc_endo_praevia_summary.csv  (overview across loci)
################################################################################

### 1) SETUP ####################################################################

required_pkgs <- c("data.table", "dplyr", "coloc", "here")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

source(here::here("config", "config.R"))

results_dir <- here::here("results")
data_dir    <- here::here("data")

rahmioglu_gw_file <- RAHMIOGLU_GW_FILE
finngen_praevia_file <- file.path(data_dir, "OUTCOME_FINNGEN", "finngen_R12_O15_PLAC_PRAEVIA")
ldsc_ref_dir <- file.path(data_dir, "LDSC_REF", "eur_w_ld_chr")
instrument_file <- file.path(results_dir, "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv")

if (!file.exists(rahmioglu_gw_file)) stop("Rahmioglu genome-wide file not found: ", rahmioglu_gw_file, " - set RAHMIOGLU_GW_FILE env var or edit config/config.R")
if (!file.exists(finngen_praevia_file)) stop("Missing file: ", finngen_praevia_file)
if (!dir.exists(ldsc_ref_dir)) stop("Missing LDSC reference dir: ", ldsc_ref_dir)
if (!file.exists(instrument_file)) stop("Missing file: ", instrument_file)

WINDOW_BP <- 500000L

### 2) BUILD RSID <-> CHR:BP (GRCh37) MAP FROM LDSC REFERENCE ##################

message("Building rsID <-> CHR:BP map from eur_w_ld_chr (HapMap3 SNPs, GRCh37)...")

ldscore_files <- list.files(ldsc_ref_dir, pattern = "^[0-9]+\\.l2\\.ldscore\\.gz$", full.names = TRUE)
stopifnot(length(ldscore_files) == 22)

rsid_map <- data.table::rbindlist(lapply(ldscore_files, function(f) {
  data.table::fread(f, select = c("CHR", "SNP", "BP"))
}))
data.table::setnames(rsid_map, c("CHR", "SNP", "BP"), c("chr", "rsid", "bp"))
message("rsID map built: ", nrow(rsid_map), " HapMap3 SNPs across 22 chromosomes.")

### 3) LOAD + ANNOTATE RAHMIOGLU GENOME-WIDE SUMSTATS WITH RSID ################

message("Loading Rahmioglu genome-wide sumstats (this is a large file, may take a minute)...")

rahm <- data.table::fread(rahmioglu_gw_file)
data.table::setnames(rahm,
  c("chromosome", "base_pair_location", "effect_allele", "other_allele",
    "effect_allele_frequency", "beta", "standard_error", "p_value"),
  c("chr", "bp", "ea", "oa", "eaf", "beta", "se", "pval"))
rahm$chr <- as.character(rahm$chr)
rsid_map$chr <- as.character(rsid_map$chr)

rahm <- merge(rahm, rsid_map, by = c("chr", "bp"), all = FALSE)
message("Rahmioglu SNPs matched to a HapMap3 rsID: ", nrow(rahm))

### 4) LOAD OUTCOME (FINNGEN PLACENTA PRAEVIA, GENOME-WIDE) ####################

message("Loading FinnGen placenta praevia genome-wide sumstats...")

finngen <- data.table::fread(
  finngen_praevia_file,
  select = c("rsids", "ref", "alt", "beta", "sebeta", "af_alt", "pval")
)
data.table::setnames(finngen,
  c("rsids", "ref", "alt", "beta", "sebeta", "af_alt", "pval"),
  c("rsid", "oa", "ea", "beta_out", "se_out", "eaf_out", "pval_out"))
# FinnGen rsids field can contain multiple semicolon-separated rsIDs per row; keep the first
finngen$rsid <- sub(";.*$", "", finngen$rsid)
finngen <- finngen[finngen$rsid != "" & !is.na(finngen$rsid), ]
finngen <- finngen[!duplicated(finngen$rsid), ]
message("FinnGen placenta praevia SNPs with usable rsID: ", nrow(finngen))

### 5) MAIN INSTRUMENT LOCI (41 SNPs) ###########################################

instrument <- data.table::fread(instrument_file)
loci <- data.frame(
  SNP = instrument$SNP,
  chr = as.character(instrument$chr.exposure),
  pos = as.integer(sub(".*:", "", instrument$chrpos))
)
message("Main instrument loci for coloc: ", nrow(loci))

### 6) PER-LOCUS COLOC.ABF #######################################################

run_locus_coloc <- function(lead_snp, chr, pos) {

  exp_region <- rahm[rahm$chr == chr & rahm$bp >= pos - WINDOW_BP & rahm$bp <= pos + WINDOW_BP, ]
  exp_region <- exp_region[!duplicated(exp_region$rsid), ]

  region <- merge(
    exp_region[, c("rsid", "ea", "oa", "eaf", "beta", "se", "pval")],
    finngen[, c("rsid", "ea", "oa", "eaf_out", "beta_out", "se_out", "pval_out")],
    by = "rsid"
  )
  if (nrow(region) < 20) {
    return(data.frame(SNP = lead_snp, chr = chr, pos = pos, n_snps = nrow(region),
                       PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA,
                       note = "too few overlapping SNPs (<20), skipped"))
  }

  # merge() suffixed the colliding "ea"/"oa" columns: .x = exposure, .y = outcome
  data.table::setnames(region, c("ea.x", "oa.x", "ea.y", "oa.y"), c("ea", "oa", "ea_out", "oa_out"))

  # harmonise alleles: flip outcome beta sign where alleles are swapped
  same      <- (region$ea == region$ea_out & region$oa == region$oa_out) |
               (region$ea == region$oa_out & region$oa == region$ea_out)
  swapped   <- region$ea == region$oa_out & region$oa == region$ea_out
  region <- region[same, ]
  region$beta_out[swapped[same]] <- -region$beta_out[swapped[same]]
  # drop remaining palindromic/ambiguous SNPs (A/T, C/G with EAF near 0.5)
  palindromic <- (region$ea == "A" & region$oa == "T") | (region$ea == "T" & region$oa == "A") |
                 (region$ea == "C" & region$oa == "G") | (region$ea == "G" & region$oa == "C")
  ambiguous <- palindromic & region$eaf > 0.4 & region$eaf < 0.6
  region <- region[!ambiguous, ]

  if (nrow(region) < 20) {
    return(data.frame(SNP = lead_snp, chr = chr, pos = pos, n_snps = nrow(region),
                       PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA,
                       note = "too few overlapping SNPs after harmonisation (<20), skipped"))
  }

  d1 <- list(snp = region$rsid, beta = region$beta, varbeta = region$se^2,
             type = "cc", MAF = pmin(region$eaf, 1 - region$eaf))
  d2 <- list(snp = region$rsid, beta = region$beta_out, varbeta = region$se_out^2,
             type = "cc", MAF = pmin(region$eaf_out, 1 - region$eaf_out))

  res <- tryCatch(
    suppressWarnings(coloc::coloc.abf(dataset1 = d1, dataset2 = d2)),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(data.frame(SNP = lead_snp, chr = chr, pos = pos, n_snps = nrow(region),
                       PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA,
                       note = "coloc.abf failed"))
  }

  s <- as.list(res$summary)
  data.frame(SNP = lead_snp, chr = chr, pos = pos, n_snps = nrow(region),
             PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf,
             PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf, note = "ok")
}

message("Running per-locus coloc.abf across ", nrow(loci), " loci...")
coloc_results <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
  message("  [", i, "/", nrow(loci), "] ", loci$SNP[i])
  run_locus_coloc(loci$SNP[i], loci$chr[i], loci$pos[i])
}))

write.csv(coloc_results, file.path(results_dir, "coloc_endo_praevia_by_locus.csv"), row.names = FALSE)

### 7) SUMMARY ###################################################################

ok <- coloc_results[coloc_results$note == "ok", ]
summary_tbl <- data.frame(
  n_loci_tested        = nrow(coloc_results),
  n_loci_analysed      = nrow(ok),
  n_loci_skipped       = nrow(coloc_results) - nrow(ok),
  n_loci_PPH4_gt_0.8   = sum(ok$PP.H4 > 0.8, na.rm = TRUE),
  n_loci_PPH4_gt_0.5   = sum(ok$PP.H4 > 0.5, na.rm = TRUE),
  max_PPH4             = if (nrow(ok) > 0) max(ok$PP.H4, na.rm = TRUE) else NA,
  locus_max_PPH4        = if (nrow(ok) > 0) ok$SNP[which.max(ok$PP.H4)] else NA,
  mean_PPH4             = if (nrow(ok) > 0) mean(ok$PP.H4, na.rm = TRUE) else NA
)
write.csv(summary_tbl, file.path(results_dir, "coloc_endo_praevia_summary.csv"), row.names = FALSE)

print(summary_tbl)
message("\n=== 05.3_coloc_endo_praevia_endoMR-PREG.R completed ===")
