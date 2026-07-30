#!/usr/bin/env Rscript
################################################################################
# Script: 05.4_ldsc_genetic_correlation_endoMR-PREG.R
# Project: endoMR-PREG
#
# REVISION LOG (revision-round1-reviewers)
#   [R2.2] Reviewer #2 ("evaluate the genetic correlation between
#   endometriosis and placenta previa (and, where power permits, other
#   placental phenotypes)... using LD score regression"). MiXeR/GNOVA
#   (polygenic overlap, also requested in R2.2) declined per coauthor
#   decision (Carolina) - not run.
#   [R4.5] Reviewer #4 ("An underpowered MR null can't speak to shared
#   etiology; LDSC rg uses genome-wide signal, is better powered, and is
#   robust to sample overlap") - requested for the principal domains of null
#   outcomes: hypertensive disorders of pregnancy (HDP), gestational diabetes
#   (GDM), and postpartum haemorrhage (PPH), and specifically named
#   "placental disorders-other" (FinnGen PLAC_DISORD) and "premature
#   placental separation" (FinnGen PLAC_PREMAT_SEPAR) as example
#   low-case-count nulls needing rg context.
#
# Method: bivariate LD score regression (Bulik-Sullivan et al. 2015),
#   implemented via GenomicSEM::munge() + GenomicSEM::ldsc() (an R-native
#   LDSC reimplementation - avoids the separate Python 2 ldsc.py toolchain).
#   LD scores / regression weights: European 1000 Genomes HapMap3 LD scores
#   (data/LDSC_REF/eur_w_ld_chr), standard reference for LDSC in
#   European-ancestry GWAS. sample.prev/population.prev are left NA for all
#   traits: rg is invariant to the liability-scale transformation, only h2
#   is affected, and h2 is not the target estimand here, so this analysis
#   does not depend on external endometriosis/placenta praevia population
#   prevalence figures.
#
# IMPORTANT SCOPE LIMITATION:
#   LDSC rg requires GENOME-WIDE summary statistics for both traits. Of the
#   three null-outcome domains requested by Reviewer #4, only postpartum
#   haemorrhage (PPH, via the Westergaard et al. meta-analysis,
#   data/OUTCOME_PPH/) has genome-wide summary statistics available.
#   Hypertensive disorders of pregnancy (hdp_subsamp) and gestational
#   diabetes (gdm_subsamp) are sourced from the MR-PREG collaboration
#   (data/OUTCOME_MR-PREG/ma_out_dat.txt), which is pre-restricted at
#   source to ~43 instrument SNPs (see [[project_reviewer_response_round1]]
#   memory / script 05.1 header) - genome-wide LDSC is NOT possible for
#   these two domains with the data available to this project. This is
#   the same "not genome-wide" limitation already documented for the
#   Koller sensitivity analysis (script 05.1) and applies here for the
#   identical reason. This script therefore reports rg for:
#     (a) endometriosis <-> placenta praevia (FinnGen R12, the positive
#         finding), and
#     (b) endometriosis <-> postpartum haemorrhage (Westergaard PPH), as
#         the one null-outcome domain with genome-wide data.
#   HDP and GDM are recorded in the output as "not_feasible" with the
#   reason above, rather than silently omitted. PLAC_DISORD and
#   PLAC_PREMAT_SEPAR ARE genome-wide (FinnGen R12, same source as
#   PLAC_PRAEVIA) and so ARE estimated here. This script therefore reports
#   rg for:
#     (a) endometriosis <-> placenta praevia (FinnGen R12, the positive
#         finding),
#     (b) endometriosis <-> placental disorders-other (FinnGen PLAC_DISORD),
#     (c) endometriosis <-> premature placental separation (FinnGen
#         PLAC_PREMAT_SEPAR), and
#     (d) endometriosis <-> postpartum haemorrhage (Westergaard PPH).
#   HDP and GDM remain "not_feasible" with the reason above.
#
# Data sources (same as script 05.3):
#   - Exposure: Rahmioglu et al. 2023, full genome-wide GRCh37
#     (My Passport/GWAS/endometriosis_Rahmioglu_2023/), rsID assigned via
#     CHR:BP join against the HapMap3 SNP positions in eur_w_ld_chr.
#   - Outcomes 1-3: FinnGen R12 O15_PLAC_PRAEVIA / PLAC_DISORD /
#     PLAC_PREMAT_SEPAR, genome-wide, native rsID.
#   - Outcome 4: Westergaard PPH meta-analysis, Postpartum_hemorrhage.txt,
#     genome-wide, native rsID; reports OR (not beta) and no SE column, so
#     beta = log(OR) and SE is back-calculated from the two-sided p-value:
#     SE = |beta| / qnorm(1 - p/2).
#   - Sample sizes: fixed TOTAL N per trait (not per-SNP effective N), as
#     required by standard LDSC practice for case-control GWAS (using
#     effective N here would bias the h2/rg regression, since the LDSC
#     model chi2_j ~ N*h2*l2_j/M + ... requires N to be the actual sample
#     size). Rahmioglu: 762,600 (60,674 cases + 701,926 controls, per
#     Rahmioglu et al. 2023 / Table 1). FinnGen R12: praevia 223,001
#     (1,815/221,186), PLAC_DISORD 221,519 (333/221,186), PLAC_PREMAT_SEPAR
#     222,061 (875/221,186) - per FinnGen R12 phenotype manifest, matching
#     Table 2 (script 06). Westergaard PPH: 331,792 (Postpartum_hemorrhage;
#     Table 2). See REVISION LOG: an earlier version of this script used
#     N_eff = 1/(2*eaf*(1-eaf)*se^2) for FinnGen/PPH and an undocumented
#     RAHM_N = 205,904 for Rahmioglu - both contradicted standard LDSC/N
#     requirements and have been corrected here.
#
# Output:
#   - results/ldsc_rg_endo_outcomes.csv
################################################################################
#
# REVISION LOG (2026-07-27, N correction)
#   [FIX] Replaced per-SNP effective N (1/(2*eaf*(1-eaf)*se^2)) with fixed
#   total N per trait for all case-control outcomes (FinnGen x3, PPH), and
#   replaced the undocumented RAHM_N=205904 with the published total N for
#   Rahmioglu et al. 2023 (762,600 = 60,674 cases + 701,926 controls).
#   Rationale: standard LDSC (and the MRlap analysis in script 05.5, which
#   shares this exact bug) requires the true total GWAS sample size, not an
#   effective/adjusted N, per LDSC/MRlap documentation. This was flagged
#   during external review of the sample-overlap response (Reviewer #4)
#   and traced to affect this script too. Results re-run after the fix;
#   compare results/ldsc_rg_endo_outcomes.csv before/after if auditing.
################################################################################

### 1) SETUP ####################################################################

required_pkgs <- c("data.table", "dplyr", "GenomicSEM", "here")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

source(here::here("config", "config.R"))

results_dir <- here::here("results")
data_dir    <- here::here("data")
munge_dir   <- here::here("data", "LDSC_MUNGE")
dir.create(munge_dir, showWarnings = FALSE)

rahmioglu_gw_file      <- RAHMIOGLU_GW_FILE
finngen_praevia_file   <- file.path(data_dir, "OUTCOME_FINNGEN", "finngen_R12_O15_PLAC_PRAEVIA")
finngen_disord_file    <- file.path(data_dir, "OUTCOME_FINNGEN", "finngen_R12_O15_PLAC_DISORD")
finngen_premsep_file   <- file.path(data_dir, "OUTCOME_FINNGEN", "finngen_R12_O15_PLAC_PREMAT_SEPAR")
pph_file               <- file.path(data_dir, "OUTCOME_PPH", "Postpartum_hemorrhage.txt")
ldsc_ref_dir           <- file.path(data_dir, "LDSC_REF", "eur_w_ld_chr")
hm3_file               <- file.path(data_dir, "LDSC_REF", "w_hm3.snplist")

if (!file.exists(rahmioglu_gw_file)) stop("Rahmioglu genome-wide file not found: ", rahmioglu_gw_file, " - set RAHMIOGLU_GW_FILE env var or edit config/config.R")
if (!file.exists(finngen_praevia_file)) stop("Missing file: ", finngen_praevia_file)
if (!file.exists(finngen_disord_file)) stop("Missing file: ", finngen_disord_file)
if (!file.exists(finngen_premsep_file)) stop("Missing file: ", finngen_premsep_file)
if (!file.exists(pph_file)) stop("Missing file: ", pph_file)
if (!dir.exists(ldsc_ref_dir)) stop("Missing LDSC reference dir: ", ldsc_ref_dir)
if (!file.exists(hm3_file)) stop("Missing file: ", hm3_file)

message("NOTE: HDP (hdp_subsamp) and GDM (gdm_subsamp) domains requested by ",
        "Reviewer #4 are NOT genome-wide (MR-PREG-collaboration-restricted ",
        "outcome files, ~43 SNPs) - rg cannot be estimated for them. See ",
        "script header for detail. Proceeding with placenta praevia + PPH only.")

### 2) RSID <-> CHR:BP (GRCh37) MAP FROM LDSC REFERENCE #########################

message("Building rsID <-> CHR:BP map from eur_w_ld_chr (HapMap3 SNPs, GRCh37)...")
ldscore_files <- list.files(ldsc_ref_dir, pattern = "^[0-9]+\\.l2\\.ldscore\\.gz$", full.names = TRUE)
stopifnot(length(ldscore_files) == 22)
rsid_map <- data.table::rbindlist(lapply(ldscore_files, function(f) {
  data.table::fread(f, select = c("CHR", "SNP", "BP"))
}))
data.table::setnames(rsid_map, c("CHR", "SNP", "BP"), c("chr", "rsid", "bp"))
rsid_map$chr <- as.character(rsid_map$chr)

### 3) PREPARE MUNGE-READY SUMSTATS FILES #######################################

RAHM_N <- 762600L  # total N, Rahmioglu et al. 2023: 60,674 cases + 701,926 controls

message("Preparing exposure (Rahmioglu) sumstats for munging...")
rahm <- data.table::fread(rahmioglu_gw_file)
data.table::setnames(rahm,
  c("chromosome", "base_pair_location", "effect_allele", "other_allele",
    "effect_allele_frequency", "beta", "standard_error", "p_value"),
  c("chr", "bp", "A1", "A2", "MAF", "effect", "se", "P"))
rahm$chr <- as.character(rahm$chr)
rahm <- merge(rahm, rsid_map, by = c("chr", "bp"), all = FALSE)
rahm$N <- RAHM_N
rahm_out <- rahm[, c("rsid", "A1", "A2", "MAF", "effect", "P", "N")]
data.table::setnames(rahm_out, "rsid", "SNP")
rahm_munge_input <- file.path(munge_dir, "rahmioglu_for_munge.tsv")
data.table::fwrite(rahm_out, rahm_munge_input, sep = "\t")
message("  Rahmioglu prepared: ", nrow(rahm_out), " SNPs -> ", rahm_munge_input)

prepare_finngen <- function(path, label, out_name, total_n) {
  message("Preparing outcome (FinnGen ", label, ") sumstats for munging...")
  fg <- data.table::fread(path,
                           select = c("rsids", "ref", "alt", "beta", "sebeta", "af_alt", "pval"))
  data.table::setnames(fg, c("rsids", "ref", "alt", "beta", "sebeta", "af_alt", "pval"),
                        c("SNP", "A2", "A1", "effect", "se", "MAF", "P"))
  fg$SNP <- sub(";.*$", "", fg$SNP)
  fg <- fg[fg$SNP != "" & !is.na(fg$SNP) & !duplicated(fg$SNP), ]
  fg$N <- total_n  # fixed total N (see REVISION LOG) - not per-SNP effective N
  fg_out <- fg[, c("SNP", "A1", "A2", "MAF", "effect", "P", "N")]
  out_path <- file.path(munge_dir, out_name)
  data.table::fwrite(fg_out, out_path, sep = "\t")
  message("  ", label, " prepared: ", nrow(fg_out), " SNPs -> ", out_path)
  message("  Total N (", label, "): ", total_n)
  out_path
}

fg_munge_input      <- prepare_finngen(finngen_praevia_file, "placenta praevia", "finngen_praevia_for_munge.tsv", 223001L)
disord_munge_input  <- prepare_finngen(finngen_disord_file, "placental disorders-other", "finngen_disord_for_munge.tsv", 221519L)
premsep_munge_input <- prepare_finngen(finngen_premsep_file, "premature placental separation", "finngen_premsep_for_munge.tsv", 222061L)

message("Preparing outcome 2 (Westergaard postpartum haemorrhage) sumstats for munging...")
pph <- data.table::fread(pph_file)
data.table::setnames(pph, c("rsName", "OA", "EA", "EAfrq", "Effect", "P"),
                      c("SNP", "A2", "A1", "MAF", "OR", "P"), skip_absent = TRUE)
pph <- pph[pph$SNP != "" & !is.na(pph$SNP) & pph$SNP != "." & !duplicated(pph$SNP), ]
pph$effect <- log(pph$OR)
pph$se <- abs(pph$effect) / qnorm(1 - pph$P / 2)
pph <- pph[is.finite(pph$se) & pph$se > 0, ]
pph$N <- 331792L  # fixed total N (Westergaard Postpartum_hemorrhage; see REVISION LOG) - not per-SNP effective N
pph_out <- pph[, c("SNP", "A1", "A2", "MAF", "effect", "P", "N")]
pph_munge_input <- file.path(munge_dir, "pph_for_munge.tsv")
data.table::fwrite(pph_out, pph_munge_input, sep = "\t")
message("  PPH prepared: ", nrow(pph_out), " SNPs -> ", pph_munge_input)
message("  Total N (PPH): ", 331792L)

### 4) MUNGE #####################################################################

trait_names  <- c("endometriosis", "placenta_praevia", "PLAC_DISORD", "PLAC_PREMAT_SEPAR", "PPH")
input_files  <- c(rahm_munge_input, fg_munge_input, disord_munge_input, premsep_munge_input, pph_munge_input)

message("\nMunging all five sumstats against HapMap3...")
GenomicSEM::munge(
  files       = input_files,
  hm3         = hm3_file,
  trait.names = trait_names,
  info.filter = 0,   # no INFO column available (not imputation-quality filtered upstream)
  maf.filter  = 0.01
)
# munge() writes <trait.names[i]>.sumstats.gz to the current working directory
munged_files <- paste0(trait_names, ".sumstats.gz")
stopifnot(all(file.exists(munged_files)))
file.copy(munged_files, munge_dir, overwrite = TRUE)
file.remove(munged_files)
munged_paths <- file.path(munge_dir, munged_files)

### 5) LDSC BIVARIATE GENETIC CORRELATION (PAIRWISE) #############################

# Run endometriosis vs each outcome as a SEPARATE 2-trait ldsc() call, rather
# than one joint 5-trait call: a single trait with a noisy/negative h2
# estimate (expected for a rare, heterogeneous "other" endpoint - see
# PLAC_DISORD below) makes GenomicSEM::ldsc() skip S_Stand/V_Stand for the
# WHOLE matrix (all(diag(S)>0) check), silently losing every pair, not just
# the problematic one. Pairwise calls isolate that failure to one outcome.

outcome_labels <- c("placenta_praevia", "PLAC_DISORD", "PLAC_PREMAT_SEPAR", "PPH")
endo_path <- munged_paths[1]
outcome_paths <- munged_paths[2:5]

run_pairwise_rg <- function(outcome_path, outcome_label) {
  out <- tryCatch(
    GenomicSEM::ldsc(
      traits          = c(endo_path, outcome_path),
      sample.prev     = c(NA, NA),
      population.prev = c(NA, NA),
      ld              = ldsc_ref_dir,
      wld             = ldsc_ref_dir,
      trait.names     = c("endometriosis", outcome_label),
      stand           = TRUE
    ),
    error = function(e) e
  )
  if (inherits(out, "error") || is.null(out$S_Stand)) {
    return(list(rg = NA, se = NA, p = NA, h2_endo = NA, h2_out = NA,
                status = "not_estimable",
                note = "LDSC could not standardise (negative h2 estimate for one trait - underpowered/heterogeneous phenotype)"))
  }
  SE_mat <- matrix(0, 2, 2)
  SE_mat[lower.tri(SE_mat, diag = TRUE)] <- sqrt(diag(out$V_Stand))
  rg <- out$S_Stand[2, 1]
  se <- SE_mat[2, 1]
  list(rg = rg, se = se, p = 2 * pnorm(-abs(rg / se)),
       h2_endo = out$S[1, 1], h2_out = out$S[2, 2],
       status = "estimated", note = "")
}

pairwise_res <- lapply(seq_along(outcome_paths), function(i) {
  message("  LDSC: endometriosis vs ", outcome_labels[i], "...")
  run_pairwise_rg(outcome_paths[i], outcome_labels[i])
})
names(pairwise_res) <- outcome_labels

results_tbl <- data.frame(
  pair        = c("endometriosis - placenta_praevia (FinnGen)",
                   "endometriosis - placental_disorders_other (FinnGen PLAC_DISORD)",
                   "endometriosis - premature_placental_separation (FinnGen PLAC_PREMAT_SEPAR)",
                   "endometriosis - PPH (Westergaard)",
                   "endometriosis - HDP (hdp_subsamp)",
                   "endometriosis - GDM (gdm_subsamp)"),
  rg          = c(sapply(pairwise_res, `[[`, "rg"), NA, NA),
  rg_se       = c(sapply(pairwise_res, `[[`, "se"), NA, NA),
  rg_pval     = c(sapply(pairwise_res, `[[`, "p"), NA, NA),
  h2_trait1   = c(sapply(pairwise_res, `[[`, "h2_endo"), NA, NA),
  h2_trait2   = c(sapply(pairwise_res, `[[`, "h2_out"), NA, NA),
  status      = c(sapply(pairwise_res, `[[`, "status"), "not_feasible", "not_feasible"),
  note        = c(sapply(pairwise_res, `[[`, "note"),
                   "outcome not genome-wide (MR-PREG-collaboration-restricted, ~43 SNPs)",
                   "outcome not genome-wide (MR-PREG-collaboration-restricted, ~43 SNPs)")
)

write.csv(results_tbl, file.path(results_dir, "ldsc_rg_endo_outcomes.csv"), row.names = FALSE)
print(results_tbl)

message("\n=== 05.4_ldsc_genetic_correlation_endoMR-PREG.R completed ===")
