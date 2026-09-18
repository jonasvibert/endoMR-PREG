#!/usr/bin/env Rscript
################################################################################
# Script: 05.1_koller_sensitivity_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   [R3.3] Reviewer #3 ("Are the authors able to comment on any similarities
#   or discrepancies between the top findings of this study and the one used
#   in their own analysis?"): re-derive an independent endometriosis
#   instrument from Koller et al. 2026 (EUR summary statistics) and re-run
#   the primary MR analysis + standard sensitivity analyses, then compare
#   against the main (Rahmioglu-based) analysis.
#
#   [R2.5] Reviewer #2, comment 5 (adenomyosis misclassification): attempt an
#   analogous MR of adenomyosis genetic liability -> placenta praevia using
#   the Koller et al. 2026 adenomyosis GWAS, to inform (not replace) the
#   discussion of possible endometriosis/adenomyosis misclassification.
#
# REVISION LOG
#   [R3.3 — MR-PREG RE-EXTRACTION UPDATE, 2026-09-16] Section 4b now reads
#   data/EXPOSURE_KOLLER_TABLES4/Updated_Koller_2026_SNP_list/, re-extracted
#   against the 86-SNP Table 4 instrument.
#
#   [R3.3 — INSTRUMENT UPDATE, 2026-09-03] Switched the ENDOMETRIOSIS
#   instrument (adenomyosis instrument below is UNCHANGED) from a
#   self-clumped set derived from the METAL-style genome-wide file
#   (SE = 1/sqrt(2pq N), Beta = Z*SE, N treated as effective N) to Koller et
#   al.'s own published Supplementary Table 4: 86 genome-wide significant
#   loci (p<5e-8), of which the paper's own cross-ancestry combined-endo
#   analysis reports n=80 as LD-independent (r2<0.1, 1000G EUR panel); the
#   remaining 6 rows are loci significant only in the EUR clinical/
#   self-reported sub-analyses, included in the same table (Type column:
#   "Clinical/Self-reported endometriosis, novel/previously reported").
#   Uses the table's real, directly-reported EUR "Endometriosis combined
#   definition" Beta/SE/P/EAF columns instead of back-calculating from Z+N —
#   more accurate and directly citable against the source table (all 86
#   rows have a non-NA EUR combined-endometriosis effect, checked
#   2026-09-03). No re-clumping performed: the paper's own loci are already
#   LD-independent, so TwoSampleMR::clump_data() / PLINK are no longer used
#   for this instrument. Genome build for the table's Chromosome/Position
#   columns is NOT stated anywhere in the paper's Methods, Data
#   Availability statement, or the associated Zenodo deposit README
#   (10.5281/zenodo.18983492) — checked directly, 2026-09-03 — so merging
#   against outcome data still uses rsID only (as already done for every
#   other outcome source in this script); Chromosome/Position/nearest gene
#   are carried through for reference/QC only, never as a join key. Source:
#   Koller D et al., Nat Genet 2026;58(5):1051–1061, Supplementary Table 4;
#   local copy data/EXPOSURE_KOLLER/Koller2026_SupplementaryTable4_significant_loci.xlsx.
#
#   [R3.3/R2.5 — SCOPE] The 3 outcome sources behind the 30 main outcomes do
#   NOT have equal SNP coverage:
#     - MR-PREG (data/OUTCOME_MR-PREG/ma_out_dat.txt) contains only 43 SNPs
#       total — it was provided by the consortium already restricted to the
#       original Rahmioglu instrument list, and cannot be re-queried for a
#       different SNP set within the revision timeframe.
#     - FinnGen R12 (data/OUTCOME_FINNGEN/*) and Westergaard PPH
#       (data/OUTCOME_PPH/*.txt) are full genome-wide summary statistics
#       (~18-20M variants each) and CAN be re-queried at the Koller
#       instrument's SNPs.
#   An initial attempt to harmonise the Koller instrument against the
#   pre-restricted results/formatted_all_pregnancy_outcomes.tsv (which only
#   carries the ~41-50 Rahmioglu-instrument SNPs) produced a spurious,
#   severely underpowered placenta praevia estimate (OR = 3.50, 95% CI
#   1.31-9.35, based on only ~4-5 of the 54 clumped Koller SNPs — 8/54
#   overlapped the restricted SNP list at all). That result is an artefact
#   of SNP-list mismatch, not a valid sensitivity estimate, and is NOT used.
#   Decision (confirmed with author): restrict the Koller sensitivity
#   analysis to the 7 outcomes sourced from FinnGen R12 + Westergaard PPH
#   (which includes placenta praevia, the outcome of primary interest to
#   Reviewer #3), re-extracting outcome associations genome-wide at the
#   Koller SNPs. The remaining 23 MR-PREG-sourced outcomes are explicitly
#   NOT tested here for lack of adequate data and are reported as such in
#   the response letter. See project memory
#   "project_reviewer_response_round1".
#
# Data source:
#   ENDOMETRIOSIS instrument (as of the 2026-09-03 update, see REVISION LOG):
#     data/EXPOSURE_KOLLER/Koller2026_SupplementaryTable4_significant_loci.xlsx
#     Sheet "Supplementary Table 4", data from row 6. Columns used: Lead SNP,
#     Chromosome, Position, Effect allele, Other allele, Nearest gene, Type,
#     and the EUR "Endometriosis combined definition" Beta/SE/P-value/EAF
#     block. Real published effect estimates — no derivation needed.
#   ADENOMYOSIS instrument (unchanged): METAL-style GWAS output, no beta/SE
#   columns — data/EXPOSURE_KOLLER/Koller2026_adenomyosis_EUR.tsv.gz
#   Columns: SNP, Allele1, Allele2, Freq1, N, Z, P.value, Direction, HetPVal
#   Beta/SE derived as:
#     SE   = 1 / sqrt(2 * Freq1 * (1 - Freq1) * N)
#     Beta = Z * SE
#   using the N column directly as effective sample size (author-confirmed
#   choice).
#
# Inputs:
#   - data/EXPOSURE_KOLLER/Koller2026_SupplementaryTable4_significant_loci.xlsx
#   - data/EXPOSURE_KOLLER/Koller2026_adenomyosis_EUR.tsv.gz
#   - data/OUTCOME_FINNGEN/finngen_R12_O15_PLAC_{PRAEVIA,DISORD,PREMAT_SEPAR}
#   - data/OUTCOME_PPH/{Antepartum_bleeding,Postpartum_hemorrhage,
#     Postpartum_hemorrhage_due_to_atony,
#     Postpartum_hemorrhage_due_to_retained_placenta}.txt
#   - data/EXPOSURE_KOLLER_TABLES4/Updated_Koller_2026_SNP_list/
#     tmp_Metanalysis_metanalyses-R3.mum.*.txt.gz (23 MoBa outcomes; see SCOPE
#     UPDATE below)
#   - results/ivw_results.csv (script 04; used for the Rahmioglu comparison)
#   - plink_mac_20241022/ (PLINK binary + 1000G EUR reference, script 01)
#
# Outputs (saved under results/):
#   - Koller_Endometriosis_Table4_EUR_snps.tsv / Koller_Adenomyosis_clumped_snps.tsv
#   - harmonised_koller_endometriosis_bpo.csv
#   - koller_ivw_results.csv / koller_egger_results.csv /
#     koller_weighted_median_results.csv / koller_all_mr_methods.csv
#   - koller_mr_heterogeneity.csv / koller_mr_pleiotropy.csv /
#     koller_mr_single_snp.csv / koller_mr_leaveoneout_snp.csv
#   - koller_vs_rahmioglu_comparison.csv (all 30 outcomes — see SCOPE UPDATE)
#   - koller_adenomyosis_placenta_praevia_mr.csv
#   - koller_adenomyosis_heterogeneity.csv / koller_adenomyosis_pleiotropy.csv
#   - koller_adenomyosis_all_outcomes_ivw.csv (IVW only, all 30 outcomes,
#     FDR-corrected; Supp Table S9)
################################################################################

### 1) SETUP ####################################################################

required_pkgs <- c("TwoSampleMR", "dplyr", "data.table", "here", "readr", "openxlsx")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

set.seed(42)

data_dir    <- here::here("data")
results_dir <- here::here("results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

koller_dir      <- file.path(data_dir, "EXPOSURE_KOLLER")
endo_table4_file <- file.path(koller_dir, "Koller2026_SupplementaryTable4_significant_loci.xlsx")
adeno_file      <- file.path(koller_dir, "Koller2026_adenomyosis_EUR.tsv.gz")

plink_dir  <- here::here("plink_mac_20241022")
plink_bin  <- file.path(plink_dir, "plink")
bfile_root <- file.path(plink_dir, "24088632", "1000G.EUR.QC")

if (!file.exists(endo_table4_file)) stop("Missing file: ", endo_table4_file)
if (!file.exists(adeno_file))       stop("Missing file: ", adeno_file)

### 2) HELPER: LOAD + DERIVE BETA/SE + FORMAT + CLUMP A KOLLER GWAS ############

prepare_koller_instrument <- function(path, phenotype_label, pval_threshold = 5e-8) {

  message("\n=== Loading Koller GWAS: ", phenotype_label, " ===")
  raw <- data.table::fread(path)

  message("Rows read: ", nrow(raw))

  raw <- raw[raw$P.value < pval_threshold, ]
  message("Genome-wide significant SNPs (p < ", pval_threshold, "): ", nrow(raw))

  raw <- raw |>
    dplyr::mutate(
      Allele1   = toupper(Allele1),
      Allele2   = toupper(Allele2),
      se        = 1 / sqrt(2 * Freq1 * (1 - Freq1) * N),
      beta      = Z * se,
      phenotype = phenotype_label
    )

  fmt <- TwoSampleMR::format_data(
    as.data.frame(raw),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    eaf_col           = "Freq1",
    effect_allele_col = "Allele1",
    other_allele_col  = "Allele2",
    pval_col          = "P.value",
    samplesize_col    = "N",
    phenotype_col     = "phenotype"
  )

  clumped <- TwoSampleMR::clump_data(
    fmt,
    clump_kb  = 10000,
    clump_r2  = 0.001,
    clump_p1  = pval_threshold,
    clump_p2  = 1,
    bfile     = bfile_root,
    plink_bin = plink_bin
  )

  message(phenotype_label, ": ", nrow(clumped), " independent SNPs after clumping.")
  clumped
}

### 3) DERIVE BOTH INSTRUMENTS (ENDOMETRIOSIS + ADENOMYOSIS) ###################
# --- [R3.3 — INSTRUMENT UPDATE, 2026-09-03] Endometriosis instrument now
#     loaded directly from Koller et al.'s published Table 4 (see REVISION
#     LOG / Data source above) instead of self-clumped from the raw
#     genome-wide file. Adenomyosis instrument below is UNCHANGED. ---

message("\n=== Loading Koller et al. 2026 Table 4 (published EUR endometriosis loci) ===")

table4_raw <- openxlsx::read.xlsx(
  endo_table4_file,
  sheet    = "Supplementary Table 4",
  startRow = 6,
  colNames = FALSE
)

# Column positions (1-indexed) confirmed against the sheet's 3-row header
# block: X1=Lead SNP, X2=Chromosome, X3=Position, X4=Effect allele,
# X5=Other allele, X6=Nearest gene, X7=Type, X28-31=EUR "Endometriosis
# combined definition" Beta/SE/P-value/EAF.
table4_endo_eur <- table4_raw |>
  dplyr::filter(grepl("^rs", X1)) |>
  dplyr::transmute(
    SNP           = X1,
    chr           = X2,
    position      = X3,
    effect_allele = toupper(X4),
    other_allele  = toupper(X5),
    nearest_gene  = X6,
    locus_type    = X7,
    beta          = as.numeric(X28),
    se            = as.numeric(X29),
    pval          = as.numeric(X30),
    eaf           = as.numeric(X31),
    phenotype     = "Endometriosis (Koller 2026, Table 4 EUR)"
  )

message("Koller Table 4 rows read: ", nrow(table4_raw),
        " | rs-identified loci: ", nrow(table4_endo_eur),
        " | with non-NA EUR combined-endometriosis effect: ",
        sum(!is.na(table4_endo_eur$beta)))

table4_endo_eur <- table4_endo_eur[!is.na(table4_endo_eur$beta), ]

koller_endo_clumped <- TwoSampleMR::format_data(
  as.data.frame(table4_endo_eur),
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "se",
  eaf_col           = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col  = "other_allele",
  pval_col          = "pval",
  phenotype_col     = "phenotype"
)

message("Endometriosis (Koller 2026, Table 4 EUR) instrument: ",
        nrow(koller_endo_clumped),
        " SNPs — published, already LD-independent (paper's own r2<0.1 ",
        "clumping against 1000G EUR); not re-clumped here.")

data.table::fwrite(
  koller_endo_clumped,
  file.path(results_dir, "Koller_Endometriosis_Table4_EUR_snps.tsv"),
  sep = "\t"
)

koller_adeno_clumped <- prepare_koller_instrument(
  adeno_file, "Adenomyosis (Koller 2026)"
)
data.table::fwrite(
  koller_adeno_clumped,
  file.path(results_dir, "Koller_Adenomyosis_clumped_snps.tsv"),
  sep = "\t"
)

snps_needed <- unique(c(koller_endo_clumped$SNP, koller_adeno_clumped$SNP))
message("\nTotal unique SNPs to look up in outcome GWAS (endo + adeno): ",
        length(snps_needed))

### 4) RE-EXTRACT OUTCOME ASSOCIATIONS GENOME-WIDE (FINNGEN + PPH ONLY) ########
# --- See REVISION LOG: MR-PREG-sourced outcomes are out of scope here. ---

# Helpers reused from 02_prepare_outcomes_endoMR-PREG.R for consistency ------
.se_from_p <- function(beta, pval) {
  p2 <- pmin(pmax(pval / 2, .Machine$double.xmin), 1 - 1e-16)
  abs(beta) / qnorm(p2, lower.tail = FALSE)
}

.to_beta <- function(effect_vec, outcome_name = "") {
  med <- suppressWarnings(median(effect_vec, na.rm = TRUE))
  if (is.finite(med) && med > 0.5 && med < 1.5) {
    message("  ", outcome_name, ": auto-detected OR format (median = ",
            sprintf('%.3f', med), ") -> converting to log OR.")
    log(effect_vec)
  } else {
    message("  ", outcome_name, ": auto-detected beta format (median = ",
            sprintf('%.3f', med), ") -> using as-is.")
    effect_vec
  }
}

message("\n=== Re-extracting FinnGen R12 outcomes at Koller SNPs (genome-wide) ===")

file_list_fg <- c(
  "finngen_R12_O15_PLAC_PRAEVIA",
  "finngen_R12_O15_PLAC_DISORD",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR"
)
common_cols_fg <- c(
  "#chrom", "pos", "ref", "alt", "rsids",
  "pval", "beta", "sebeta",
  "af_alt", "af_alt_cases", "af_alt_controls"
)

process_fg_koller <- function(base) {
  path <- file.path(data_dir, "OUTCOME_FINNGEN", base)
  if (!file.exists(path)) stop("File not found: ", path)

  message("Reading: ", base, " ...")
  raw <- data.table::fread(path, select = common_cols_fg, showProgress = FALSE)

  df <- raw |>
    dplyr::filter(rsids %in% snps_needed) |>
    dplyr::transmute(
      outcome = base,
      rsid    = rsids,
      a1      = alt,
      a2      = ref,
      eaf     = as.numeric(af_alt),
      beta    = as.numeric(beta),
      SE      = as.numeric(sebeta),
      pval    = as.numeric(pval),
      Units   = "logOR",
      n       = NA_integer_
    ) |>
    as.data.frame()

  message("  Matched ", nrow(df), " / ", length(snps_needed), " SNPs for ", base)

  TwoSampleMR::format_data(
    dat               = df,
    type              = "outcome",
    snp_col           = "rsid",
    beta_col          = "beta",
    se_col            = "SE",
    effect_allele_col = "a1",
    other_allele_col  = "a2",
    eaf_col           = "eaf",
    pval_col          = "pval",
    phenotype_col     = "outcome",
    samplesize_col    = "n",
    units_col         = "Units"
  ) |>
    dplyr::mutate(outcome = base, id.outcome = base)
}

fg_koller_outcomes <- dplyr::bind_rows(lapply(file_list_fg, process_fg_koller))

message("\n=== Re-extracting Westergaard PPH outcomes at Koller SNPs (genome-wide) ===")

file_list_pph <- c(
  "Antepartum_bleeding",
  "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta"
)
common_cols_pph <- c(
  "Chr", "PosB38", "Marker", "rsName", "OA", "EA",
  "EAfrq", "Effect", "P", "Phet", "I2", "nCoh"
)

process_pph_koller <- function(base) {
  path <- file.path(data_dir, "OUTCOME_PPH", paste0(base, ".txt"))
  if (!file.exists(path)) stop("File not found: ", path)

  message("Reading: ", base, " ...")
  raw <- data.table::fread(
    path, col.names = common_cols_pph, data.table = FALSE, showProgress = FALSE
  )
  if (nrow(raw) > 0 && all(tolower(as.character(raw[1, ])) == tolower(common_cols_pph))) {
    raw <- raw[-1, , drop = FALSE]
  }

  df <- raw |>
    dplyr::filter(rsName %in% snps_needed) |>
    dplyr::transmute(
      outcome = base,
      rsid    = rsName,
      a1      = EA,
      a2      = OA,
      eaf     = suppressWarnings(as.numeric(EAfrq)),
      beta    = .to_beta(suppressWarnings(as.numeric(Effect)), base),
      pval    = suppressWarnings(as.numeric(P))
    ) |>
    dplyr::mutate(SE = .se_from_p(beta, pval), Units = "logOR", n = NA_integer_) |>
    as.data.frame()

  message("  Matched ", nrow(df), " / ", length(snps_needed), " SNPs for ", base)

  TwoSampleMR::format_data(
    dat               = df,
    type              = "outcome",
    snp_col           = "rsid",
    beta_col          = "beta",
    se_col            = "SE",
    effect_allele_col = "a1",
    other_allele_col  = "a2",
    eaf_col           = "eaf",
    pval_col          = "pval",
    phenotype_col     = "outcome",
    samplesize_col    = "n",
    units_col         = "Units"
  ) |>
    dplyr::mutate(outcome = base, id.outcome = base)
}

pph_koller_outcomes <- dplyr::bind_rows(lapply(file_list_pph, process_pph_koller))

### 4b) LOAD MR-PREG (KOLLER-MATCHED) OUTCOME DATA FOR THE 23 MOBA OUTCOMES ####
# --- [R3.3/R2.5 — SCOPE UPDATE, 2026-07-27] The MR-PREG collaboration (Tuck
#     Seng) has now re-extracted outcome associations at the Koller instrument
#     SNPs for the 23 outcomes previously out of scope (see REVISION LOG
#     above). These 23 outcomes are added here, extending the Koller
#     sensitivity comparison from 7/30 to all 30 main outcomes. ---

message("\n=== Loading MR-PREG (Koller-matched) outcome data for the 23 MoBa outcomes ===")

mrpreg_koller_dir <- file.path(
  data_dir, "EXPOSURE_KOLLER_TABLES4", "Updated_Koller_2026_SNP_list"
)

file_list_mrpreg <- c(
  "rup_memb", "pretb_all", "vpretb_all", "ga_all", "sga", "lbw_all", "hbw_all",
  "lga", "zbw_all", "lowapgar1", "lowapgar5", "nicu", "sb_subsamp",
  "anaemia_preg_all", "gdm_subsamp", "gh_subsamp", "hdp_subsamp", "pe_subsamp",
  "depr_subsamp", "induction", "posttb_all", "el_cs", "em_cs"
)

process_mrpreg_koller <- function(base) {
  path <- file.path(mrpreg_koller_dir, paste0("tmp_Metanalysis_metanalyses-R3.mum.", base, ".txt.gz"))
  if (!file.exists(path)) stop("File not found: ", path)

  message("Reading: ", base, " ...")
  raw <- data.table::fread(path, showProgress = FALSE)

  df <- raw |>
    dplyr::filter(SNP %in% snps_needed) |>
    dplyr::transmute(
      outcome  = base,
      rsid     = SNP,
      a1       = effect_allele,
      a2       = other_allele,
      eaf      = as.numeric(eaf),
      beta     = as.numeric(beta),
      SE       = as.numeric(se),
      pval     = as.numeric(pval),
      n        = as.numeric(samplesize),
      ncase    = if ("ncase" %in% names(raw)) as.numeric(ncase) else NA_real_,
      ncontrol = if ("ncontrol" %in% names(raw)) as.numeric(ncontrol) else NA_real_
    ) |>
    as.data.frame()

  message("  Matched ", nrow(df), " / ", length(snps_needed), " SNPs for ", base)

  TwoSampleMR::format_data(
    dat                = df,
    type               = "outcome",
    snp_col            = "rsid",
    beta_col           = "beta",
    se_col             = "SE",
    effect_allele_col  = "a1",
    other_allele_col   = "a2",
    eaf_col            = "eaf",
    pval_col           = "pval",
    phenotype_col      = "outcome",
    samplesize_col     = "n",
    ncase_col          = "ncase",
    ncontrol_col       = "ncontrol"
  ) |>
    dplyr::mutate(outcome = base, id.outcome = base)
}

mrpreg_koller_outcomes <- dplyr::bind_rows(lapply(file_list_mrpreg, process_mrpreg_koller))

outcome_dat_koller <- dplyr::bind_rows(fg_koller_outcomes, pph_koller_outcomes, mrpreg_koller_outcomes)

vars_keep_koller <- c(file_list_fg, file_list_pph, file_list_mrpreg)

message("\nTotal re-extracted outcome rows (30 outcomes): ", nrow(outcome_dat_koller))

### 5) HARMONISE ENDOMETRIOSIS INSTRUMENT AGAINST THE 7 OUTCOMES ###############

message("\n=== Harmonising Koller endometriosis instrument (all 30 outcomes) ===")

koller_harm <- TwoSampleMR::harmonise_data(
  exposure_dat = koller_endo_clumped,
  outcome_dat  = outcome_dat_koller,
  action       = 2
)
koller_harm <- koller_harm[koller_harm$mr_keep == TRUE, ]

message("Harmonised Koller dataset: ", nrow(koller_harm), " rows across ",
        dplyr::n_distinct(koller_harm$outcome), " outcomes (of ", length(vars_keep_koller), " attempted).")
print(table(koller_harm$outcome))

write.csv(
  koller_harm,
  file.path(results_dir, "harmonised_koller_endometriosis_bpo.csv"),
  row.names = FALSE
)

### 6) MAIN MR ANALYSES + FDR ####################################################

message("\n=== Main MR analyses (Koller endometriosis instrument, 30 outcomes) ===")

koller_ivw   <- TwoSampleMR::mr(koller_harm, method_list = "mr_ivw")
koller_egger <- TwoSampleMR::mr(koller_harm, method_list = "mr_egger_regression")
koller_wm    <- TwoSampleMR::mr(koller_harm, method_list = "mr_weighted_median")
koller_all   <- TwoSampleMR::mr(koller_harm)

koller_ivw$qval   <- p.adjust(koller_ivw$pval,   method = "fdr")
koller_egger$qval <- p.adjust(koller_egger$pval, method = "fdr")
koller_wm$qval    <- p.adjust(koller_wm$pval,    method = "fdr")

write.csv(koller_ivw,   file.path(results_dir, "koller_ivw_results.csv"),            row.names = FALSE)
write.csv(koller_egger, file.path(results_dir, "koller_egger_results.csv"),          row.names = FALSE)
write.csv(koller_wm,    file.path(results_dir, "koller_weighted_median_results.csv"), row.names = FALSE)
write.csv(koller_all,   file.path(results_dir, "koller_all_mr_methods.csv"),          row.names = FALSE)

### 7) STANDARD SENSITIVITY ANALYSES #############################################

message("\n=== Standard sensitivity analyses (Koller endometriosis instrument) ===")

koller_het    <- TwoSampleMR::mr_heterogeneity(koller_harm)
koller_plt    <- TwoSampleMR::mr_pleiotropy_test(koller_harm)
koller_single <- TwoSampleMR::mr_singlesnp(koller_harm)
koller_loo    <- TwoSampleMR::mr_leaveoneout(koller_harm)

write.csv(koller_het,    file.path(results_dir, "koller_mr_heterogeneity.csv"),   row.names = FALSE)
write.csv(koller_plt,    file.path(results_dir, "koller_mr_pleiotropy.csv"),      row.names = FALSE)
write.csv(koller_single, file.path(results_dir, "koller_mr_single_snp.csv"),      row.names = FALSE)
write.csv(koller_loo,    file.path(results_dir, "koller_mr_leaveoneout_snp.csv"), row.names = FALSE)

### 8) COMPARISON WITH THE MAIN (RAHMIOGLU) ANALYSIS #############################
# --- [R3.3] Side-by-side comparison requested by Reviewer #3. Originally
#     restricted to the 7 FinnGen/PPH outcomes with adequate Koller SNP
#     coverage; now covers all 30 main outcomes (see SCOPE UPDATE above). ---

message("\n=== Comparing Koller vs. Rahmioglu instrument (IVW, 30 outcomes) ===")

rahmioglu_ivw_file <- file.path(results_dir, "ivw_results.csv")
if (!file.exists(rahmioglu_ivw_file)) {
  stop("Rahmioglu IVW results not found: ", rahmioglu_ivw_file,
       ". Run 04_main_analyses_endoMR-PREG.R first.")
}
rahmioglu_ivw <- read.csv(rahmioglu_ivw_file, stringsAsFactors = FALSE)
rahmioglu_ivw <- rahmioglu_ivw[rahmioglu_ivw$outcome %in% vars_keep_koller, ]

fmt_or <- function(b, se) {
  sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
}

koller_cmp <- koller_ivw |>
  dplyr::transmute(
    outcome,
    nsnp_koller = nsnp,
    OR_koller   = fmt_or(b, se),
    pval_koller = pval,
    qval_koller = qval
  )

rahmioglu_cmp <- rahmioglu_ivw |>
  dplyr::transmute(
    outcome,
    nsnp_rahmioglu = nsnp,
    OR_rahmioglu   = fmt_or(b, se),
    pval_rahmioglu = as.numeric(pval),
    qval_rahmioglu = as.numeric(qval)
  )

comparison <- dplyr::full_join(rahmioglu_cmp, koller_cmp, by = "outcome") |>
  dplyr::arrange(outcome)

write.csv(
  comparison,
  file.path(results_dir, "koller_vs_rahmioglu_comparison.csv"),
  row.names = FALSE
)
print(comparison)

praevia_cmp <- comparison |> dplyr::filter(outcome == "finngen_R12_O15_PLAC_PRAEVIA")
message(sprintf(
  "\nPlacenta praevia - Rahmioglu: OR %s (p = %.2e) | Koller: OR %s (p = %.2e)",
  praevia_cmp$OR_rahmioglu, praevia_cmp$pval_rahmioglu,
  praevia_cmp$OR_koller,    praevia_cmp$pval_koller
))

### 9) ADENOMYOSIS -> PLACENTA PRAEVIA (MISCLASSIFICATION CHECK) #################
# --- [R2.5] Reviewer #2, comment 5: adenomyosis misclassification ---

message("\n=== Adenomyosis instrument (Koller 2026) -> placenta praevia ===")

praevia_outcome_dat <- outcome_dat_koller[
  outcome_dat_koller$outcome == "finngen_R12_O15_PLAC_PRAEVIA",
]

adeno_harm <- TwoSampleMR::harmonise_data(
  exposure_dat = koller_adeno_clumped,
  outcome_dat  = praevia_outcome_dat,
  action       = 2
)
adeno_harm <- adeno_harm[adeno_harm$mr_keep == TRUE, ]

message("Adenomyosis instrument SNPs available for placenta praevia after harmonisation: ",
        nrow(adeno_harm), " (of ", nrow(koller_adeno_clumped), " clumped instrument SNPs).")

if (nrow(adeno_harm) >= 3) {
  adeno_mr <- TwoSampleMR::mr(
    adeno_harm,
    method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
  )
} else if (nrow(adeno_harm) >= 2) {
  message("Fewer than 3 SNPs: running IVW only (Egger/WM require >= 3).")
  adeno_mr <- TwoSampleMR::mr(adeno_harm, method_list = "mr_ivw")
} else {
  message("Fewer than 2 SNPs available: MR not performed (underpowered).")
  adeno_mr <- data.frame()
}

write.csv(
  adeno_mr,
  file.path(results_dir, "koller_adenomyosis_placenta_praevia_mr.csv"),
  row.names = FALSE
)

# --- [R2.5 / Supp Table S7] Heterogeneity + Egger intercept for the
#     adenomyosis -> placenta praevia exploratory analysis (only ever run for
#     the main Koller endometriosis instrument before; added for completeness
#     of Supplementary Table S7). ---
if (nrow(adeno_harm) >= 3) {
  adeno_het <- TwoSampleMR::mr_heterogeneity(adeno_harm)
  adeno_plt <- TwoSampleMR::mr_pleiotropy_test(adeno_harm)
} else {
  adeno_het <- data.frame()
  adeno_plt <- data.frame()
}
write.csv(adeno_het, file.path(results_dir, "koller_adenomyosis_heterogeneity.csv"), row.names = FALSE)
write.csv(adeno_plt, file.path(results_dir, "koller_adenomyosis_pleiotropy.csv"),   row.names = FALSE)

if (nrow(adeno_mr) > 0) {
  ivw_row <- adeno_mr[adeno_mr$method %in%
                        c("Inverse variance weighted", "Inverse variance weighted (fixed effects)"), ]
  if (nrow(ivw_row) > 0) {
    message(sprintf(
      "Adenomyosis -> placenta praevia, IVW: OR = %s, p = %.3f, nsnp = %d",
      fmt_or(ivw_row$b[1], ivw_row$se[1]), ivw_row$pval[1], ivw_row$nsnp[1]
    ))
  }
}

### 9b) ADENOMYOSIS -> ALL 30 OUTCOMES (IVW SCREEN, SUPP TABLE S9) ###############
# --- Exploratory screen across all 30 outcomes with the adenomyosis
#     instrument. IVW only, given the limited number of instruments (6 SNPs)
#     and the exploratory nature of this analysis; pleiotropy-robust methods
#     are reserved for placenta praevia specifically (section 9 above), the
#     outcome of primary interest for this misclassification check. ---

message("\n=== Adenomyosis instrument (Koller 2026) -> all 30 outcomes (IVW only) ===")

adeno_harm_all <- TwoSampleMR::harmonise_data(
  exposure_dat = koller_adeno_clumped,
  outcome_dat  = outcome_dat_koller,
  action       = 2
)
adeno_harm_all <- adeno_harm_all[adeno_harm_all$mr_keep == TRUE, ]

# ga_all (gestational age) and zbw_all (birthweight z-score) are continuous;
# OR/LCL/UCL are not meaningful for them and are left NA (see script 06,
# continuous_outcomes list, for the beta-scale display of these two rows).
continuous_outcomes_adeno <- c("ga_all", "zbw_all")

adeno_ivw_all <- TwoSampleMR::mr(adeno_harm_all, method_list = "mr_ivw") |>
  dplyr::mutate(
    is_binary = !(outcome %in% continuous_outcomes_adeno),
    OR   = dplyr::if_else(is_binary, exp(b), NA_real_),
    LCL  = dplyr::if_else(is_binary, exp(b - 1.96 * se), NA_real_),
    UCL  = dplyr::if_else(is_binary, exp(b + 1.96 * se), NA_real_),
    qval = p.adjust(pval, method = "fdr")
  ) |>
  dplyr::select(-is_binary) |>
  dplyr::arrange(pval)

write.csv(
  adeno_ivw_all,
  file.path(results_dir, "koller_adenomyosis_all_outcomes_ivw.csv"),
  row.names = FALSE
)
message("Adenomyosis IVW across ", nrow(adeno_ivw_all), " outcomes saved.")

message("\n=== 05.1_koller_sensitivity_endoMR-PREG.R completed ===")
message("NOTE: Koller sensitivity comparison now covers all 30/30 main outcomes ",
        "(FinnGen + PPH + MR-PREG); see SCOPE UPDATE in REVISION LOG.")
