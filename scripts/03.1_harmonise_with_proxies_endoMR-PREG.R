#!/usr/bin/env Rscript
################################################################################
# Script: 03.1_harmonise_with_proxies_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   Proxy-aware extension of script 03. Harmonises Rahmioglu endometriosis
#   instruments with pregnancy outcome GWAS, substituting missing instruments
#   with LD proxies. Does NOT overwrite script 03 outputs.
#
# ── Proxy strategy ───────────────────────────────────────────────────────────
#   For each missing instrument, PLINK (--r2 --ld-snps) identifies the
#   highest-r² proxy in the 1000G EUR reference (r² ≥ PROXY_R2_MIN, within
#   ±PROXY_WINDOW_KB). Substitution occurs only if that single top proxy is
#   present in the local outcome file; no fallback to lower-ranked proxies is
#   attempted. The proxy rsID is relabelled as the target SNP ID so that
#   harmonise_data() merges correctly on the exposure identifier.
#   Original exposure effects are never modified.
#
# ── Allele harmonisation ──────────────────────────────────────────────────────
#   harmonise_data(action = 2): aligns alleles; palindromic SNPs with
#   intermediate MAF (0.42–0.58) are flagged and excluded. Proxy alleles are
#   treated identically to direct-match alleles.
#
# ── Key parameters ────────────────────────────────────────────────────────────
#   PROXY_R2_MIN    = 0.8    (Burgess et al. 2013 threshold)
#   PROXY_WINDOW_KB = 500    (±500 kb search window)
#   Reference panel : 1000G.EUR.QC (same as clumping in script 01; offline)
#
# ── Outputs (results/) ───────────────────────────────────────────────────────
#   harmonised_rahmioglu_bpo.rds          — main harmonised dataset
#   harmonised_rahmioglu_bpo.csv          — same, plain text
#   proxy_summary.csv                     — per-outcome proxy usage
#   instrument_retention_comparison.csv   — before vs after SNP counts
#
# ── Dependencies ─────────────────────────────────────────────────────────────
#   Run 01_select_instruments_endoMR-PREG.R and 02_prepare_outcomes_endoMR-PREG.R
#   PLINK v1.9 binary and 1000G.EUR.QC reference (already in project)
################################################################################


### 1) Setup ###################################################################

required_pkgs <- c(
  "TwoSampleMR",   # harmonise_data(), format_data()
  "dplyr",
  "data.table",
  "here",
  "progress",
  "tibble"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ── Proxy search parameters ───────────────────────────────────────────────────
# r² ≥ 0.8 is the standard MR threshold for proxy validity (Burgess et al. 2013)
# Population: EUR (1000 Genomes European superpopulation, as used in clumping)
# Window: 500 kb either side of the target SNP
PROXY_R2_MIN    <- 0.8    # Minimum r² to accept a proxy
PROXY_WINDOW_KB <- 500    # Search window in kilobases

# ── Directories ───────────────────────────────────────────────────────────────
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")

# ── PLINK paths (same as script 01) ──────────────────────────────────────────
plink_bin   <- file.path(project_dir, "plink_mac_20241022", "plink")
bfile_root  <- file.path(project_dir, "plink_mac_20241022", "24088632", "1000G.EUR.QC")

if (!file.exists(plink_bin)) {
  stop("PLINK binary not found at: ", plink_bin)
}
if (!file.exists(paste0(bfile_root, ".bim"))) {
  stop("1000G reference not found at: ", bfile_root)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

message("=================================================================")
message(" endoMR-PREG: Proxy-aware harmonisation (script 03.1)")
message(" Parameters: r² ≥ ", PROXY_R2_MIN,
        " | pop = EUR (1000G local) | window = ±", PROXY_WINDOW_KB, " kb")
message("=================================================================\n")


################################################################################
# 2) Exposure data                                                              #
################################################################################

message("=== [1/6] Loading exposure instruments ===")

exp_file_candidates <- c(
  file.path(results_dir, "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv"),
  file.path(results_dir, "endometriosis_clumped_snps.tsv")
)
exp_file <- exp_file_candidates[file.exists(exp_file_candidates)][1]

if (is.na(exp_file)) {
  stop(
    "Exposure file not found. Run 01_select_instruments_endoMR-PREG.R first.\n",
    "Expected one of:\n  ", paste(exp_file_candidates, collapse = "\n  ")
  )
}

exposure_dat  <- data.table::fread(exp_file, data.table = FALSE)
exposure_snps <- unique(exposure_dat$SNP)

message("  File:        ", basename(exp_file))
message("  Rows:        ", nrow(exposure_dat))
message("  Unique SNPs: ", length(exposure_snps))


################################################################################
# 3) Helper functions                                                           #
################################################################################

message("\n=== [2/6] Initialising PLINK-based proxy search functions ===")

#' Find LD proxies for a SET of SNPs using PLINK --r2 --ld-snps.
#'
#' Runs a single PLINK call for all target SNPs simultaneously, then for each
#' target keeps the best proxy (highest r²), excluding the SNP itself.
#' Uses the same 1000G.EUR.QC reference as instrument clumping in script 01.
#'
#' @param snps  Character vector of rsIDs to proxy.
#' @param r2    Numeric. Minimum r² (default: PROXY_R2_MIN).
#' @param kb    Numeric. Search window in kb (default: PROXY_WINDOW_KB).
#' @return data.frame(target_snp, proxy_snp, proxy_r2); 0 rows if nothing found.
batch_query_proxies <- function(snps,
                                 r2 = PROXY_R2_MIN,
                                 kb = PROXY_WINDOW_KB) {
  empty <- data.frame(target_snp = character(), proxy_snp = character(),
                      proxy_r2   = numeric(), stringsAsFactors = FALSE)

  if (length(snps) == 0) return(empty)

  message("    PLINK LD proxy search [", length(snps), " SNP(s) | r²≥", r2,
          " | window=", kb, "kb | ref=1000G.EUR.QC] ...")

  # --ld-snps takes a comma-separated list of rsIDs directly (no file support)
  tmp_dir    <- tempdir()
  out_prefix <- file.path(tmp_dir, "proxy_out")

  snps_arg <- paste(snps, collapse = ",")

  # Build PLINK command
  cmd <- paste(
    shQuote(plink_bin),
    "--bfile",        shQuote(bfile_root),
    "--ld-snps",      snps_arg,
    "--ld-window-kb", kb,
    "--ld-window-r2", r2,
    "--ld-window",    99999,        # no variant-count limit within window
    "--r2",
    "--out",          shQuote(out_prefix),
    "--noweb",
    "--silent"
  )

  ret <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ld_file <- paste0(out_prefix, ".ld")
  if (!file.exists(ld_file) || file.size(ld_file) == 0) {
    message("    WARNING: PLINK produced no .ld output (SNPs may not be in reference).")
    return(empty)
  }

  # Parse .ld output: fixed-width whitespace-separated
  ld <- tryCatch(
    data.table::fread(ld_file, data.table = FALSE),
    error = function(e) { message("    WARNING: could not parse .ld file: ", e$message); NULL }
  )
  if (is.null(ld) || nrow(ld) == 0) return(empty)

  names(ld) <- trimws(names(ld))

  # Exclude self-hits (SNP_A == SNP_B, r²=1 by definition)
  ld <- ld[ld$SNP_A != ld$SNP_B, , drop = FALSE]
  if (nrow(ld) == 0) return(empty)

  # For each target SNP, keep BEST proxy (highest r²)
  best <- ld %>%
    dplyr::group_by(SNP_A) %>%
    dplyr::slice_max(order_by = R2, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      target_snp = SNP_A,
      proxy_snp  = SNP_B,
      proxy_r2   = as.numeric(R2)
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  message("    Proxies found: ", nrow(best), " / ", length(snps))
  best
}


#' Extract exposure instruments from a local outcome data.frame,
#' substituting missing instruments with the highest-r² LD proxy from the
#' 1000G EUR reference, if that proxy is present in the local file.
#' Only the single top-ranked proxy per instrument is attempted.
#'
#' The function performs five passes:
#'   1. Direct rsID match (exposure_snps ∩ outcome_snps)
#'   2. PLINK-based LD proxy search (offline; 1000G EUR reference)
#'   3. Check top proxy availability in the local file
#'   4. Extract proxy rows, relabel rsID as target (original exposure SNP ID)
#'   5. Combine, deduplicate (direct match wins over proxy if both present)
#'
#' @param df_full    data.frame. Complete outcome file (ALL SNPs, not filtered).
#' @param rsid_col   Character.  Name of the rsID column in df_full.
#' @param outcome_nm Character.  Outcome label for console messages.
#' @param proxy_tbl  Optional data.frame (target_snp, proxy_snp, proxy_r2).
#'                   If supplied, the PLINK query step is skipped (re-use a
#'                   pre-computed proxy table — avoids redundant PLINK calls).
#' @return Named list:
#'   $data      — data.frame of instruments (+proxies) with tracking columns:
#'                proxy_used (logical), target_snp, proxy_snp, proxy_r2
#'   $proxy_log — data.frame recording all proxy candidates and their outcome
#'                availability (outcome, target_snp, proxy_snp, proxy_r2,
#'                in_outcome)
extract_instruments_with_proxies <- function(df_full,
                                              rsid_col,
                                              outcome_nm,
                                              proxy_tbl = NULL) {

  all_snps_in_file <- unique(df_full[[rsid_col]])

  # ── Pass 1: Direct match ───────────────────────────────────────────────────
  df_direct    <- df_full[df_full[[rsid_col]] %in% exposure_snps, , drop = FALSE]
  present_snps <- unique(df_direct[[rsid_col]])
  missing_snps <- setdiff(exposure_snps, present_snps)

  message("  [", outcome_nm, "]")
  message("    Direct match:   ", length(present_snps), " / ", length(exposure_snps),
          "  |  Missing: ", length(missing_snps))

  # Annotate direct-match rows with proxy tracking columns
  df_direct$proxy_used <- FALSE
  df_direct$target_snp <- df_direct[[rsid_col]]
  df_direct$proxy_snp  <- NA_character_
  df_direct$proxy_r2   <- NA_real_

  if (length(missing_snps) == 0) {
    proxy_log <- data.frame(
      outcome    = outcome_nm, target_snp = character(),
      proxy_snp  = character(), proxy_r2   = numeric(),
      in_outcome = logical(), stringsAsFactors = FALSE
    )
    return(list(data = df_direct, proxy_log = proxy_log))
  }

  # ── Pass 2: LD proxy query (or accept pre-computed table) ─────────────────
  if (is.null(proxy_tbl)) {
    proxy_candidates <- batch_query_proxies(missing_snps)
  } else {
    # Use pre-computed table, restricting to missing_snps for this outcome
    proxy_candidates <- proxy_tbl[proxy_tbl$target_snp %in% missing_snps, ,
                                  drop = FALSE]
  }

  # Build full proxy log (one row per missing SNP, regardless of outcome status)
  proxy_log_all <- data.frame(
    outcome    = outcome_nm,
    target_snp = missing_snps,
    proxy_snp  = NA_character_,
    proxy_r2   = NA_real_,
    in_outcome = FALSE,
    stringsAsFactors = FALSE
  )
  if (nrow(proxy_candidates) > 0) {
    match_idx <- match(proxy_candidates$target_snp, proxy_log_all$target_snp)
    valid_idx <- !is.na(match_idx)
    proxy_log_all$proxy_snp[match_idx[valid_idx]] <- proxy_candidates$proxy_snp[valid_idx]
    proxy_log_all$proxy_r2[match_idx[valid_idx]]  <- proxy_candidates$proxy_r2[valid_idx]
  }

  if (nrow(proxy_candidates) == 0) {
    message("    No proxy candidates returned by PLINK.")
    return(list(data = df_direct, proxy_log = proxy_log_all))
  }

  # ── Pass 3: Check proxy availability in local file ─────────────────────────
  proxy_candidates$in_outcome <- proxy_candidates$proxy_snp %in% all_snps_in_file
  proxy_log_all$in_outcome[
    match(proxy_candidates$target_snp, proxy_log_all$target_snp)
  ] <- proxy_candidates$in_outcome

  proxy_found <- proxy_candidates[proxy_candidates$in_outcome, , drop = FALSE]

  message("    Proxies found in local file: ",
          nrow(proxy_found), " / ", nrow(proxy_candidates))

  if (nrow(proxy_found) == 0) {
    return(list(data = df_direct, proxy_log = proxy_log_all))
  }

  # ── Pass 4: Extract proxy rows from the full file, relabel as target ───────
  # Relabelling the proxy rsID as the target SNP ID ensures that
  # harmonise_data() merges correctly on the exposure SNP identifier.
  # Proxy alleles will still be correctly harmonised by TwoSampleMR.
  proxy_rows_list <- lapply(seq_len(nrow(proxy_found)), function(i) {
    tgt <- proxy_found$target_snp[i]
    prx <- proxy_found$proxy_snp[i]
    r2v <- proxy_found$proxy_r2[i]

    rows <- df_full[df_full[[rsid_col]] == prx, , drop = FALSE]
    if (nrow(rows) == 0) return(NULL)

    # Relabel: proxy SNP becomes the target (original exposure instrument ID)
    rows[[rsid_col]] <- tgt
    rows$proxy_used  <- TRUE
    rows$target_snp  <- tgt
    rows$proxy_snp   <- prx
    rows$proxy_r2    <- r2v
    rows
  })
  proxy_rows <- dplyr::bind_rows(proxy_rows_list)

  # ── Pass 5: Combine + deduplicate (direct match takes precedence) ──────────
  combined <- dplyr::bind_rows(df_direct, proxy_rows) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(rsid_col))) %>%
    dplyr::arrange(proxy_used) %>%   # FALSE (direct) sorts before TRUE (proxy)
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    as.data.frame(stringsAsFactors = FALSE)

  n_proxy_sub <- sum(combined$proxy_used, na.rm = TRUE)
  message("    Proxy substitutions made:   ", n_proxy_sub)
  message("    Instruments retained total: ",
          nrow(combined), " / ", length(exposure_snps))

  list(data = combined, proxy_log = proxy_log_all)
}


################################################################################
# 4) Load baseline SNP counts (from script 02, for comparison)                 #
################################################################################

message("\n=== [3/6] Loading baseline outcome data (script 02 outputs) ===")

out_file_base <- file.path(results_dir, "formatted_all_pregnancy_outcomes.tsv")
if (!file.exists(out_file_base)) {
  stop(
    "Baseline outcome file not found. ",
    "Run 02_prepare_outcomes_endoMR-PREG.R first.\n",
    out_file_base
  )
}

outcome_dat_base <- data.table::fread(out_file_base, data.table = FALSE)

before_counts <- outcome_dat_base %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarise(n_snp_before = dplyr::n_distinct(SNP), .groups = "drop") %>%
  as.data.frame()

message("  Baseline counts loaded: ", nrow(before_counts), " outcomes.")
message("  Exposure instrument set: ", length(exposure_snps), " SNPs.")


################################################################################
# 5) Process each outcome source with proxy substitution                       #
################################################################################

message("\n=== [4/6] Extracting outcome data with proxy substitution ===")

all_proxy_logs <- list()   # proxy log accumulator (one entry per outcome)

# ── Shared helpers ────────────────────────────────────────────────────────────

# Derive SE from two-sided p-value (used where SE is absent from raw PPH files)
.se_from_p <- function(beta, pval) {
  p2 <- pmin(pmax(pval / 2, .Machine$double.xmin), 1 - 1e-16)
  abs(beta) / qnorm(p2, lower.tail = FALSE)
}

# Auto-detect OR vs log-OR format and convert to log-OR where needed
.to_log_or <- function(effect_vec, nm = "") {
  med <- suppressWarnings(median(effect_vec, na.rm = TRUE))
  if (is.finite(med) && med > 0.5 && med < 1.5) {
    message("    [", nm, "] OR format detected (median = ",
            round(med, 3), ") → converting to log OR.")
    log(effect_vec)
  } else {
    effect_vec
  }
}


### 5a) MR-PREG outcomes -------------------------------------------------------
# The MR-PREG file contains multiple phenotypes.
# Proxy search is done ONCE on the full file, then applied per phenotype.

message("\n--- MR-PREG adverse pregnancy outcomes ---")

mrpreg_path  <- file.path(data_dir, "OUTCOME_MR-PREG", "ma_out_dat.txt")
raw_lines    <- readLines(mrpreg_path)
header_line  <- raw_lines[1]
data_clean   <- raw_lines[-1][!grepl("\\.png$", raw_lines[-1])]
mrpreg_raw   <- read.table(
  text             = c(header_line, data_clean),
  header           = TRUE,
  stringsAsFactors = FALSE
)

# Query proxies ONCE for SNPs missing from the full MR-PREG file
all_mrpreg_snps  <- unique(mrpreg_raw$SNP)
missing_mrpreg   <- setdiff(exposure_snps, all_mrpreg_snps)
message("  MR-PREG global missing SNPs (all phenotypes): ", length(missing_mrpreg))
mrpreg_proxy_tbl <- batch_query_proxies(missing_mrpreg)   # single PLINK call

mrpreg_phenotypes <- unique(mrpreg_raw$Phenotype)
message("  Processing ", length(mrpreg_phenotypes), " MR-PREG phenotypes ...")

mrpreg_results <- lapply(mrpreg_phenotypes, function(ph) {
  df_ph <- mrpreg_raw[mrpreg_raw$Phenotype == ph, , drop = FALSE]
  res   <- extract_instruments_with_proxies(
    df_full    = df_ph,
    rsid_col   = "SNP",
    outcome_nm = ph,
    proxy_tbl  = mrpreg_proxy_tbl   # re-use pre-computed table — no extra PLINK call
  )
  all_proxy_logs[[ph]] <<- res$proxy_log
  res$data
})
mrpreg_combined <- dplyr::bind_rows(mrpreg_results)

# Cache proxy metadata BEFORE format_data() strips custom columns
proxy_meta_mrpreg <- mrpreg_combined %>%
  dplyr::select(SNP, Phenotype, proxy_used, target_snp, proxy_snp, proxy_r2) %>%
  dplyr::distinct(SNP, Phenotype, .keep_all = TRUE)

mr_preg_dat <- TwoSampleMR::format_data(
  dat               = mrpreg_combined,
  type              = "outcome",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "se",
  eaf_col           = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col  = "other_allele",
  pval_col          = "pval",
  samplesize_col    = "samplesize",
  chr_col           = "chr",
  pos_col           = "pos_b38",
  phenotype_col     = "Phenotype",
  id_col            = "StudyID"
)

# Re-attach proxy metadata (format_data does not guarantee custom col retention)
mr_preg_dat <- dplyr::left_join(
  mr_preg_dat,
  proxy_meta_mrpreg %>% dplyr::rename(outcome = Phenotype),
  by = c("SNP", "outcome")
)
mr_preg_dat <- mr_preg_dat %>%
  dplyr::select(-dplyr::any_of(c("chr.outcome", "pos.outcome")))


### 5b) PPH outcomes (Westergaard) --------------------------------------------

message("\n--- Postpartum haemorrhage outcomes (Westergaard et al.) ---")

file_list_pph <- c(
  "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta",
  "Antepartum_bleeding"
  # Early_bleeding_with_any_outcome and Early_bleeding_ending_in_live_birth excluded:
  # analysis restricted to 2nd/3rd trimester outcomes
)

common_cols_pph <- c(
  "Chr", "PosB38", "Marker", "rsName", "OA", "EA",
  "EAfrq", "Effect", "P", "Phet", "I2", "nCoh"
)

process_pph_with_proxies <- function(base) {
  path <- file.path(data_dir, "OUTCOME_PPH", paste0(base, ".txt"))

  raw <- data.table::fread(
    path, col.names = common_cols_pph,
    data.table = FALSE, showProgress = FALSE
  )
  # Drop any repeated header row present in some PPH files
  if (nrow(raw) > 0 &&
      all(tolower(as.character(raw[1, ])) == tolower(common_cols_pph))) {
    raw <- raw[-1, , drop = FALSE]
  }

  # Build the FULL data.frame (all SNPs) needed for proxy row extraction
  df_full <- raw %>%
    dplyr::transmute(
      outcome = base,
      rsid    = rsName,
      chr     = gsub("^chr", "", as.character(Chr)),
      pos     = suppressWarnings(as.integer(PosB38)),
      a1      = EA,
      a2      = OA,
      eaf     = suppressWarnings(as.numeric(EAfrq)),
      beta    = .to_log_or(suppressWarnings(as.numeric(Effect)), base),
      pval    = suppressWarnings(as.numeric(P))
    ) %>%
    dplyr::mutate(
      SE    = .se_from_p(beta, pval),
      Units = "logOR",
      Gene  = NA_character_,
      n     = NA_integer_
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  res <- extract_instruments_with_proxies(
    df_full    = df_full,
    rsid_col   = "rsid",
    outcome_nm = base
  )
  all_proxy_logs[[base]] <<- res$proxy_log

  # Cache proxy metadata before format_data()
  proxy_meta <- res$data %>%
    dplyr::select(rsid, proxy_used, target_snp, proxy_snp, proxy_r2) %>%
    dplyr::distinct(rsid, .keep_all = TRUE)

  mr_out <- TwoSampleMR::format_data(
    dat               = res$data,
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
    units_col         = "Units",
    gene_col          = "Gene",
    chr_col           = "chr",
    pos_col           = "pos"
  )
  mr_out$outcome    <- base
  mr_out$id.outcome <- base

  dplyr::left_join(mr_out, proxy_meta, by = c("SNP" = "rsid"))
}

message("Processing ", length(file_list_pph), " PPH outcomes ...")
res_pph_proxied <- lapply(file_list_pph, process_pph_with_proxies)
names(res_pph_proxied) <- file_list_pph
pph_dat <- dplyr::bind_rows(res_pph_proxied)


### 5c) FinnGen R12 outcomes --------------------------------------------------

message("\n--- FinnGen R12 outcomes ---")

file_list_fg <- c(
  "finngen_R12_O15_PLAC_PRAEVIA",
  "finngen_R12_O15_PLAC_DISORD",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",
  "finngen_R12_O15_PREG_ECTOP",
  "finngen_R12_N14_FEMALEINFERT"
)

common_cols_fg <- c(
  "#chrom", "pos", "ref", "alt", "rsids",
  "pval", "beta", "sebeta",
  "af_alt", "af_alt_cases", "af_alt_controls"
)

process_fg_with_proxies <- function(base) {
  path <- file.path(data_dir, "OUTCOME_FINNGEN", base)

  raw <- data.table::fread(path, select = common_cols_fg, showProgress = FALSE)

  df_full <- raw %>%
    dplyr::transmute(
      outcome = base,
      rsid    = rsids,
      chr     = as.character(`#chrom`),
      pos     = pos,
      a1      = alt,
      a2      = ref,
      eaf     = suppressWarnings(as.numeric(af_alt)),
      beta    = suppressWarnings(as.numeric(beta)),
      SE      = suppressWarnings(as.numeric(sebeta)),
      pval    = suppressWarnings(as.numeric(pval)),
      Units   = "logOR",
      Gene    = NA_character_,
      n       = NA_integer_
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  res <- extract_instruments_with_proxies(
    df_full    = df_full,
    rsid_col   = "rsid",
    outcome_nm = base
  )
  all_proxy_logs[[base]] <<- res$proxy_log

  proxy_meta <- res$data %>%
    dplyr::select(rsid, proxy_used, target_snp, proxy_snp, proxy_r2) %>%
    dplyr::distinct(rsid, .keep_all = TRUE)

  mr_out <- TwoSampleMR::format_data(
    dat               = res$data,
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
    units_col         = "Units",
    gene_col          = "Gene"
  )
  mr_out$outcome    <- base
  mr_out$id.outcome <- base

  dplyr::left_join(mr_out, proxy_meta, by = c("SNP" = "rsid"))
}

message("Processing ", length(file_list_fg), " FinnGen outcomes ...")
res_fg_proxied <- lapply(file_list_fg, process_fg_with_proxies)
names(res_fg_proxied) <- file_list_fg
fg_dat <- dplyr::bind_rows(res_fg_proxied)


### 5d) Combine all outcomes --------------------------------------------------

all_outcomes_proxied <- dplyr::bind_rows(mr_preg_dat, pph_dat, fg_dat)

message("\n  Combined proxied outcome dataset:")
message("    Rows:          ", nrow(all_outcomes_proxied))
message("    Unique outcomes:", dplyr::n_distinct(all_outcomes_proxied$outcome))
message("    Unique SNPs:    ", dplyr::n_distinct(all_outcomes_proxied$SNP))
message("    Proxy rows:     ",
        sum(all_outcomes_proxied$proxy_used, na.rm = TRUE))


################################################################################
# 6) Harmonisation                                                              #
################################################################################

message("\n=== [5/6] Harmonising (action = 2) ===")
message("  action = 2: align alleles; flag palindromic SNPs with MAF 0.42–0.58.")
message("  Exposure effects are NEVER modified — only outcome alleles are flipped.")
message("  Proxy alleles are treated identically to direct-match alleles.")

dat_harmonised <- TwoSampleMR::harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat  = all_outcomes_proxied,
  action       = 2
)

# ── Re-attach proxy tracking metadata (robust re-join) ────────────────────────
# harmonise_data() may or may not propagate custom outcome_dat columns.
# If it does, a subsequent left_join creates .x/.y name conflicts.
# Solution: explicitly drop any existing proxy cols first, then re-join cleanly.
proxy_tracking <- all_outcomes_proxied %>%
  dplyr::select(SNP, outcome, proxy_used, target_snp, proxy_snp, proxy_r2) %>%
  dplyr::distinct(SNP, outcome, .keep_all = TRUE)

dat_harmonised <- dat_harmonised %>%
  dplyr::select(-dplyr::any_of(c("proxy_used", "target_snp", "proxy_snp", "proxy_r2")))

dat_harmonised <- dplyr::left_join(
  dat_harmonised,
  proxy_tracking,
  by = c("SNP", "outcome")
)

dat_harmonised$proxy_used <- ifelse(
  is.na(dat_harmonised$proxy_used), FALSE, dat_harmonised$proxy_used
)

# ── Apply mr_keep filter ──────────────────────────────────────────────────────
dat_final <- dat_harmonised[dat_harmonised$mr_keep == TRUE, ]

message("  Rows retained after mr_keep: ", nrow(dat_final))

# Print per-outcome SNP counts
outcome_counts_after <- dat_final %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarise(n_snp_after = dplyr::n_distinct(SNP), .groups = "drop") %>%
  as.data.frame()

message("\n  SNPs per outcome (post-harmonisation with proxies):")
print(outcome_counts_after, row.names = FALSE)


################################################################################
# 7) Summary tables                                                             #
################################################################################

message("\n=== [6/6] Generating summary tables ===")

### 7a) proxy_summary.csv ------------------------------------------------------
# One row per outcome:
#   outcome           — outcome label
#   n_snp_original    — number of exposure instruments (always = length(exposure_snps))
#   n_snp_retained    — distinct SNPs retained after harmonisation
#   n_proxies_used    — how many of the retained rows used a proxy
#   pct_retained      — percentage of instruments retained
#   proxy_snps        — semicolon-separated list of proxy SNPs used

proxy_summary <- dat_final %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarise(
    n_snp_original = length(exposure_snps),
    n_snp_retained = dplyr::n_distinct(SNP),
    n_proxies_used = sum(proxy_used, na.rm = TRUE),
    pct_retained   = round(100 * dplyr::n_distinct(SNP) / length(exposure_snps), 1),
    proxy_snps     = paste(
      unique(proxy_snp[!is.na(proxy_snp) & proxy_used]),
      collapse = "; "
    ),
    .groups = "drop"
  )

write.csv(
  proxy_summary,
  file.path(results_dir, "proxy_summary.csv"),
  row.names = FALSE
)
message("  Saved: proxy_summary.csv")
print(proxy_summary %>% dplyr::select(-proxy_snps), n = Inf)


### 7b) instrument_retention_comparison.csv ------------------------------------
# Side-by-side comparison of SNP counts before (script 02, no proxies) and
# after (this script, with proxies) harmonisation.

retention_comparison <- dplyr::full_join(
  before_counts,
  outcome_counts_after,
  by = "outcome"
) %>%
  dplyr::mutate(
    n_snp_before         = ifelse(is.na(n_snp_before), 0L, n_snp_before),
    n_snp_after          = ifelse(is.na(n_snp_after),  0L, n_snp_after),
    n_total_instruments  = length(exposure_snps),
    snp_gain             = n_snp_after - n_snp_before,
    pct_before           = round(100 * n_snp_before / length(exposure_snps), 1),
    pct_after            = round(100 * n_snp_after  / length(exposure_snps), 1),
    pct_gain             = round(100 * snp_gain     / length(exposure_snps), 1)
  ) %>%
  dplyr::arrange(dplyr::desc(snp_gain), outcome) %>%
  dplyr::select(
    outcome, n_total_instruments,
    n_snp_before, pct_before,
    n_snp_after,  pct_after,
    snp_gain, pct_gain
  )

write.csv(
  retention_comparison,
  file.path(results_dir, "instrument_retention_comparison.csv"),
  row.names = FALSE
)
message("  Saved: instrument_retention_comparison.csv")
print(as.data.frame(retention_comparison), row.names = FALSE)


################################################################################
# 8) Save harmonised dataset                                                   #
################################################################################

saveRDS(
  dat_final,
  file.path(results_dir, "harmonised_rahmioglu_bpo.rds")
)
message("\n  Saved: harmonised_rahmioglu_bpo.rds")

write.csv(
  dat_final,
  file.path(results_dir, "harmonised_rahmioglu_bpo.csv"),
  row.names = FALSE
)
message("  Saved: harmonised_rahmioglu_bpo.csv")


################################################################################
# 9) Session information                                                        #
################################################################################

message("\n=================================================================")
message(" COMPLETED: endoMR-PREG proxy-aware harmonisation")
message("=================================================================")
message("\nOutputs in: ", results_dir)
message("  harmonised_rahmioglu_bpo.rds          (main analysis object)")
message("  harmonised_rahmioglu_bpo.csv          (same, plain text)")
message("  proxy_summary.csv                     (per-outcome proxy usage)")
message("  instrument_retention_comparison.csv   (before vs after SNP counts)")
message("\nSession info:")
print(sessionInfo())


# If any outcome source becomes available on IEU OpenGWAS, use
# TwoSampleMR::extract_outcome_data(proxies = TRUE, rsq = PROXY_R2_MIN) instead.
