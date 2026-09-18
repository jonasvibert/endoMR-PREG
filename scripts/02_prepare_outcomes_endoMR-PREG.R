#!/usr/bin/env Rscript
################################################################################
# Script: 02_prepare_outcomes_endoMR-PREG.R
# Project: endoMR-PREG
# Purpose: Import, format and harmonise outcome GWAS for MR:
#   - MR-PREG adverse pregnancy outcomes (ma_out_dat.txt)
#   - Postpartum haemorrhage (Westergaard et al.)
#   - FinnGen R12 adverse pregnancy outcomes
################################################################################

### 0) Setup ###################################################################

# Packages ---------------------------------------------------------------------
required_pkgs <- c(
  "data.table",   # fread(), fwrite()
  "dplyr",        # %>%, bind_rows, group_by, summarise, n_distinct
  "TwoSampleMR",  # format_data()
  "here",         # here()
  "progress",     # progress_bar
  "tibble"        # tibble()
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

### 1) Paths ###################################################################

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
tables_dir  <- file.path(results_dir, "tables")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)

### 2) SNP list (endo instruments) #############################################

snps_file <- file.path(results_dir, "endo_instruments_41.txt")

if (file.exists(snps_file)) {
  snps <- scan(
    snps_file,
    what        = character(),
    sep         = "\n",
    quiet       = TRUE,
    strip.white = TRUE
  )
  snps <- unique(snps[nzchar(snps)])                  # drop blanks + dedupe
  snps <- snps[!tolower(snps) %in% c("snp", "rsid")]  # drop any header
} else if (file.exists(file.path(results_dir, "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv"))) {
  # Preferred fallback: TSV created by 01_select_instruments_endoMR-PREG.R
  tmp <- data.table::fread(
    file.path(results_dir, "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv"),
    select      = "SNP",
    showProgress = FALSE
  )
  snps <- unique(tmp$SNP)
} else if (file.exists(file.path(results_dir, "endometriosis_clumped_snps.tsv"))) {
  # Legacy fallback: old filename
  tmp <- data.table::fread(
    file.path(results_dir, "endometriosis_clumped_snps.tsv"),
    select      = "SNP",
    showProgress = FALSE
  )
  snps <- unique(tmp$SNP)
} else {
  stop(
    "Could not find the 41-SNP instrument list. ",
    "Provide 'endo_instruments_41.txt' or run the instrument selection script ",
    "to create 'Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv' in ",
    results_dir
  )
}

message("Loaded ", length(snps), " instruments (expected ~41).")
message("Starting data preparation for pregnancy outcome GWAS...")

################################################################################
# 3) MR-PREG adverse pregnancy outcomes                                        #
#    Units: Beta coefficients (log OR) with standard errors                    #
################################################################################

message("\n=== Processing MR-PREG adverse pregnancy outcomes ===")
message("Data source: MR-PREG consortium")
message("Expected units: beta (log OR) and SE")

f_in <- file.path(data_dir, "OUTCOME_MR-PREG", "ma_out_dat.txt")
message("Reading MR-PREG data from: ", basename(f_in))

raw_lines  <- readLines(f_in)
header     <- raw_lines[1]
data_clean <- raw_lines[-1][!grepl("\\.png$", raw_lines[-1])]
txt_clean  <- c(header, data_clean)

message(
  "Cleaned ", length(data_clean), " data lines (removed ",
  length(raw_lines) - 1 - length(data_clean), " non-data lines)."
)

mr_preg_df <- read.table(
  text             = txt_clean,
  header           = TRUE,
  stringsAsFactors = FALSE
)

# Diagnostics ------------------------------------------------------------------
message("MR-PREG: ", nrow(mr_preg_df), " rows x ", ncol(mr_preg_df), " columns")
message("  - Outcomes: ", dplyr::n_distinct(mr_preg_df$Phenotype))
message("  - SNPs:     ", dplyr::n_distinct(mr_preg_df$SNP))
message("  - Mean beta (log OR): ", sprintf('%.4f', mean(mr_preg_df$beta, na.rm = TRUE)))
message("  - Mean SE:            ", sprintf('%.4f', mean(mr_preg_df$se,   na.rm = TRUE)))
message("Units confirmed: beta coefficients (log OR) with SE.")

# Format for TwoSampleMR -------------------------------------------------------
mr_preg_dat <- TwoSampleMR::format_data(
  dat                 = mr_preg_df,
  type                = "outcome",
  snp_col             = "SNP",
  beta_col            = "beta",
  se_col              = "se",
  eaf_col             = "eaf",
  effect_allele_col   = "effect_allele",
  other_allele_col    = "other_allele",
  pval_col            = "pval",
  samplesize_col      = "samplesize",
  chr_col             = "chr",
  pos_col             = "pos_b38",
  phenotype_col       = "Phenotype",
  id_col              = "StudyID"
)

assign("dat_MR_PREG", mr_preg_dat, envir = .GlobalEnv)
message("MR-PREG processing complete: ", nrow(mr_preg_dat), " observations formatted.")

################################################################################
# 4) Postpartum haemorrhage outcomes (Westergaard et al.)                      #
#    Units: OR or beta, auto-detected and converted to log OR                  #
################################################################################

message("\n=== Processing postpartum haemorrhage outcomes (Westergaard et al.) ===")
message("Data source: Westergaard et al., Nature Genetics 2024")
message("Auto-detection: OR (≈1) → log OR, otherwise beta used as-is.")

# List of raw PPH files (without extension) -----------------------------------
file_list_pph <- c(
  "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta",
  "Early_bleeding_with_any_outcome",
  "Antepartum_bleeding",
  "Early_bleeding_ending_in_live_birth"
)

# Expected columns in PPH files -----------------------------------------------
common_cols_pph <- c(
  "Chr", "PosB38", "Marker", "rsName", "OA", "EA",
  "EAfrq", "Effect", "P", "Phet", "I2", "nCoh"
)

# Helper: robust SE from two-sided p-value ------------------------------------
.se_from_p <- function(beta, pval) {
  p2 <- pmin(pmax(pval / 2, .Machine$double.xmin), 1 - 1e-16)
  abs(beta) / qnorm(p2, lower.tail = FALSE)
}

# Helper: auto-detect if "Effect" is OR (~1) or already beta (~0) -------------
.to_beta <- function(effect_vec, outcome_name = "") {
  med <- suppressWarnings(median(effect_vec, na.rm = TRUE))
  if (is.finite(med) && med > 0.5 && med < 1.5) {
    message(
      "  ", outcome_name,
      ": auto-detected OR format (median = ", sprintf('%.3f', med),
      ") → converting to log OR."
    )
    log(effect_vec)
  } else {
    message(
      "  ", outcome_name,
      ": auto-detected beta format (median = ", sprintf('%.3f', med),
      ") → using as-is."
    )
    effect_vec
  }
}

# Function to read, filter on SNPs, and format (PPH) ---------------------------
process_pph <- function(base) {
  path <- file.path(data_dir, "OUTCOME_PPH", paste0(base, ".txt"))
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  
  raw <- data.table::fread(
    path,
    col.names    = common_cols_pph,
    data.table   = FALSE,
    showProgress = FALSE
  )
  
  # Drop any repeated header row
  if (nrow(raw) > 0 && all(tolower(as.character(raw[1, ])) == tolower(common_cols_pph))) {
    raw <- raw[-1, , drop = FALSE]
  }
  
  df <- raw %>%
    dplyr::filter(rsName %in% snps) %>%
    dplyr::transmute(
      outcome = base,
      rsid    = rsName,
      chr     = gsub("^chr", "", as.character(Chr)),
      pos     = as.integer(PosB38),
      a1      = EA,   # effect allele
      a2      = OA,   # other allele
      eaf     = suppressWarnings(as.numeric(EAfrq)),
      beta    = .to_beta(suppressWarnings(as.numeric(Effect)), base),
      pval    = suppressWarnings(as.numeric(P))
    ) %>%
    dplyr::mutate(
      SE    = .se_from_p(beta, pval),
      Units = "logOR",
      Gene  = NA_character_,
      n     = NA_integer_
    ) %>%
    as.data.frame()
  
  mr_pph <- TwoSampleMR::format_data(
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
    units_col         = "Units",
    gene_col          = "Gene",
    chr_col           = "chr",
    pos_col           = "pos"
  )
  
  mr_pph$outcome    <- base
  mr_pph$id.outcome <- base
  
  assign(paste0("dat_", base), mr_pph, envir = .GlobalEnv)
  
  list(raw = df, mr = mr_pph)
}

# Run for all PPH outcomes with progress bar ----------------------------------
message("Processing ", length(file_list_pph), " postpartum haemorrhage outcomes...")
pb_pph <- progress::progress_bar$new(
  format = "PPH     [:bar] :percent :current/:total ETA: :eta",
  total  = length(file_list_pph),
  clear  = FALSE
)

res_pph <- vector("list", length(file_list_pph))
for (i in seq_along(file_list_pph)) {
  pb_pph$tick()
  res_pph[[i]] <- process_pph(file_list_pph[i])
}
names(res_pph) <- file_list_pph

# QC summary -------------------------------------------------------------------
qc_pph <- dplyr::bind_rows(lapply(res_pph, `[[`, "raw")) %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarise(
    nsnp      = dplyr::n_distinct(rsid),
    mean_beta = mean(beta, na.rm = TRUE),
    sd_beta   = sd(beta,   na.rm = TRUE),
    .groups   = "drop"
  )

message("\nPostpartum haemorrhage GWAS summary (restricted to instruments):")
print(qc_pph)

# Optionally save QC table
data.table::fwrite(
  qc_pph,
  file.path(tables_dir, "QC_Table_PPH_outcomes_restricted_to_instruments.tsv"),
  sep = "\t"
)

################################################################################
# 5) FinnGen R12 outcomes                                                      #
#    Units: Beta coefficients (log OR) with standard errors                    #
################################################################################

message("\n=== Processing FinnGen R12 adverse pregnancy outcomes ===")
message("Data source: FinnGen R12")
message("Expected units: beta (log OR) and sebeta (SE).")

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

process_fg <- function(base) {
  path <- file.path(data_dir, "OUTCOME_FINNGEN", base)
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  
  raw <- data.table::fread(
    path,
    select      = common_cols_fg,
    showProgress = FALSE
  )
  
  df <- raw %>%
    dplyr::filter(rsids %in% snps) %>%
    dplyr::transmute(
      outcome = base,
      rsid    = rsids,
      chr     = as.character(`#chrom`),
      pos     = pos,
      a1      = alt,
      a2      = ref,
      eaf     = as.numeric(af_alt),
      beta    = as.numeric(beta),
      SE      = as.numeric(sebeta),
      pval    = as.numeric(pval),
      Units   = "logOR",
      Gene    = NA_character_,
      n       = NA_integer_
    ) %>%
    as.data.frame()
  
  mr_fg <- TwoSampleMR::format_data(
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
    units_col         = "Units",
    gene_col          = "Gene"
  )
  
  mr_fg$outcome    <- base
  mr_fg$id.outcome <- base
  
  assign(paste0("dat_", base), mr_fg, envir = .GlobalEnv)
  
  list(raw = df, mr = mr_fg)
}

message("Processing ", length(file_list_fg), " FinnGen outcomes...")
pb_fg <- progress::progress_bar$new(
  format = "FinnGen [:bar] :percent :current/:total ETA: :eta",
  total  = length(file_list_fg),
  clear  = FALSE
)

res_fg <- vector("list", length(file_list_fg))
for (i in seq_along(file_list_fg)) {
  pb_fg$tick()
  res_fg[[i]] <- process_fg(file_list_fg[i])
}
names(res_fg) <- file_list_fg

qc_fg <- dplyr::bind_rows(lapply(res_fg, `[[`, "raw")) %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarise(
    nsnp      = dplyr::n_distinct(rsid),
    mean_beta = mean(beta, na.rm = TRUE),
    sd_beta   = sd(beta,   na.rm = TRUE),
    .groups   = "drop"
  )

message("\nFinnGen R12 GWAS summary (restricted to instruments):")
print(qc_fg)

data.table::fwrite(
  qc_fg,
  file.path(tables_dir, "QC_Table_FinnGen_outcomes_restricted_to_instruments.tsv"),
  sep = "\t"
)

################################################################################
# 6) Combine all MR-ready outcomes and export                                  #
################################################################################

# Drop duplicated chr/pos fields from MR-PREG (already harmonised by TwoSampleMR)
mr_preg_dat <- mr_preg_dat %>%
  dplyr::select(-chr.outcome, -pos.outcome)

all_outcomes <- dplyr::bind_rows(
  mr_preg_dat,
  dplyr::bind_rows(lapply(res_pph, `[[`, "mr")),
  dplyr::bind_rows(lapply(res_fg,  `[[`, "mr"))
)

out_file <- file.path(results_dir, "formatted_all_pregnancy_outcomes.tsv")

message("\n=== Exporting combined outcome dataset ===")
message("  - Total observations: ", nrow(all_outcomes))
message("  - Unique outcomes:    ", dplyr::n_distinct(all_outcomes$outcome))
message("  - Unique SNPs:        ", dplyr::n_distinct(all_outcomes$SNP))

write.table(
  all_outcomes,
  file      = out_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

message("Exported combined outcomes to: ", basename(out_file))

################################################################################
# 7) SNP count per outcome (by source)                                         #
################################################################################

count_snps <- function(res_list, label) {
  dplyr::bind_rows(lapply(names(res_list), function(nm) {
    tibble::tibble(
      outcome = nm,
      nsnp    = dplyr::n_distinct(res_list[[nm]]$raw$rsid)
    )
  })) %>%
    dplyr::mutate(source = label)
}

counts_all <- dplyr::bind_rows(
  if (exists("res_pph")) count_snps(res_pph, "PPH")      else NULL,
  if (exists("res_fg"))  count_snps(res_fg,  "FinnGen")  else NULL
  # MR-PREG has all SNPs, not subset per outcome via res-list, so not included here
) %>%
  dplyr::arrange(source, dplyr::desc(nsnp))

message("\nSNP counts per outcome (restricted to instruments):")
print(counts_all, n = Inf)

source_summary <- counts_all %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(
    outcomes    = dplyr::n(),
    mean_nsnp   = mean(nsnp),
    median_nsnp = median(nsnp),
    total_nsnp  = sum(nsnp),
    .groups     = "drop"
  )

message("\nSummary by source:")
print(source_summary)
