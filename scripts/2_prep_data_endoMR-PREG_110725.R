#!/usr/bin/env Rscript
# 2_prep_out_endoMR-PREG.R
###############################################################################
# Import, format and harmonize all outcome summary stats for MR:
#  • MR-PREG adverse pregnancy outcomes (ma_out_dat.txt)
#  • Postpartum hemorrhage (Westergaard et al.)
#  • FinnGen R12 adverse pregnancy outcomes
###############################################################################

# 0) Libraries ---------------------------------------------------------------
library(data.table)   # fread(), fwrite()
library(dplyr)        # %>%
library(readr)        # read_table()
library(TwoSampleMR)  # format_data()
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# 1) Paths ------------------------------------------------------------------
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# 2) SNP list ---------------------------------------------------------------
snps <- scan(file.path(data_dir, "rsid.txt"), what = character())

# ─────────────────────────────────────────────────────────────────────────────
# 3) MR-PREG adverse pregnancy outcomes ------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
f_in       <- file.path(data_dir, "ma_out_dat.txt")
raw_lines  <- readLines(f_in)
header     <- raw_lines[1]
data_clean <- raw_lines[-1][!grepl("\\.png$", raw_lines[-1])]
txt_clean  <- c(header, data_clean)

mr_preg_df <- read.table(
  text            = txt_clean,
  header          = TRUE,
  stringsAsFactors = FALSE
)

# Format for TwoSampleMR
mr_preg_dat <- format_data(
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

# ─────────────────────────────────────────────────────────────────────────────
# 4) Postpartum Hemorrhage Outcomes (Westergaard et al.) --------------------
# ─────────────────────────────────────────────────────────────────────────────

# List of PPH outcome files (without extension, like FinnGen style)
file_list_pph <- c(
  "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta",
  "Early_bleeding_with_any_outcome",
  "Antepartum_bleeding",
  "Early_bleeding_ending_in_live_birth"
)

# Expected columns in the Westergaard files
common_cols_pph <- c(
  "Chr","PosB38","Marker","rsName","OA","EA",
  "EAfrq","Effect","P","Phet","I2","nCoh"
)

# Function to read, filter on SNPs, and format (analogous to process_fg)
process_pph <- function(base) {
  raw <- readr::read_table(
    file.path(data_dir, paste0(base, ".txt")),
    col_names = common_cols_pph,
    show_col_types = FALSE
  )[-1, ]  # drop repeated header row
  
  # Filter for selected SNPs
  df <- raw %>%
    dplyr::filter(rsName %in% snps) %>%
    dplyr::transmute(
      outcome       = base,
      rsid          = rsName,
      chr           = as.character(Chr),
      pos           = as.integer(PosB38),
      effect_allele = EA,
      other_allele  = OA,
      eaf           = as.numeric(EAfrq),
      beta          = as.numeric(Effect),
      pval          = as.numeric(P)
    ) %>%
    dplyr::mutate(
      se = abs(beta) / qnorm(pval/2, lower.tail = FALSE)
    ) %>%
    as.data.frame()
  
  # Format for MR
  mr_pph <- format_data(
    dat               = df,
    type              = "outcome",
    snp_col           = "rsid",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "chr",
    pos_col           = "pos",
    phenotype_col     = "outcome",
    id_col            = "outcome"
  )
  
  # Ensure consistent naming (as in FinnGen)
  mr_pph$outcome    <- base
  mr_pph$id.outcome <- base
  
  # Save in global environment (as in FinnGen)
  assign(paste0("dat_", base), mr_pph, envir = .GlobalEnv)
  
  list(raw = df, mr = mr_pph)
}

# Run for all PPH outcomes
res_pph <- lapply(file_list_pph, process_pph)
names(res_pph) <- file_list_pph

# QC summary 
qc_pph <- dplyr::bind_rows(lapply(res_pph, `[[`, "raw")) %>%
  dplyr::group_by(outcome) %>%
  dplyr::summarize(
    nsnp      = dplyr::n(),
    mean_beta = mean(beta, na.rm = TRUE),
    sd_beta   = sd(beta, na.rm = TRUE),
    .groups   = "drop"
  )
print(qc_pph)


# ─────────────────────────────────────────────────────────────────────────────
# 5) FinnGen R12 Outcomes -----------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# List of raw (unfiltered) FinnGen files
file_list_fg <- c(
  "finngen_R12_O15_PLAC_PRAEVIA",
  "finngen_R12_O15_PLAC_DISORD",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",
  "finngen_R12_O15_PREG_ECTOP",
  "finngen_R12_N14_FEMALEINFERT"
)

# Expected columns in FinnGen files
common_cols_fg <- c(
  "#chrom","pos","ref","alt","rsids",
  "pval","beta","sebeta",
  "af_alt","af_alt_cases","af_alt_controls"
)

# Function to read, filter on SNPs, and format
process_fg <- function(base) {
  raw <- fread(
    file.path(data_dir, base),
    select = common_cols_fg
  )
  
  # Filter for selected SNPs
  df <- raw %>%
    filter(rsids %in% snps) %>%
    transmute(
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
  
  # Format for MR
  mr_fg <- format_data(
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
  
  # Ensure consistent naming
  mr_fg$outcome    <- base
  mr_fg$id.outcome <- base
  
  # Save in global environment
  assign(paste0("dat_", base), mr_fg, envir = .GlobalEnv)
  
  list(raw = df, mr = mr_fg)
}

# Run for all FinnGen outcomes
res_fg <- lapply(file_list_fg, process_fg)
names(res_fg) <- file_list_fg

# QC summary
qc_fg <- bind_rows(lapply(res_fg, `[[`, "raw")) %>%
  group_by(outcome) %>%
  summarize(
    nsnp      = n(),
    mean_beta = mean(beta, na.rm = TRUE),
    sd_beta   = sd(beta, na.rm = TRUE)
  )
print(qc_fg)


# ─────────────────────────────────────────────────────────────────────────────
# 6) Combine all MR-ready outcomes & Export -------------------------------
# ─────────────────────────────────────────────────────────────────────────────
mr_preg_dat <- mr_preg_dat %>%
  select(-chr.outcome, -pos.outcome)

all_outcomes <- bind_rows(
  mr_preg_dat,
  bind_rows(lapply(res_pph, `[[`, "mr")),
  bind_rows(lapply(res_fg,  `[[`, "mr"))
)

out_file <- file.path(results_dir, "formatted_all_pregnancy_outcomes.tsv")
write.table(
  all_outcomes,
  file      = out_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)


# ─────────────────────────────────────────────────────────────────────────────
# 7) outcome summary table for MR datasets -------------------------------
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(stringr)

## ---- 1) MR-PREG ----
case_col <- names(mr_preg_dat)[str_detect(names(mr_preg_dat), "^n(case|_case)\\.?outcome$")] %>% dplyr::first()
ctrl_col <- names(mr_preg_dat)[str_detect(names(mr_preg_dat), "^n(control|_control)\\.?outcome$")] %>% dplyr::first()

mrpreg_summary <- mr_preg_dat %>%
  group_by(outcome) %>%
  summarise(
    nsnp         = n_distinct(SNP),
    N_total      = suppressWarnings(as.integer(median(samplesize.outcome, na.rm = TRUE))),
    N_total_min  = suppressWarnings(as.integer(min(samplesize.outcome, na.rm = TRUE))),
    N_total_max  = suppressWarnings(as.integer(max(samplesize.outcome, na.rm = TRUE))),
    cases        = if (!is.null(case_col)) suppressWarnings(as.integer(median(.data[[case_col]], na.rm = TRUE))) else NA_integer_,
    controls     = if (!is.null(ctrl_col)) suppressWarnings(as.integer(median(.data[[ctrl_col]], na.rm = TRUE))) else NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(source = "MR-PREG") %>%
  relocate(source)

## ---- 2) Westergaard PPH ----
pph_summary <- bind_rows(lapply(res_pph, `[[`, "raw")) %>%
  group_by(outcome) %>%
  summarise(
    nsnp    = n_distinct(rsid),
    cohorts = suppressWarnings(as.integer(median(n, na.rm = TRUE))),
    N_total  = NA_integer_,
    cases    = NA_integer_,
    controls = NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(source = "Westergaard (PPH)") %>%
  relocate(source)

## ---- 3) FinnGen ----
finngen_summary <- bind_rows(lapply(res_fg, `[[`, "raw")) %>%
  group_by(outcome) %>%
  summarise(
    nsnp = n_distinct(rsid),
    N_total  = NA_integer_,
    cases    = NA_integer_,
    controls = NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(source = "FinnGen R12") %>%
  relocate(source)

## ---- 4) Combine all ----
outcome_summary <- bind_rows(
  mrpreg_summary %>% select(source, outcome, nsnp, N_total, cases, controls),
  pph_summary %>% select(source, outcome, nsnp, N_total, cases, controls),
  finngen_summary %>% select(source, outcome, nsnp, N_total, cases, controls)
) %>%
  arrange(source, outcome)

## ---- 5) Print directly in console ----
print(outcome_summary, n = nrow(outcome_summary))

## ---- 6) Print ranges for manuscript text ----
mrpreg_range <- mrpreg_summary %>%
  summarise(
    minN = min(N_total_min, na.rm = TRUE),
    maxN = max(N_total_max, na.rm = TRUE)
  ) %>%
  mutate(
    range_text = paste0("N = ", scales::comma(minN), " to N = ", scales::comma(maxN))
  ) %>%
  pull(range_text)

cat("\n--- Section 2.1 MR-PREG ---\n")
cat("Total sample sizes ranged from ", mrpreg_range, ".\n", sep = "")

cat("\n--- Section 2.2 FinnGen ---\n")
cat("Sample sizes to be added from FinnGen metadata.\n")

cat("\n--- Section 2.3 Westergaard PPH ---\n")
cat("Provide subtype-specific totals or use overall study totals (overall PPH ≈ 175,000).\n")

