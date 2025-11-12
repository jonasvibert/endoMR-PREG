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

# Load progress bar library
if (!requireNamespace("progress", quietly = TRUE)) {
  message("Installing progress package for progress bars...")
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  install.packages("progress")
}
library(progress)

# Load forest plot libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  install.packages("ggplot2")
}
if (!requireNamespace("forestplot", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  install.packages("forestplot")
}
if (!requireNamespace("metafor", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  install.packages("metafor")
}
library(ggplot2)
library(forestplot)
library(metafor)

# 1) Paths ------------------------------------------------------------------
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# 2) SNP list ---------------------------------------------------------------
snps_file <- file.path(results_dir, "endo_instruments_41.txt")

if (file.exists(snps_file)) {
  snps <- scan(snps_file,
               what = character(),
               sep = "\n",
               quiet = TRUE,
               strip.white = TRUE)
  snps <- unique(snps[nzchar(snps)])                 # drop blanks + dedupe
  snps <- snps[!tolower(snps) %in% c("snp", "rsid")] # drop any header
} else if (file.exists(file.path(results_dir, "endometriosis_clumped_snps.tsv"))) {
  # fallback: read from the TSV saved by the selection script
  tmp <- data.table::fread(file.path(results_dir, "endometriosis_clumped_snps.tsv"),
                           select = "SNP", showProgress = FALSE)
  snps <- unique(tmp$SNP)
} else {
  stop("Could not find the 41-SNP list. Run the selection script first, ",
       "or place 'endo_instruments_41.txt' in ", results_dir)
}

message("Loaded ", length(snps), " instruments (expected ~41).")
message("Starting data preparation for pregnancy outcomes...")
# head(snps)

# ─────────────────────────────────────────────────────────────────────────────
# 3) MR-PREG adverse pregnancy outcomes ------------------------------------
#    Units: Direct beta coefficients (log odds ratios) with standard errors
# ─────────────────────────────────────────────────────────────────────────────

message("\n=== Processing MR-PREG adverse pregnancy outcomes ===")
message("Data source: MR-PREG consortium")
message("Expected units: Beta (log OR) and SE")
f_in       <- file.path(data_dir, "OUTCOME_MR-PREG", "ma_out_dat.txt")
message("Reading MR-PREG data from: ", basename(f_in))
raw_lines  <- readLines(f_in)
header     <- raw_lines[1]
data_clean <- raw_lines[-1][!grepl("\\.png$", raw_lines[-1])]
txt_clean  <- c(header, data_clean)
message("Cleaned ", length(data_clean), " data lines (removed ", length(raw_lines) - 1 - length(data_clean), " non-data lines)")

mr_preg_df <- read.table(
  text            = txt_clean,
  header          = TRUE,
  stringsAsFactors = FALSE
)

# Diagnostic messages for MR-PREG data
message("MR-PREG data dimensions: ", nrow(mr_preg_df), " rows x ", ncol(mr_preg_df), " columns")
message("Available outcomes: ", n_distinct(mr_preg_df$Phenotype))
message("SNPs in data: ", n_distinct(mr_preg_df$SNP))
message("Mean beta (log OR): ", sprintf("%.4f", mean(mr_preg_df$beta, na.rm = TRUE)))
message("Mean SE: ", sprintf("%.4f", mean(mr_preg_df$se, na.rm = TRUE)))
message("Units confirmed: Beta coefficients (log OR) with SE")

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
message("MR-PREG processing complete: ", nrow(mr_preg_dat), " observations formatted")

# ─────────────────────────────────────────────────────────────────────────────
# 4) Postpartum Hemorrhage Outcomes (Westergaard et al.)
#    Units: Effect sizes (OR or beta) auto-detected and converted to log OR
# ─────────────────────────────────────────────────────────────────────────────

message("\n=== Processing Postpartum Hemorrhage outcomes ===")
message("Data source: Westergaard et al. Nature Genetics 2024")
message("Auto-detection: OR (≈1.0) → log OR, or direct beta if already log scale")

# List of raw PPH files (without extension)
file_list_pph <- c(
  "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta",
  "Early_bleeding_with_any_outcome",
  "Antepartum_bleeding",
  "Early_bleeding_ending_in_live_birth"
)

# Expected columns in PPH files
common_cols_pph <- c(
  "Chr","PosB38","Marker","rsName","OA","EA",
  "EAfrq","Effect","P","Phet","I2","nCoh"
)

# Helper: robust SE from two-sided p-value
.se_from_p <- function(beta, pval) {
  p2 <- pmin(pmax(pval/2, .Machine$double.xmin), 1 - 1e-16)
  abs(beta) / qnorm(p2, lower.tail = FALSE)
}

# Helper: auto-detect if "Effect" is an OR (~1) or already beta (~0)
.to_beta <- function(effect_vec, outcome_name = "") {
  med <- suppressWarnings(median(effect_vec, na.rm = TRUE))
  if (is.finite(med) && med > 0.5 && med < 1.5) {
    message("  ", outcome_name, ": Auto-detected OR format (median = ", sprintf("%.3f", med), ") → converting to log OR")
    log(effect_vec)  # OR → beta
  } else {
    message("  ", outcome_name, ": Auto-detected beta format (median = ", sprintf("%.3f", med), ") → using as-is")
    effect_vec      # already beta
  }
}

# Function to read, filter on SNPs, and format (PPH)
process_pph <- function(base) {
  path <- file.path(data_dir, "OUTCOME_PPH", paste0(base, ".txt"))
  if (!file.exists(path)) stop("File not found: ", path)
  
  # Read full file; assign column names explicitly
  raw <- data.table::fread(
    path,
    col.names     = common_cols_pph,
    data.table    = FALSE,
    showProgress  = FALSE
  )
  
  # Drop any repeated header row 
  if (nrow(raw) > 0 && all(tolower(as.character(raw[1, ])) == tolower(common_cols_pph))) {
    raw <- raw[-1, , drop = FALSE]
  }
  
  # Filter for selected SNPs
  df <- raw %>%
    filter(rsName %in% snps) %>%
    transmute(
      outcome = base,
      rsid    = rsName,
      chr     = gsub("^chr", "", as.character(Chr)),
      pos     = as.integer(PosB38),
      a1      = EA,              # effect allele (EA)
      a2      = OA,              # other allele (OA)
      eaf     = suppressWarnings(as.numeric(EAfrq)),
      beta    = .to_beta(suppressWarnings(as.numeric(Effect)), base),
      pval    = suppressWarnings(as.numeric(P))
    ) %>%
    mutate(
      SE      = .se_from_p(beta, pval),
      Units   = "logOR",
      Gene    = NA_character_,
      n       = NA_integer_
    ) %>%
    as.data.frame()
  
  # Format for MR
  mr_pph <- format_data(
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
  
  # Ensure consistent naming
  mr_pph$outcome    <- base
  mr_pph$id.outcome <- base
  
  # Save in global environment
  assign(paste0("dat_", base), mr_pph, envir = .GlobalEnv)
  
  list(raw = df, mr = mr_pph)
}

# Run for all PPH outcomes with progress bar
message("Processing ", length(file_list_pph), " PPH outcomes...")
pb_pph <- progress_bar$new(
  format = "PPH [:bar] :percent :current/:total ETA: :eta",
  total = length(file_list_pph), 
  clear = FALSE
)

res_pph <- list()
for (i in seq_along(file_list_pph)) {
  pb_pph$tick()
  res_pph[[i]] <- process_pph(file_list_pph[i])
}
names(res_pph) <- file_list_pph

# QC summary
qc_pph <- bind_rows(lapply(res_pph, `[[`, "raw")) %>%
  group_by(outcome) %>%
  summarize(
    nsnp      = n(),
    mean_beta = mean(beta, na.rm = TRUE),
    sd_beta   = sd(beta, na.rm = TRUE),
    .groups   = "drop"
  )
print(qc_pph)



# ─────────────────────────────────────────────────────────────────────────────
# 5) FinnGen R12 Outcomes -----------------------------------------------------
#    Units: Direct beta coefficients (log odds ratios) with standard errors
# ─────────────────────────────────────────────────────────────────────────────

message("\n=== Processing FinnGen R12 adverse pregnancy outcomes ===")
message("Data source: FinnGen R12 release")
message("Units: Direct beta (log OR) and sebeta (SE)")

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
    file.path(data_dir,"OUTCOME_FINNGEN", base),
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

# Run for all FinnGen outcomes with progress bar
message("Processing ", length(file_list_fg), " FinnGen outcomes...")
pb_fg <- progress_bar$new(
  format = "FinnGen [:bar] :percent :current/:total ETA: :eta",
  total = length(file_list_fg), 
  clear = FALSE
)

res_fg <- list()
for (i in seq_along(file_list_fg)) {
  pb_fg$tick()
  res_fg[[i]] <- process_fg(file_list_fg[i])
}
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
message("\n=== Exporting combined dataset ===")
message("Total observations: ", nrow(all_outcomes))
message("Unique outcomes: ", n_distinct(all_outcomes$outcome))
message("Unique SNPs: ", n_distinct(all_outcomes$SNP))

write.table(
  all_outcomes,
  file      = out_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("Exported to: ", basename(out_file))

# ─────────────────────────────────────────────────────────────────────────────
# 7) SNP count per outcome_dat -------------------------------
# ─────────────────────────────────────────────────────────────────────────────
count_snps <- function(res_list, label) {
  bind_rows(lapply(names(res_list), function(nm) {
    tibble(outcome = nm, nsnp = n_distinct(res_list[[nm]]$raw$rsid))
  })) %>% mutate(source = label)
}

counts_all <- bind_rows(
  if (exists("res_pph"))    count_snps(res_pph,    "PPH")    else NULL,
  if (exists("res_fg"))     count_snps(res_fg,     "FinnGen")else NULL,
  if (exists("res_mrpreg")) count_snps(res_mrpreg, "MR-PREG")else NULL
) %>% arrange(source, desc(nsnp))

print(counts_all, n = Inf)

# Optional summary by source
source_summary <- counts_all %>%
  group_by(source) %>%
  summarise(outcomes = n(),
            mean_nsnp = mean(nsnp),
            median_nsnp = median(nsnp),
            total_nsnp = sum(nsnp),
            .groups = "drop")

print(source_summary)

