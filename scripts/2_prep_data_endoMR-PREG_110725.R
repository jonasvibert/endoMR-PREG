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
snps <- c(
  "rs10917151","rs71575922","rs3858429","rs11674184","rs1903068","rs1971256",
  "rs10122243","rs6456259","rs10090060","rs2421985","rs12320196","rs12441483",
  "rs12030576","rs10828249","rs1451383","rs73633307","rs4540228","rs10757277",
  "rs7907732","rs2226158","rs10860864","rs7924571","rs66683298","rs12549438",
  "rs7214750","rs11756073","rs56090796","rs57281976","rs3803042","rs2510770",
  "rs55909142","rs17053711","rs1430787","rs507666","rs2967684","rs7334326",
  "rs2946160","rs10983311","rs1352889","rs2036754","rs6435157"
)

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
  id_col              = "Phenotype"
)

assign("dat_MR_PREG", mr_preg_dat, envir = .GlobalEnv)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Postpartum Hemorrhage Outcomes (Westergaard et al.) --------------------
# ─────────────────────────────────────────────────────────────────────────────
file_list_pph   <- c(
  "Postpartum_hemorrhage_filtered.txt",
  "Postpartum_hemorrhage_due_to_atony_filtered.txt",
  "Postpartum_hemorrhage_due_to_retained_placenta_filtered.txt",
  "Early_bleeding_with_any_outcome_filtered.txt",
  "Antepartum_bleeding_filtered.txt",
  "Early_bleeding_ending_in_live_birth_filtered.txt"
)
common_cols_pph <- c(
  "Chr","PosB38","Marker","rsName","OA","EA",
  "EAfrq","Effect","P","Phet","I2","nCoh"
)

process_pph <- function(fname) {
  base <- tools::file_path_sans_ext(fname)
  raw  <- read_table(
    file.path(data_dir, fname),
    col_names = common_cols_pph
  )[-1, ]
  
  df_pph <- raw %>%
    filter(rsName %in% snps) %>%
    transmute(
      outcome = base,
      rsid    = rsName,
      beta    = log(as.numeric(Effect)),
      SE      = abs(beta) / qnorm(as.numeric(P)/2, lower.tail = FALSE),
      a1      = EA,
      a2      = OA,
      eaf     = as.numeric(EAfrq),
      pval    = as.numeric(P),
      Units   = "logOR",
      Gene    = NA_character_,
      n       = as.integer(nCoh)
    ) %>%
    as.data.frame()
  
  mr_pph <- format_data(
    dat               = df_pph,
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
  mr_pph$outcome    <- base
  mr_pph$id.outcome <- base
  
  assign(paste0("dat_", base), mr_pph, envir = .GlobalEnv)
  
  list(raw = df_pph, mr = mr_pph)
}

res_pph <- lapply(file_list_pph, process_pph)
names(res_pph) <- tools::file_path_sans_ext(file_list_pph)

# QC PPH
qc_pph <- bind_rows(lapply(res_pph, `[[`, "raw")) %>%
  group_by(outcome) %>%
  summarize(
    nsnp      = n(),
    mean_beta = mean(beta, na.rm = TRUE),
    sd_beta   = sd(beta, na.rm = TRUE)
  )
print(qc_pph)


# ─────────────────────────────────────────────────────────────────────────────
# 5) FinnGen R12 Outcomes ---------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
file_list_fg   <- c(
  "finngen_R12_O15_PLAC_PRAEVIA_filtered",
  "finngen_R12_O15_PLAC_DISORD_filtered",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR_filtered",
  "finngen_R12_O15_PREG_ECTOP_filtered",
  "finngen_R12_N14_FEMALEINFERT_filtered"
)
common_cols_fg <- c(
  "#chrom","pos","ref","alt","rsids",
  "pval","beta","sebeta",
  "af_alt","af_alt_cases","af_alt_controls"
)

process_fg <- function(base) {
  raw <- fread(file.path(data_dir, base),
               select = common_cols_fg)
  df  <- raw %>%
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
  mr_fg$outcome    <- base
  mr_fg$id.outcome <- base
  
  assign(paste0("dat_", base), mr_fg, envir = .GlobalEnv)
  
  list(raw = df, mr = mr_fg)
}

res_fg <- lapply(file_list_fg, process_fg)
names(res_fg) <- file_list_fg

# QC FinnGen
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
message("✅ Exported ", nrow(all_outcomes), " rows to ", out_file)
