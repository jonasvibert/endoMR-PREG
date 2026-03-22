################################################################################
# Script: 01_select_instruments_endoMR-PREG.R
# Project: endoMR-PREG
#
# Purpose:
#   Select, map and clump SNP instruments for endometriosis using the
#   Rahmioglu et al. GWAS (Supp. Table 33 – "Supp33.Top10K-SNPs")
#   extracted from: NIHMS1873427-Supplementary_Materials.xlsx
################################################################################

### 0) Setup ###################################################################

# Clear workspace and set numeric display
rm(list = ls())
options(digits = 10)

# Packages ---------------------------------------------------------------------
required_pkgs <- c(
  "TwoSampleMR",
  "MRInstruments",
  "knitr",
  "ggplot2",
  "dplyr",
  "tidyr",
  "openxlsx",
  "data.table",
  "here",
  "readxl"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

# Directories ------------------------------------------------------------------
data_dir    <- here::here("data")
results_dir <- here::here("results")
scripts_dir <- here::here("scripts")
plots_dir   <- file.path(results_dir, "plots instruments selection")
tables_dir  <- file.path(results_dir, "tables instruments selection")

dir.create(data_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(scripts_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)

################################################################################
# 1) Load Endometriosis GWAS - Rahmioglu data                                  #
################################################################################

# Load the sheet with the top 10,000 SNPs (Supp. Table 33, including 23andMe)
expdat_raw <- readxl::read_excel(
  here::here("data", "EXPOSURE_RAHMIGLU", "NIHMS1873427-Supplementary_Materials.xlsx"),
  sheet = "Supp33.Top10K-SNPs",
  skip  = 1
)

print(colnames(expdat_raw))
head(expdat_raw)

################################################################################
# 2) Format exposure data for TwoSampleMR                                      #
################################################################################

# Clean column names
names(expdat_raw) <- gsub(" ", ".", names(expdat_raw))
names(expdat_raw) <- gsub("\\(hg.19\\)", "Position.hg19", names(expdat_raw))

expdat_raw$phenotype <- "Endometriosis (Rahmioglu)"

expdat <- TwoSampleMR::format_data(
  expdat_raw,
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "Effect",
  se_col            = "Standard.Error",
  eaf_col           = "Effective.Allele.Frequency",
  effect_allele_col = "Effective.Allele",
  other_allele_col  = "Non-effective.Allele",
  pval_col          = "P-value",
  samplesize_col    = "Sample.Size",
  chr_col           = "Chromosome",
  pos_col           = "Position.hg19",
  phenotype_col     = "phenotype",
  id_col            = "id.exposure"
)

# Clean possible suffix artefacts
colnames(expdat) <- gsub("\\.id\\.exposure", ".exposure", colnames(expdat))
colnames(expdat)

################################################################################
# 3) Map chr:position to rsIDs using 1000G (bim file)                          #
################################################################################

# PLINK / 1000G paths ----------------------------------------------------------
plink_dir  <- here::here("plink_mac_20241022")
plink_bin  <- file.path(plink_dir, "plink")
bfile_root <- file.path(plink_dir, "24088632", "1000G.EUR.QC")
bim_file   <- paste0(bfile_root, ".bim")

# Load .bim
bim_data <- data.table::fread(bim_file, header = FALSE)
data.table::setnames(
  bim_data,
  c("CHR", "SNP_bim", "CM", "BP", "A1", "A2")
)
bim_data[, chrpos := paste0("chr", CHR, ":", BP)]

# Create chrpos in exposure data
expdat$chrpos <- expdat$SNP

# Remove unnecessary duplicate columns before merging
expdat_clean <- expdat %>%
  dplyr::select(-dplyr::matches("^id\\.exposure$"))

# Merge on chrpos to map rsIDs
expdat_mapped <- merge(
  expdat_clean,
  bim_data[, c("chrpos", "SNP_bim")],
  by    = "chrpos",
  all.x = TRUE
)

# Replace SNP column by rsID and keep only mapped SNPs
expdat_mapped$SNP <- expdat_mapped$SNP_bim

expdat_ready <- expdat_mapped %>%
  dplyr::filter(!is.na(SNP))

cat("Number of SNPs mapped with an rsID:", nrow(expdat_ready), "\n")
head(expdat_ready$SNP)

################################################################################
# 4) LD-based clumping using PLINK and 1000G                                   #
################################################################################

clumped <- TwoSampleMR::clump_data(
  expdat_ready,
  clump_kb  = 10000,
  clump_r2  = 0.001,
  clump_p1  = 5e-8,
  clump_p2  = 1,
  bfile     = bfile_root,
  plink_bin = plink_bin
)

clumped_snps <- clumped$SNP
final_instruments <- expdat_ready[expdat_ready$SNP %in% clumped_snps, ]

nrow(clumped)
head(clumped)

# Save formatted exposure dataset (raw instruments before clumping)
write.table(
  expdat_ready,
  file      = file.path(data_dir, "Exposure_Endometriosis_Rahmioglu_endo_dat.txt"),
  quote     = FALSE,
  row.names = FALSE,
  sep       = "\t"
)

# Save unique rsID list
rsid <- unique(expdat_ready$SNP)
write.table(
  rsid,
  file      = file.path(data_dir, "Exposure_Endometriosis_Rahmioglu_rsid.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

# Save clumped SNPs list
data.table::fwrite(
  clumped,
  file.path(results_dir, "Exposure_Endometriosis_Rahmioglu_clumped_snps.tsv"),
  sep = "\t"
)

################################################################################
# 5) Instrument strength: F-statistics and variance explained                  #
################################################################################

# 5.1 Retrieve SNP statistics for clumped SNPs --------------------------------
stats <- expdat_ready %>%
  dplyr::select(
    SNP,
    beta.exposure,
    se.exposure,
    pval.exposure,
    eaf.exposure,
    samplesize.exposure
  )

clumped2 <- clumped %>%
  dplyr::left_join(stats, by = "SNP")

names(clumped2)

clumped2 %>%
  dplyr::select(
    SNP,
    beta.exposure = beta.exposure.y,
    se.exposure   = se.exposure.y,
    eaf.exposure  = eaf.exposure.y,
    samplesize.exposure = samplesize.exposure.y
  ) %>%
  head()

# 5.2 Simple F-statistic (beta / SE)^2 ----------------------------------------
clumped2 <- clumped2 %>%
  dplyr::mutate(
    beta.exposure       = dplyr::coalesce(beta.exposure.y,       beta.exposure.x),
    se.exposure         = dplyr::coalesce(se.exposure.y,         se.exposure.x),
    eaf.exposure        = dplyr::coalesce(eaf.exposure.y,        eaf.exposure.x),
    samplesize.exposure = dplyr::coalesce(samplesize.exposure.y, samplesize.exposure.x),
    F_simple            = (beta.exposure / se.exposure)^2
  )

clumped2 %>%
  dplyr::select(SNP, beta.exposure, se.exposure, F_simple) %>%
  head()

# 5.3 Approximate per-SNP R² and exact F-statistic ----------------------------
clumped2 <- clumped2 %>%
  dplyr::mutate(
    R2 = (2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2) /
      (2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2 +
         se.exposure^2 * samplesize.exposure),
    F_exact = R2 * (samplesize.exposure - 2) / (1 - R2)
  )

head(
  clumped2[, c("SNP", "beta.exposure", "se.exposure", "F_simple", "R2", "F_exact")],
  10
)

cat("Mean F_simple =", round(mean(clumped2$F_simple, na.rm = TRUE), 2), "\n")
cat("Mean F_exact  =", round(mean(clumped2$F_exact,  na.rm = TRUE), 2), "\n")

################################################################################
# 6) Variance explained by selected SNPs                                       #
################################################################################

clumped2 <- clumped2 %>%
  dplyr::mutate(
    R2_i = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
    F_i  = R2_i * (samplesize.exposure - 2) / (1 - R2_i)
  )

total_R2 <- sum(clumped2$R2_i, na.rm = TRUE)

k <- nrow(clumped2)
N <- as.integer(median(clumped2$samplesize.exposure, na.rm = TRUE))

F_multi <- (total_R2 / k) * ((N - k - 1) / (1 - total_R2))
mean_F  <- mean(clumped2$F_i, na.rm = TRUE)

cat(sprintf("Total variance explained (R²)        = %.6g\n", total_R2))
cat(sprintf("Multi-SNP F-statistic               = %.3f\n", F_multi))
cat(sprintf("Mean per-SNP F-statistic            = %.3f\n", mean_F))

################################################################################
# 7) Visualisations (saved with explicit figure names)                         #
################################################################################

# A) Bar plot: SNP count per chromosome ---------------------------------------
chr_summary <- expdat_ready %>%
  dplyr::filter(SNP %in% clumped$SNP) %>%
  dplyr::mutate(chr = sub(":.*", "", chrpos)) %>%
  dplyr::mutate(chr = factor(chr, levels = paste0("chr", 1:22))) %>%
  dplyr::count(chr) %>%
  dplyr::filter(n > 0)

p_chr <- ggplot(chr_summary, aes(x = chr, y = n)) +
  geom_col() +
  labs(
    title = "Endometriosis instruments: SNPs per chromosome",
    x     = "Chromosome",
    y     = "Number of SNPs"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plots_dir, "Figure_chr_distribution_endometriosis_instruments.png"),
  plot     = p_chr,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# B) Histogram of effect allele frequency (EAF) --------------------------------
p_eaf <- ggplot(
  clumped2 %>% dplyr::filter(SNP %in% clumped$SNP),
  aes(x = eaf.exposure)
) +
  geom_histogram(bins = 20) +
  labs(
    title = "Endometriosis instruments: distribution of effect allele frequency",
    x     = "Effect allele frequency",
    y     = "Number of SNPs"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plots_dir, "Figure_EAF_distribution_endometriosis_instruments.png"),
  plot     = p_eaf,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# C) Forest plot of β estimates with 95% CI ------------------------------------
forest_df <- clumped2 %>%
  dplyr::mutate(
    lower = beta.exposure - 1.96 * se.exposure,
    upper = beta.exposure + 1.96 * se.exposure
  ) %>%
  dplyr::arrange(beta.exposure) %>%
  dplyr::mutate(SNP = factor(SNP, levels = SNP))

p_forest <- ggplot(forest_df, aes(x = beta.exposure, y = SNP)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  labs(
    title = "Endometriosis instruments: SNP effect estimates",
    x     = expression(beta~"(95% CI)"),
    y     = "SNP"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plots_dir, "Figure_forest_endometriosis_instruments.png"),
  plot     = p_forest,
  width    = 7,
  height   = 6,
  dpi      = 300
)

# D) Volcano plot: β vs –log10(p-value) ---------------------------------------
volc_df <- clumped2 %>%
  dplyr::mutate(logp = -log10(pval.exposure.x))

p_volcano <- ggplot(volc_df, aes(x = beta.exposure, y = logp)) +
  geom_point() +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  labs(
    title = "Endometriosis instruments: volcano plot",
    x     = expression(beta),
    y     = expression(-log[10]~"(p-value)")
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plots_dir, "Figure_volcano_endometriosis_instruments.png"),
  plot     = p_volcano,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# E) Scatterplot: F_simple vs R² ----------------------------------------------
p_F_R2 <- ggplot(clumped2, aes(x = R2, y = F_simple)) +
  geom_point() +
  labs(
    title = "Endometriosis instruments: F-statistic vs R²",
    x     = expression(R^2),
    y     = expression(F["simple"])
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plots_dir, "Figure_Fstat_vs_R2_endometriosis_instruments.png"),
  plot     = p_F_R2,
  width    = 7,
  height   = 5,
  dpi      = 300
)

################################################################################
# 8) Tables                                          #
################################################################################

# Helper: coalesce over candidate columns -------------------------------------
coalesce_into <- function(df, new_name, candidates) {
  df %>%
    dplyr::mutate(
      "{new_name}" := dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(candidates)))
    )
}

# 1) Build unified columns safely ---------------------------------------------
clumped2 <- clumped2 %>%
  coalesce_into("beta.exposure",       c("beta.exposure",       "beta.exposure.y",       "beta.exposure.x")) %>%
  coalesce_into("se.exposure",         c("se.exposure",         "se.exposure.y",         "se.exposure.x")) %>%
  coalesce_into("eaf.exposure",        c("eaf.exposure",        "eaf.exposure.y",        "eaf.exposure.x")) %>%
  coalesce_into("pval.exposure",       c("pval.exposure",       "pval.exposure.y",       "pval.exposure.x")) %>%
  coalesce_into("samplesize.exposure", c("samplesize.exposure", "samplesize.exposure.y", "samplesize.exposure.x"))

# 2) Ensure per-SNP R² and F_i are present ------------------------------------
if (!all(c("R2_i", "F_i") %in% names(clumped2))) {
  clumped2 <- clumped2 %>%
    dplyr::mutate(
      R2_i = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
      F_i  = R2_i * (samplesize.exposure - 2) / (1 - R2_i)
    )
}

# 3) Build Supplementary Table S1 ---------------------------------------------
supp_table <- clumped2 %>%
  dplyr::select(
    SNP,
    Beta          = beta.exposure,
    SE            = se.exposure,
    EAF           = eaf.exposure,
    `P-value`     = pval.exposure,
    `R²`          = R2_i,
    `F-statistic` = F_i
  ) %>%
  dplyr::mutate(
    `P-value` = formatC(`P-value`, format = "e", digits = 2)
  )

# 4) Export to Excel (explicit supplementary table name) ----------------------
openxlsx::write.xlsx(
  supp_table,
  file.path(tables_dir, "Endometriosis_Instruments.xlsx"),
  rowNames = FALSE
)
saveRDS(clumped2, file.path(results_dir, "clumped2_with_stats.rds"))
message("Saved clumped2 object to: ", file.path(results_dir, "clumped2_with_stats.rds"))

cat("Saved: ", file.path(tables_dir, "Endometriosis_Instruments.xlsx"), "\n")
