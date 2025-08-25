############################################################################
#                                                                          #
#   1/ Select SNPs & clump from endometriosis GWAS - Rahmioglu dataset     #
#                                                                          #
############################################################################

############################################################################
#                                SETUP                                     #
############################################################################

## Clear the workspace
  rm(list = ls())

## Set digits for numeric display
  options(digits = 10)

## Load necessary packages
  #install.packages("remotes")
  #remotes::install_github("MRCIEU/TwoSampleMR")
  # remotes::install_github("MRCIEU/MRInstruments")
  
  library(TwoSampleMR)
  library(MRInstruments)
  library(knitr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(data.table)

## Define project directories
  if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
  library(here)
  data_dir    <- here("data")
  results_dir <- here("results")
  scripts_dir <- here("scripts")
  plots_dir   <- here("plots")

## Create directories 
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(scripts_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

############################################################################
#                   Load Endometriosis GWAS - Rahmioglu data               #
############################################################################
  
## Load top 10,000 SNPs from supplementary Table 33 (including 23andMe)
  readxl::excel_sheets(here("data", "EXPOSURE_RAHMIGLU", "NIHMS1873427-Supplementary_Materials.xlsx"))
  
  # Load the sheet with the top 10,000 SNPs
  expdat_raw <- readxl::read_excel(
    here("data", "EXPOSURE_RAHMIGLU", "NIHMS1873427-Supplementary_Materials.xlsx"),
    sheet = "Supp33.Top10K-SNPs",
    skip = 1
  )
  
  # Check column names and view first rows
  print(colnames(expdat_raw))
  head(expdat_raw)
  
## Format data for TwoSampleMR
  
  # Rename columns to avoid issues with spaces and special characters
  names(expdat_raw) <- gsub(" ", ".", names(expdat_raw))
  names(expdat_raw) <- gsub("\\(hg.19\\)", "Position.hg19", names(expdat_raw))
  
  expdat_raw$phenotype <- "Endometriosis (Rahmioglu)"
  expdat <- format_data(
    expdat_raw,
    type               = "exposure",                 
    snp_col            = "SNP",
    beta_col           = "Effect",
    se_col             = "Standard.Error",
    eaf_col            = "Effective.Allele.Frequency",
    effect_allele_col  = "Effective.Allele",
    other_allele_col   = "Non-effective.Allele",
    pval_col           = "P-value",
    samplesize_col     = "Sample.Size",
    chr_col            = "Chromosome",              
    pos_col            = "Position.hg19",        
    phenotype_col      = "phenotype",               
    id_col             = "id.exposure"
  )

############################################################################
#                      LD-Based Clumping Using PLINK and 1000G            #
############################################################################
colnames(expdat) <- gsub("\\.id\\.exposure", ".exposure", colnames(expdat))
colnames(expdat)
# path to the PLINK directory
plink_dir   <- here("plink_mac_20241022")
# PLINK executable
plink_path  <- here("plink_mac_20241022","plink")
# prefix for the PLINK files (without .bed/.bim/.fam extensions)
bfile_path  <- here("plink_mac_20241022","24088632")
# now derive each file on the fly
bed_file <- paste0(bfile_path, ".bed")
bim_file <- paste0(bfile_path, ".bim")
fam_file <- paste0(bfile_path, ".fam")


## Map SNP chr:position to rsIDs 
# Load .bim
bim_file <- file.path(bfile_path, "1000G.EUR.QC.bim")
bim_data <- data.table::fread(bim_file, header = FALSE)
setnames(bim_data, c("CHR", "SNP_bim", "CM", "BP", "A1", "A2"))
bim_data[, chrpos := paste0("chr", CHR, ":", BP)]

# Create chrpos in expdat.
  expdat$chrpos <- expdat$SNP

# Remove unnecessary duplicate columns before merging
  expdat_clean <- expdat %>% select(-matches("^id\\.exposure$"))
# Merge
  expdat_mapped <- merge(expdat_clean, bim_data[, c("chrpos", "SNP_bim")], by = "chrpos", all.x = TRUE)

# Replace the SNP column with the rsID  
  expdat_mapped$SNP <- expdat_mapped$SNP_bim

# Filter only SNPs with an rsID
  expdat_ready <- expdat_mapped %>% filter(!is.na(SNP))

# Check how many SNPs are mapped
  cat("Number of SNPs mapped with an rsID:", nrow(expdat_ready), "\n")
  head(expdat_ready$SNP)

## Perform LD-based clumping using PLINK and 1000G reference
  clumped <- clump_data(
    expdat_ready,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = 5e-8,
    clump_p2 = 1,
    bfile = file.path(plink_dir, "24088632", "1000G.EUR.QC"),
    plink_bin = plink_path
  )

  clumped_snps <- clumped$SNP

# Filter final dataset
  final_instruments <- expdat_ready[expdat_ready$SNP %in% clumped_snps, ]

# Check the number of SNPs and preview the first few rows of the clumped dataset
  nrow(clumped)
  head(clumped)

## Save the formatted exposure dataset 
  write.table(expdat_ready, file = file.path(data_dir, "endo_dat.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

#Save and unique list of rsIDs
  rsid <- unique(expdat_ready$SNP)
  write.table(rsid, file = file.path(data_dir, "rsid.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

#Save table with list of rsIDs
  data.table::fwrite(
    clumped,
    file.path(results_dir, "endometriosis_clumped_snps.tsv"),
    sep = "\t"
  )

  
  ############################################################################
  #                      LD-Based Clumping Using PLINK and UKbiobank        #
  ############################################################################
  # 2.1. Read the file (no header)
  # plink_path <- "/Users/jonasvibert/Desktop/Stata/endoMR-PREG/plink_mac_20241022/plink"
  # bfile_path <- "/Users/jonasvibert/Desktop/Stata/endoMR-PREG/plink_mac_20241022/24088632/1000G.EUR.QC"
  # endo <- fread("/Users/jonasvibert/Desktop/Stata/endoMR-PREG/plink_mac_20241022/endo_dat_10000snps.txt", header = FALSE)
   
  # 2.2 clean sheet
  # endo <- endo[-1, ]
  # setnames(endo, 
  #old = names(endo),
  #new = c("chrpos","CHR","BP","A1","A2","EAF","beta.exposure","se.exposure","pval.exposure","N","id.exposure","mr_keep","reported","units","SNP"))
# num_cols <- c("CHR","BP","EAF","beta.exposure","se.exposure","pval.exposure","N")
# endo[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
  
  # 2.2.1 Display table 
  # head(endo, 3)
  
  # check again
  # sum(is.na(endo$CHR))  # should be zero
  
  # 2.3. Create the minimal PLINK‐clump input
  #clump_input <- endo %>%
  #  transmute( SNP = SNP, CHR = as.integer(CHR),  BP  = as.integer(BP), P   = pval.exposure)
  
  # 2.4. Write to a tab‐delimited text file
  # write.table(clump_input, file      = "clump_input.txt", sep       = "\t", quote     = FALSE, row.names = FALSE,col.names = TRUE)
  
  # 7) Run PLINK for clumping
  # out_prefix <- file.path(results_dir, "plink_clump_ukbb")
  # system2(plink_path,args = c( "--bfile", bfile_path,"--clump", "clump_input.txt", "--clump-field", "P", "--clump-snp-field", "SNP", "--clump-p1", "5e-8", "--clump-p2", "1","--clump-r2", "0.001", "--clump-kb", "10000",  "--out", out_prefix))
  
  # 8) Load the clumping results and report how many SNPs were kept
  #library(data.table)
  # clumped <- fread(paste0(out_prefix, ".clumped"), header = TRUE)
  # message("Number of SNPs retained after clumping: ", nrow(clumped))
  
  # 9) Extract only the SNP column
  # snps_table <- clumped[, .(SNP)]
  # kable(snps_table,caption = "Table of 41 SNPs retained after clumping")
  
  
  # compare SNP from Clumped_SNPs_1000G and endo10k_clumped.clumped
  # list2 <- c("rs10917151","rs71575922","rs3858429","rs11674184","rs1903068","rs1971256","rs10122243","rs6456259","rs10090060","rs2421985","rs12320196","rs12441483", "rs12030576","rs10828249","rs1451383","rs73633307","rs4540228","rs10757277", "rs7907732","rs2226158","rs10860864","rs7924571","rs66683298","rs12549438",  "rs7214750","rs11756073","rs56090796","rs57281976","rs3803042","rs2510770", "rs55909142","rs17053711","rs1430787","rs507666","rs2967684","rs7334326", "rs2946160","rs10983311","rs1352889","rs2036754","rs6435157")
  
  # Compare with your clumped SNPs:
  # list1 <- clumped$SNP
  
  # Are they exactly the same set?
  #setequal(list1, list2)
  #> TRUE
  
  # To see any mismatches:
  #setdiff(list1, list2)  # SNPs in list1 but not in list2
  #setdiff(list2, list1)  # SNPs in list2 but not in list1
  
  
  ###########################
  # Calculate pseudo F‐stat #
  ###########################
  
  # 1) Retrieve SNPs and their statistics (including EAF and sample size)
  stats <- expdat_ready %>%
    select(SNP, beta.exposure, se.exposure, pval.exposure, eaf.exposure, samplesize.exposure)
  
  # 2) Perform a join to add these statistics to the clumped results
  library(dplyr)
  clumped2 <- clumped %>%
    left_join(stats, by = "SNP")
  
  # 3) Verify that the columns are present
  names(clumped2)
  clumped2 %>%
    select(
      SNP,
      beta.exposure = beta.exposure.y,
      se.exposure   = se.exposure.y,
      eaf.exposure  = eaf.exposure.y,
      samplesize.exposure = samplesize.exposure.y
    ) %>%
    head()

  # Create or replace simple F‐statistic: (β / SE)²
  clumped2 <- clumped2 %>%
    mutate(
      beta.exposure       = coalesce(beta.exposure.y, beta.exposure.x),
      se.exposure         = coalesce(se.exposure.y,   se.exposure.x),
      eaf.exposure        = coalesce(eaf.exposure.y,  eaf.exposure.x),
      samplesize.exposure = coalesce(samplesize.exposure.y, samplesize.exposure.x),
      F_simple            = (beta.exposure / se.exposure)^2
    )
  clumped2 <- clumped2 %>%
    mutate(F_simple = (beta.exposure / se.exposure)^2)
  
  # 4) Inspect F_simple
  clumped2 %>%
    select(SNP, beta.exposure, se.exposure, F_simple) %>%
    head()
  
  # 5) Approximate per‐SNP R²
  clumped2 <- clumped2 %>%
    mutate(
      R2 = (2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2) /
        (2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2 +
           se.exposure^2 * samplesize.exposure)
    )
  
  # 6) Exact F‐statistic with df1 = 1, df2 = N – 2
  clumped2 <- clumped2 %>%
    mutate(F_exact = R2 * (samplesize.exposure - 2) / (1 - R2))
  
  # 7) Quick check: first 10 rows
  head(
    clumped2[, c("SNP", "beta.exposure", "se.exposure", "F_simple", "R2", "F_exact")],
    10
  )
  
  # 8) Mean F‐statistics
  cat("Mean F_simple =", round(mean(clumped2$F_simple, na.rm = TRUE), 2), "\n")
  cat("Mean F_exact  =", round(mean(clumped2$F_exact,  na.rm = TRUE), 2), "\n")
  
  ###############################################################################
  # Variance Explained by Selected SNPs
  ###############################################################################
  
  # We quantify how much of the variability in our exposure
  # (“Endometriosis (Rahmioglu)”) is captured by the SNPs retained after
  # LD‐clumping. We compute:
  #
  # 1. Per‐SNP R²:       R2_i = 2 * EAF_i * (1 − EAF_i) * β_i²
  # 2. Total R²:         sum(R2_i) across all selected SNPs
  # 3. Multi‐SNP F‑stat: F = (R²/k) * ((N − k − 1)/(1 − R²))
  # 4. Mean per‑SNP F‑statistic as a measure of individual instrument strength
  #
  # We start from the data frame `clumped2` as defined above, which contains
  # the columns:
  #   – beta.exposure       (β_i)
  #   – se.exposure         (SE_i)
  #   – eaf.exposure        (EAF_i)
  #   – samplesize.exposure (N, identical across SNPs)
  
  # 1. Compute per‐SNP R² and exact per‐SNP F‑statistics
  clumped2 <- clumped2 %>%
    mutate(
      R2_i = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2,
      F_i  = R2_i * (samplesize.exposure - 2) / (1 - R2_i)
    )
  
  # 2. Sum R² across all SNPs to get total explained variance
  total_R2 <- sum(clumped2$R2_i, na.rm = TRUE)
  
  # 3. Calculate multi‐SNP F‑statistic
  k       <- nrow(clumped2)                     # number of instruments
  N <- as.integer(median(clumped2$samplesize.exposure, na.rm = TRUE))
  F_multi <- (total_R2 / k) * ((N - k - 1) / (1 - total_R2))
  
  # 4. Compute mean per‐SNP F‑statistic
  mean_F  <- mean(clumped2$F_i, na.rm = TRUE)
  
  # 5. Report results
  cat(sprintf(
    "Total variance explained (R²)        = %.6g\n", total_R2
  ))
  cat(sprintf(
    "Multi‑SNP F‑statistic               = %.3f\n", F_multi
  ))
  cat(sprintf(
    "Mean per‑SNP F‑statistic            = %.3f\n", mean_F
  ))
  
  
  
