#!/usr/bin/env Rscript
###############################################################################
# endoMR-PREG Configuration File
# Centralized configuration for all analysis scripts
###############################################################################

# Project metadata
PROJECT_NAME <- "endoMR-PREG"
PROJECT_VERSION <- "1.0.0"
ANALYSIS_DATE <- Sys.Date()

# Required R version
REQUIRED_R_VERSION <- "4.0.0"

# Package versions (for reproducibility)
REQUIRED_PACKAGES <- list(
  TwoSampleMR = "0.5.6",
  MRInstruments = NULL,  # Latest version
  MRPRESSO = NULL,       # Latest from GitHub
  dplyr = "1.1.0",
  ggplot2 = "3.4.0",
  data.table = "1.14.0",
  here = "1.0.1",
  readr = "2.1.0",
  openxlsx = "4.2.5",
  knitr = "1.42",
  tidyr = "1.3.0",
  metafor = "4.0.0",
  forestplot = "3.1.0",
  progress = "1.2.2",
  forcats = "1.0.0",
  scales = "1.2.1"
)

# Script-specific package requirements
REQUIRED_PACKAGES_PREP <- c("data.table", "dplyr", "readr", "TwoSampleMR", "here", "progress", "ggplot2", "forestplot", "metafor")
REQUIRED_PACKAGES_HARMONISE <- c("TwoSampleMR", "dplyr", "data.table", "here")
REQUIRED_PACKAGES_MR <- c("TwoSampleMR", "dplyr", "here", "MRPRESSO")
REQUIRED_PACKAGES_PLOTS <- c("TwoSampleMR", "ggplot2", "dplyr", "here", "forcats", "scales", "metafor")

# Expected number of SNPs (for validation)
EXPECTED_N_SNPS <- 41

# Directory structure
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}
library(here)

PROJECT_ROOT <- here::here()
DATA_DIR <- here::here("data")
RESULTS_DIR <- here::here("results")
PLOTS_DIR <- here::here("plots")
SCRIPTS_DIR <- here::here("scripts")
CONFIG_DIR <- here::here("config")
LOGS_DIR <- here::here("logs")

# Create directories if they don't exist
REQUIRED_DIRS <- c(DATA_DIR, RESULTS_DIR, PLOTS_DIR, SCRIPTS_DIR, CONFIG_DIR, LOGS_DIR)
for (dir in REQUIRED_DIRS) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
}

# Data file paths
EXPOSURE_FILE <- file.path(DATA_DIR, "EXPOSURE_RAHMIGLU", "NIHMS1873427-Supplementary_Materials.xlsx")
SNP_LIST_FILE <- file.path(RESULTS_DIR, "endometriosis_clumped_snps.tsv")
HARMONISED_FILE <- file.path(RESULTS_DIR, "harmonised_rahmioglu_bpo.csv")

# PLINK configuration
PLINK_DIR <- here::here("plink_mac_20241022")
PLINK_PATH <- file.path(PLINK_DIR, "plink")
BFILE_PATH <- file.path(PLINK_DIR, "24088632")
REFERENCE_PANEL <- file.path(BFILE_PATH, "1000G.EUR.QC")

# MR analysis parameters
MR_PARAMS <- list(
  # Clumping parameters
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 5e-8,
  clump_p2 = 1,
  
  # MR-PRESSO parameters
  presso_nb_distribution = 5000,  # Increased from 1000
  presso_signif_threshold = 0.05,
  
  # Significance thresholds
  alpha_nominal = 0.05,
  alpha_bonferroni = NULL,  # Will be calculated as 0.05/n_outcomes
  alpha_fdr = 0.05,
  
  # Minimum SNPs for analysis
  min_snps_mr = 3,
  min_snps_sensitivity = 4
)

# Plot parameters
PLOT_PARAMS <- list(
  # Standard dimensions
  scatter_width = 7,
  scatter_height = 7,
  forest_width = 10,
  forest_height = 8,
  funnel_width = 7,
  funnel_height = 7,
  
  # High-resolution settings
  dpi = 300,
  
  # Color palettes
  significance_colors = c(
    "NS" = "grey70",
    "p < 0.05" = "#2b8cbe",
    "FDR < 0.05" = "#e6550d", 
    "Bonferroni < 0.05" = "#d62728"
  ),
  
  mr_method_colors = c(
    "Inverse variance weighted" = "#E31A1C",
    "MR Egger" = "#1F78B4",
    "Weighted median" = "#33A02C"
  )
)

# External URLs (with fallbacks)
EXTERNAL_URLS <- list(
  finngen_r12_manifest = "https://zenodo.org/records/15815805/files/finngen_R12_manifest.csv",
  finngen_r12_base = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/"
)

# Data source file lists
PPH_FILES <- c(
  "Postpartum_hemorrhage",
  "Postpartum_hemorrhage_due_to_atony",
  "Postpartum_hemorrhage_due_to_retained_placenta",
  "Early_bleeding_with_any_outcome",
  "Antepartum_bleeding",
  "Early_bleeding_ending_in_live_birth"
)

FINNGEN_FILES <- c(
  "finngen_R12_O15_PLAC_PRAEVIA",
  "finngen_R12_O15_PLAC_DISORD",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",
  "finngen_R12_O15_PREG_ECTOP",
  "finngen_R12_N14_FEMALEINFERT"
)

# Path collections for validation
PATHS <- list(
  data = DATA_DIR,
  results = RESULTS_DIR,
  plots = PLOTS_DIR,
  scripts = SCRIPTS_DIR,
  config = CONFIG_DIR,
  logs = LOGS_DIR
)

REQUIRED_PATHS_PREP <- list(
  project_root = PROJECT_ROOT,
  data_dir = DATA_DIR,
  results_dir = RESULTS_DIR
)

# Logging configuration
LOGGING <- list(
  console = TRUE,
  file = TRUE,
  level = "INFO"
)

# Comprehensive outcome mappings with clean, publication-ready names
OUTCOME_LABELS <- list(
  # ========================================================================
  # PREGNANCY COMPLICATIONS & DELIVERY
  # ========================================================================
  
  # Cesarean delivery
  el_cs = "Elective cesarean delivery",
  em_cs = "Emergency cesarean delivery", 
  cs = "Cesarean delivery (any)",
  
  # Preterm birth spectrum
  pretb_all = "Preterm birth",
  pretb_subsamp = "Preterm birth (replication)",
  vpretb_all = "Very preterm birth (<32 weeks)",
  posttb_all = "Post-term birth (>42 weeks)",
  
  # Membrane rupture
  rup_memb = "Prelabor rupture of membranes",
  
  # Labor complications
  induction = "Labor induction",
  
  # ========================================================================
  # PLACENTAL DISORDERS
  # ========================================================================
  
  # FinnGen placental outcomes
  "finngen_R12_O15_PLAC_PRAEVIA" = "Placenta previa",
  "finngen_R12_O15_PLAC_PREMAT_SEPAR" = "Placental abruption", 
  "finngen_R12_O15_PLAC_DISORD" = "Placental disorders (other)",
  "finngen_R12_O15_PREG_ECTOP" = "Ectopic pregnancy",
  
  # ========================================================================
  # HEMORRHAGIC COMPLICATIONS
  # ========================================================================
  
  # Postpartum hemorrhage spectrum
  "Postpartum_hemorrhage" = "Postpartum hemorrhage",
  "Postpartum_hemorrhage_due_to_atony" = "PPH due to uterine atony",
  "Postpartum_hemorrhage_due_to_retained_placenta" = "PPH due to retained placenta",
  
  # Antepartum bleeding
  "Antepartum_bleeding" = "Antepartum bleeding",
  "Early_bleeding_with_any_outcome" = "Early pregnancy bleeding",
  "Early_bleeding_ending_in_live_birth" = "Early bleeding with live birth",
  
  # ========================================================================
  # HYPERTENSIVE DISORDERS
  # ========================================================================
  
  pe_subsamp = "Preeclampsia",
  gh_subsamp = "Gestational hypertension", 
  hdp_subsamp = "Hypertensive disorders of pregnancy",
  
  # ========================================================================
  # METABOLIC COMPLICATIONS
  # ========================================================================
  
  gdm_subsamp = "Gestational diabetes",
  
  # ========================================================================
  # NAUSEA & VOMITING
  # ========================================================================
  
  nvp_sev_all = "Severe nausea and vomiting",
  nvp_sev_subsamp = "Severe nausea and vomiting (replication)",
  hyp = "Hyperemesis gravidarum",
  
  # ========================================================================
  # PREGNANCY LOSS
  # ========================================================================
  
  misc_subsamp = "Miscarriage", 
  s_misc_subsamp = "Sporadic miscarriage",
  r_misc_subsamp = "Recurrent miscarriage",
  sb_subsamp = "Stillbirth",
  
  # ========================================================================
  # FERTILITY & CONCEPTION
  # ========================================================================
  
  "finngen_R12_N14_FEMALEINFERT" = "Female infertility",
  
  # ========================================================================
  # FETAL GROWTH & DEVELOPMENT
  # ========================================================================
  
  # Birth weight spectrum
  lbw_all = "Low birth weight (<2500g)",
  hbw_all = "High birth weight (>4000g)", 
  zbw_all = "Birth weight (standardized)",
  
  # Gestational age assessment
  ga_all = "Gestational age",
  ga_subsamp = "Gestational age (replication)",
  
  # Size for gestational age
  sga = "Small for gestational age",
  lga = "Large for gestational age",
  
  # ========================================================================
  # NEONATAL OUTCOMES
  # ========================================================================
  
  # Apgar scores
  apgar1 = "Apgar score (1 minute)", 
  apgar5 = "Apgar score (5 minutes)",
  lowapgar1 = "Low Apgar score (1 minute)",
  lowapgar5 = "Low Apgar score (5 minutes)",
  
  # Neonatal care
  nicu = "NICU admission",
  
  # ========================================================================
  # MATERNAL HEALTH CONDITIONS
  # ========================================================================
  
  anaemia_preg_all = "Pregnancy anemia",
  depr_subsamp = "Postpartum depression",
  
  # ========================================================================
  # BREASTFEEDING OUTCOMES  
  # ========================================================================
  
  bf_ini = "Breastfeeding initiation",
  bf_est = "Exclusive breastfeeding",
  bf_dur_4c = "Breastfeeding ≥4 months", 
  bf_sus = "Breastfeeding duration"
)

# Add filtered versions
OUTCOME_LABELS_FILTERED <- OUTCOME_LABELS
for (outcome in names(OUTCOME_LABELS)) {
  filtered_name <- paste0(outcome, "_filtered")
  OUTCOME_LABELS_FILTERED[[filtered_name]] <- OUTCOME_LABELS[[outcome]]
}

# Outcome categories for grouping in tables and figures
OUTCOME_CATEGORIES <- list(
  "Delivery complications" = c("el_cs", "em_cs", "cs", "induction"),
  
  "Preterm birth & timing" = c("pretb_all", "pretb_subsamp", "vpretb_all", "posttb_all", 
                               "ga_all", "ga_subsamp"),
  
  "Placental disorders" = c("finngen_R12_O15_PLAC_PRAEVIA", "finngen_R12_O15_PLAC_PREMAT_SEPAR", 
                           "finngen_R12_O15_PLAC_DISORD", "finngen_R12_O15_PREG_ECTOP"),
  
  "Hemorrhagic complications" = c("Postpartum_hemorrhage", "Postpartum_hemorrhage_due_to_atony",
                                 "Postpartum_hemorrhage_due_to_retained_placenta", "Antepartum_bleeding",
                                 "Early_bleeding_with_any_outcome", "Early_bleeding_ending_in_live_birth"),
  
  "Hypertensive disorders" = c("pe_subsamp", "gh_subsamp", "hdp_subsamp"),
  
  "Metabolic conditions" = c("gdm_subsamp"),
  
  "Nausea & vomiting" = c("nvp_sev_all", "nvp_sev_subsamp", "hyp"),
  
  "Pregnancy loss" = c("misc_subsamp", "s_misc_subsamp", "r_misc_subsamp", "sb_subsamp"),
  
  "Fertility" = c("finngen_R12_N14_FEMALEINFERT"),
  
  "Fetal growth" = c("lbw_all", "hbw_all", "zbw_all", "sga", "lga"),
  
  "Neonatal outcomes" = c("apgar1", "apgar5", "lowapgar1", "lowapgar5", "nicu"),
  
  "Maternal health" = c("anaemia_preg_all", "depr_subsamp", "rup_memb"),
  
  "Breastfeeding" = c("bf_ini", "bf_est", "bf_dur_4c", "bf_sus")
)

# Create reverse mapping (outcome -> category)
OUTCOME_TO_CATEGORY <- list()
for (category in names(OUTCOME_CATEGORIES)) {
  for (outcome in OUTCOME_CATEGORIES[[category]]) {
    OUTCOME_TO_CATEGORY[[outcome]] <- category
    # Also add filtered versions
    OUTCOME_TO_CATEGORY[[paste0(outcome, "_filtered")]] <- category
  }
}

# Data source mapping
DATA_SOURCE_MAPPING <- list(
  # MR-PREG outcomes (most outcomes)
  mrpreg = c("pretb_all", "pretb_subsamp", "vpretb_all", "posttb_all", "el_cs", "em_cs", "cs",
             "ga_all", "ga_subsamp", "lbw_all", "hbw_all", "zbw_all", "sga", "lga",
             "apgar1", "apgar5", "lowapgar1", "lowapgar5", "nicu", "pe_subsamp", "gh_subsamp",
             "hdp_subsamp", "gdm_subsamp", "nvp_sev_all", "nvp_sev_subsamp", "hyp",
             "misc_subsamp", "s_misc_subsamp", "r_misc_subsamp", "sb_subsamp",
             "anaemia_preg_all", "depr_subsamp", "rup_memb", "induction",
             "bf_ini", "bf_est", "bf_dur_4c", "bf_sus"),
  
  # FinnGen outcomes
  finngen = c("finngen_R12_N14_FEMALEINFERT", "finngen_R12_O15_PLAC_PRAEVIA", 
             "finngen_R12_O15_PLAC_PREMAT_SEPAR", "finngen_R12_O15_PLAC_DISORD", 
             "finngen_R12_O15_PREG_ECTOP"),
  
  # PPH/Westergaard outcomes
  pph = c("Postpartum_hemorrhage", "Postpartum_hemorrhage_due_to_atony",
          "Postpartum_hemorrhage_due_to_retained_placenta", "Antepartum_bleeding",
          "Early_bleeding_with_any_outcome", "Early_bleeding_ending_in_live_birth")
)

# Create reverse mapping (outcome -> data source)
OUTCOME_TO_SOURCE <- list()
for (source in names(DATA_SOURCE_MAPPING)) {
  for (outcome in DATA_SOURCE_MAPPING[[source]]) {
    OUTCOME_TO_SOURCE[[outcome]] <- source
    # Also add filtered versions
    OUTCOME_TO_SOURCE[[paste0(outcome, "_filtered")]] <- source
  }
}

# Priority outcomes for main figures (most clinically relevant)
PRIORITY_OUTCOMES <- c(
  "pretb_all",                                    # Preterm birth
  "finngen_R12_O15_PLAC_PRAEVIA",                # Placenta previa
  "finngen_R12_N14_FEMALEINFERT",                # Female infertility
  "el_cs",                                        # Elective cesarean
  "finngen_R12_O15_PLAC_PREMAT_SEPAR",          # Placental abruption
  "rup_memb",                                     # PROM
  "pe_subsamp",                                   # Preeclampsia
  "Postpartum_hemorrhage",                        # PPH
  "misc_subsamp",                                 # Miscarriage
  "lbw_all",                                      # Low birth weight
  "sga"                                           # Small for GA
)

# File patterns for different data sources
FILE_PATTERNS <- list(
  pph_files = c(
    "Postpartum_hemorrhage",
    "Postpartum_hemorrhage_due_to_atony", 
    "Postpartum_hemorrhage_due_to_retained_placenta",
    "Early_bleeding_with_any_outcome",
    "Antepartum_bleeding",
    "Early_bleeding_ending_in_live_birth"
  ),
  
  finngen_files = c(
    "finngen_R12_O15_PLAC_PRAEVIA",
    "finngen_R12_O15_PLAC_DISORD",
    "finngen_R12_O15_PLAC_PREMAT_SEPAR",
    "finngen_R12_O15_PREG_ECTOP",
    "finngen_R12_N14_FEMALEINFERT"
  )
)

# Column mappings for different data sources
COLUMN_MAPPINGS <- list(
  pph_columns = c(
    "Chr", "PosB38", "Marker", "rsName", "OA", "EA",
    "EAfrq", "Effect", "P", "Phet", "I2", "nCoh"
  ),
  
  finngen_columns = c(
    "#chrom", "pos", "ref", "alt", "rsids",
    "pval", "beta", "sebeta",
    "af_alt", "af_alt_cases", "af_alt_controls"
  )
)

# Logging configuration
LOGGING <- list(
  level = "INFO",  # DEBUG, INFO, WARN, ERROR
  console = TRUE,
  file = TRUE,
  max_file_size = "10MB"
)

cat("Configuration loaded successfully for endoMR-PREG v", PROJECT_VERSION, "\n")
cat("Project root:", PROJECT_ROOT, "\n")
cat("Analysis date:", as.character(ANALYSIS_DATE), "\n")