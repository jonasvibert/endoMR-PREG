#!/usr/bin/env Rscript
###############################################################################
# endoMR-PREG Configuration File
# Centralized configuration for all analysis scripts
###############################################################################

# Project metadata
PROJECT_NAME <- "endoMR-PREG"
PROJECT_VERSION <- "1.0.0"
ANALYSIS_DATE <- Sys.Date()

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

# Full genome-wide Rahmioglu et al. 2023 sumstats (GCST90205183, GRCh37).
# Not redistributed with this repo (too large, external license). Set the
# RAHMIOGLU_GW_FILE env var to override, or edit the default below.
RAHMIOGLU_GW_FILE <- Sys.getenv(
  "RAHMIOGLU_GW_FILE",
  unset = "/Volumes/My Passport/GWAS/endometriosis_Rahmioglu_2023/GCST90205183_buildGRCh37.tsv"
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

# Logging configuration
LOGGING <- list(
  level = "INFO",  # DEBUG, INFO, WARN, ERROR
  console = TRUE,
  file = TRUE
)

cat("Configuration loaded successfully for endoMR-PREG v", PROJECT_VERSION, "\n")
cat("Project root:", PROJECT_ROOT, "\n")
cat("Analysis date:", as.character(ANALYSIS_DATE), "\n")
