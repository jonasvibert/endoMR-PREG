#!/usr/bin/env Rscript
###############################################################################
# endoMR-PREG Utility Functions
# Logging helpers used by scripts/06_tables_endoMR-PREG.R
###############################################################################

# Load configuration
if (!exists("PROJECT_ROOT")) {
  source(here::here("config", "config.R"))
}

###############################################################################
# LOGGING FUNCTIONS
###############################################################################

#' Initialize logging for the current script
#' @param script_name Name of the current script
#' @param log_level Logging level (DEBUG, INFO, WARN, ERROR)
init_logging <- function(script_name, log_level = "INFO") {

  # Create logs directory if it doesn't exist
  if (!dir.exists(LOGS_DIR)) {
    dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
  }

  # Set up log file
  log_file <- file.path(LOGS_DIR, paste0(script_name, "_", Sys.Date(), ".log"))

  # Initialize global logging variables
  assign(".log_file", log_file, envir = .GlobalEnv)
  assign(".log_level", log_level, envir = .GlobalEnv)
  assign(".script_name", script_name, envir = .GlobalEnv)

  log_info("=== Starting script:", script_name, "===")
  log_info("Log file:", log_file)
  log_info("R version:", R.version.string)
  log_info("Working directory:", getwd())
}

#' Log message with timestamp and level
#' @param level Log level
#' @param ... Messages to log
log_message <- function(level, ...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  script <- if(exists(".script_name")) get(".script_name") else "unknown"
  message_text <- paste(..., collapse = " ")

  log_entry <- sprintf("[%s] [%s] [%s] %s", timestamp, level, script, message_text)

  # Print to console if enabled
  if (LOGGING$console) {
    cat(log_entry, "\n")
  }

  # Write to file if enabled and log file exists
  if (LOGGING$file && exists(".log_file")) {
    cat(log_entry, "\n", file = get(".log_file"), append = TRUE)
  }
}

#' Log info message
log_info <- function(...) log_message("INFO", ...)

#' Log warning message
log_warn <- function(...) log_message("WARN", ...)

#' Log error message
log_error <- function(...) log_message("ERROR", ...)

###############################################################################
# PACKAGE MANAGEMENT
###############################################################################

#' Load required packages safely
#' @param packages Vector of package names to load
load_packages <- function(packages) {

  log_info("Loading packages:", paste(packages, collapse = ", "))

  for (pkg in packages) {
    success <- tryCatch({
      library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      TRUE
    }, error = function(e) {
      log_error("Failed to load package", pkg, ":", e$message)
      FALSE
    })

    if (!success) {
      stop("Critical package loading failure: ", pkg)
    }
  }

  log_info("All packages loaded successfully")
}

###############################################################################
# SCRIPT SETUP
###############################################################################

#' Set up logging for a script
#' @param script_name Name of the script for logging
setup_logging <- function(script_name) {
  init_logging(script_name, log_level = "INFO")
}
