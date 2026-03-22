#!/usr/bin/env Rscript
###############################################################################
# endoMR-PREG Utility Functions
# Common functions used across multiple analysis scripts
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

#' Log debug message
log_debug <- function(...) {
  if (exists(".log_level") && get(".log_level") == "DEBUG") {
    log_message("DEBUG", ...)
  }
}

###############################################################################
# PACKAGE MANAGEMENT FUNCTIONS
###############################################################################

#' Check and install required packages with version control
#' @param packages Named list of packages with versions (NULL for latest)
#' @param repos CRAN repository URL
check_and_install_packages <- function(packages = REQUIRED_PACKAGES, 
                                     repos = "https://cloud.r-project.org") {
  
  log_info("Checking required packages...")
  
  for (pkg_name in names(packages)) {
    required_version <- packages[[pkg_name]]
    
    # Check if package is installed
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      log_info("Installing package:", pkg_name)
      
      # Special handling for GitHub packages
      if (pkg_name == "MRPRESSO") {
        if (!requireNamespace("devtools", quietly = TRUE)) {
          install.packages("devtools", repos = repos)
        }
        devtools::install_github("rondolab/MR-PRESSO")
      } else {
        install.packages(pkg_name, repos = repos)
      }
    } else {
      # Check version if specified
      if (!is.null(required_version)) {
        installed_version <- packageVersion(pkg_name)
        if (installed_version < required_version) {
          log_warn("Package", pkg_name, "version", installed_version, 
                  "is older than required", required_version)
        }
      }
      log_debug("Package", pkg_name, "is available")
    }
  }
  
  log_info("Package check completed")
}

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
# FILE VALIDATION FUNCTIONS
###############################################################################

#' Check if required files exist
#' @param files Vector of file paths to check
#' @param stop_on_missing Whether to stop execution if files are missing
validate_files <- function(files, stop_on_missing = TRUE) {
  
  log_info("Validating", length(files), "required files...")
  
  missing_files <- c()
  
  for (file_path in files) {
    if (!file.exists(file_path)) {
      missing_files <- c(missing_files, file_path)
      log_error("Missing file:", file_path)
    } else {
      log_debug("Found file:", file_path)
    }
  }
  
  if (length(missing_files) > 0) {
    if (stop_on_missing) {
      stop("Missing required files: ", paste(missing_files, collapse = ", "))
    } else {
      log_warn("Missing files found but continuing:", paste(missing_files, collapse = ", "))
      return(FALSE)
    }
  }
  
  log_info("All required files validated")
  return(TRUE)
}

#' Check if required variables exist in environment
#' @param variables Vector of variable names to check
#' @param envir Environment to check (default: global)
#' @param stop_on_missing Whether to stop execution if variables are missing
validate_variables <- function(variables, envir = .GlobalEnv, stop_on_missing = TRUE) {
  
  log_info("Validating", length(variables), "required variables...")
  
  missing_vars <- c()
  
  for (var_name in variables) {
    if (!exists(var_name, envir = envir)) {
      missing_vars <- c(missing_vars, var_name)
      log_error("Missing variable:", var_name)
    } else {
      log_debug("Found variable:", var_name)
    }
  }
  
  if (length(missing_vars) > 0) {
    if (stop_on_missing) {
      stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
    } else {
      log_warn("Missing variables found but continuing:", paste(missing_vars, collapse = ", "))
      return(FALSE)
    }
  }
  
  log_info("All required variables validated")
  return(TRUE)
}

###############################################################################
# DATA PROCESSING FUNCTIONS
###############################################################################

#' Safely read CSV/TSV files with error handling
#' @param file_path Path to the file
#' @param sep Separator (auto-detected if NULL)
#' @param ... Additional arguments passed to read function
safe_read_file <- function(file_path, sep = NULL, ...) {
  
  log_info("Reading file:", basename(file_path))
  
  # Validate file exists
  if (!file.exists(file_path)) {
    log_error("File not found:", file_path)
    return(NULL)
  }
  
  # Auto-detect separator if not provided
  if (is.null(sep)) {
    file_ext <- tools::file_ext(file_path)
    sep <- switch(tolower(file_ext),
                  "csv" = ",",
                  "tsv" = "\t",
                  "txt" = "\t",
                  "\t")  # default to tab
  }
  
  # Try to read the file
  data <- tryCatch({
    if (sep == ",") {
      readr::read_csv(file_path, show_col_types = FALSE, ...)
    } else {
      readr::read_delim(file_path, delim = sep, show_col_types = FALSE, ...)
    }
  }, error = function(e) {
    log_error("Failed to read file", file_path, ":", e$message)
    return(NULL)
  })
  
  if (!is.null(data)) {
    log_info("Successfully read", nrow(data), "rows and", ncol(data), "columns")
  }
  
  return(data)
}

#' Apply outcome labels with fallback
#' @param data Data frame with outcome column
#' @param outcome_col Name of the outcome column
#' @param labels Named vector of outcome labels
apply_outcome_labels <- function(data, outcome_col = "outcome", labels = OUTCOME_LABELS_FILTERED) {
  
  if (!outcome_col %in% names(data)) {
    log_warn("Outcome column", outcome_col, "not found in data")
    return(data)
  }
  
  # Create outcome_full column
  data$outcome_full <- ifelse(
    data[[outcome_col]] %in% names(labels),
    labels[data[[outcome_col]]],
    data[[outcome_col]]
  )
  
  # Log mapping statistics
  mapped <- sum(data[[outcome_col]] %in% names(labels))
  total <- nrow(data)
  log_info("Mapped", mapped, "of", total, "outcomes to readable labels")
  
  return(data)
}

#' Calculate significance levels with multiple correction
#' @param pvalues Vector of p-values
#' @param alpha_nominal Nominal significance level
#' @param methods Vector of correction methods to apply
calculate_significance <- function(pvalues, alpha_nominal = 0.05, 
                                 methods = c("bonferroni", "fdr")) {
  
  result <- data.frame(
    pval = pvalues,
    sig_nominal = pvalues < alpha_nominal
  )
  
  # Apply corrections
  for (method in methods) {
    adj_pvals <- p.adjust(pvalues, method = method)
    result[[paste0("padj_", method)]] <- adj_pvals
    result[[paste0("sig_", method)]] <- adj_pvals < alpha_nominal
  }
  
  # Combined significance column
  result$significance <- ifelse(
    result$sig_bonferroni, "Bonferroni < 0.05",
    ifelse(result$sig_fdr, "FDR < 0.05",
           ifelse(result$sig_nominal, "p < 0.05", "NS"))
  )
  
  return(result)
}

#' Add outcome categories and data sources
#' @param data Data frame with outcome column
#' @param outcome_col Name of the outcome column
add_outcome_metadata <- function(data, outcome_col = "outcome") {
  
  if (!outcome_col %in% names(data)) {
    log_warn("Outcome column", outcome_col, "not found in data")
    return(data)
  }
  
  # Add category information
  data$outcome_category <- sapply(data[[outcome_col]], function(x) {
    category <- OUTCOME_TO_CATEGORY[[x]]
    if (is.null(category)) "Other" else category
  })
  
  # Add data source information
  data$data_source <- sapply(data[[outcome_col]], function(x) {
    source <- OUTCOME_TO_SOURCE[[x]]
    if (is.null(source)) {
      # Try to infer from outcome name
      if (grepl("^finngen", x)) {
        "FinnGen"
      } else if (grepl("hemorrhage|bleeding", x, ignore.case = TRUE)) {
        "Westergaard (PPH)" 
      } else {
        "MR-PREG"
      }
    } else {
      switch(source,
             "mrpreg" = "MR-PREG",
             "finngen" = "FinnGen", 
             "pph" = "Westergaard (PPH)",
             source)
    }
  })
  
  # Add priority flag
  data$is_priority <- data[[outcome_col]] %in% PRIORITY_OUTCOMES
  
  log_info("Added outcome metadata for", nrow(data), "observations")
  
  return(data)
}

#' Format outcome data for publication tables
#' @param data MR results data frame
#' @param include_metadata Whether to include category/source columns
format_publication_table <- function(data, include_metadata = TRUE) {
  
  log_info("Formatting data for publication table")
  
  # Apply outcome labels
  data <- apply_outcome_labels(data)
  
  # Add metadata if requested
  if (include_metadata) {
    data <- add_outcome_metadata(data)
  }
  
  # Format numeric columns
  if ("b" %in% names(data)) {
    data$OR <- exp(data$b)
    data$OR_LCL <- exp(data$b - 1.96 * data$se)
    data$OR_UCL <- exp(data$b + 1.96 * data$se)
    data$OR_CI <- sprintf("%.2f (%.2f–%.2f)", data$OR, data$OR_LCL, data$OR_UCL)
  }
  
  # Format p-values
  if ("pval" %in% names(data)) {
    data$P_formatted <- ifelse(data$pval < 0.001, 
                              formatC(data$pval, format = "e", digits = 2),
                              sprintf("%.3f", data$pval))
  }
  
  # Calculate significance levels
  if ("pval" %in% names(data)) {
    sig_data <- calculate_significance(data$pval)
    data <- cbind(data, sig_data[, !names(sig_data) %in% "pval"])
  }
  
  log_info("Publication table formatted successfully")
  
  return(data)
}

#' Create grouped forest plot data
#' @param mr_results MR results with formatted outcomes
#' @param group_by Grouping variable ("outcome_category", "data_source", etc.)
create_grouped_plot_data <- function(mr_results, group_by = "outcome_category") {
  
  log_info("Creating grouped plot data by", group_by)
  
  # Apply labels and metadata
  plot_data <- mr_results %>%
    apply_outcome_labels() %>%
    add_outcome_metadata()
  
  # Filter for IVW method if multiple methods present
  if ("method" %in% names(plot_data)) {
    plot_data <- plot_data %>%
      dplyr::filter(grepl("Inverse variance weighted", method, ignore.case = TRUE))
  }
  
  # Calculate effect sizes
  if ("b" %in% names(plot_data)) {
    plot_data <- plot_data %>%
      dplyr::mutate(
        OR = exp(b),
        OR_LCL = exp(b - 1.96 * se),
        OR_UCL = exp(b + 1.96 * se)
      )
  }
  
  # Order by group and effect size
  if (group_by %in% names(plot_data)) {
    plot_data <- plot_data %>%
      dplyr::arrange(.data[[group_by]], OR) %>%
      dplyr::mutate(
        outcome_ordered = factor(outcome_full, levels = unique(outcome_full))
      )
  }
  
  log_info("Grouped plot data created with", nrow(plot_data), "observations")
  
  return(plot_data)
}

#' Filter outcomes by priority or category
#' @param data Data frame with outcome information
#' @param filter_type Type of filter ("priority", "category", "source")
#' @param filter_values Specific values to filter (for category/source)
filter_outcomes <- function(data, filter_type = "priority", filter_values = NULL) {
  
  # Add metadata if not present
  if (!"is_priority" %in% names(data)) {
    data <- add_outcome_metadata(data)
  }
  
  filtered_data <- switch(filter_type,
    "priority" = data %>% dplyr::filter(is_priority),
    "category" = {
      if (is.null(filter_values)) {
        log_warn("No categories specified for filtering")
        data
      } else {
        data %>% dplyr::filter(outcome_category %in% filter_values)
      }
    },
    "source" = {
      if (is.null(filter_values)) {
        log_warn("No sources specified for filtering")
        data
      } else {
        data %>% dplyr::filter(data_source %in% filter_values)
      }
    },
    data  # default: return unchanged
  )
  
  log_info("Filtered outcomes:", filter_type, "- kept", nrow(filtered_data), "of", nrow(data), "outcomes")
  
  return(filtered_data)
}

###############################################################################
# MR-SPECIFIC FUNCTIONS
###############################################################################

#' Run MR analysis with comprehensive error handling
#' @param harmonised_data Harmonised exposure-outcome data
#' @param methods Vector of MR methods to run
run_mr_analysis <- function(harmonised_data, methods = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")) {
  
  log_info("Running MR analysis with", length(methods), "methods")
  
  # Validate input data
  required_cols <- c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
  if (!all(required_cols %in% names(harmonised_data))) {
    log_error("Missing required columns for MR analysis")
    return(NULL)
  }
  
  # Filter out invalid data
  valid_data <- harmonised_data %>%
    dplyr::filter(
      !is.na(beta.exposure), !is.na(se.exposure),
      !is.na(beta.outcome), !is.na(se.outcome),
      se.exposure > 0, se.outcome > 0,
      is.finite(beta.exposure), is.finite(se.exposure),
      is.finite(beta.outcome), is.finite(se.outcome)
    )
  
  log_info("Filtered to", nrow(valid_data), "valid observations from", nrow(harmonised_data))
  
  # Run MR analysis
  mr_results <- tryCatch({
    TwoSampleMR::mr(valid_data, method_list = methods)
  }, error = function(e) {
    log_error("MR analysis failed:", e$message)
    return(NULL)
  })
  
  if (!is.null(mr_results)) {
    log_info("MR analysis completed successfully")
    # Add significance calculations
    mr_results <- mr_results %>%
      dplyr::group_by(outcome) %>%
      dplyr::do({
        sig_data <- calculate_significance(.$pval)
        dplyr::bind_cols(., sig_data[, !names(sig_data) %in% "pval"])
      }) %>%
      dplyr::ungroup()
  }
  
  return(mr_results)
}

#' Export results to CSV with consistent formatting
#' @param data Data frame to export
#' @param filename Output filename (without path)
#' @param output_dir Output directory (default: RESULTS_DIR)
export_results <- function(data, filename, output_dir = RESULTS_DIR) {
  
  if (is.null(data) || nrow(data) == 0) {
    log_warn("No data to export for", filename)
    return(FALSE)
  }
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  output_path <- file.path(output_dir, filename)
  
  # Apply outcome labels if outcome column exists
  if ("outcome" %in% names(data)) {
    data <- apply_outcome_labels(data)
  }
  
  # Format numeric columns
  numeric_cols <- sapply(data, is.numeric)
  data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
    ifelse(abs(x) < 0.001 & x != 0, formatC(x, format = "e", digits = 3), round(x, 4))
  })
  
  # Export
  success <- tryCatch({
    readr::write_csv(data, output_path)
    log_info("Exported results to:", filename)
    TRUE
  }, error = function(e) {
    log_error("Failed to export", filename, ":", e$message)
    FALSE
  })
  
  return(success)
}

###############################################################################
# PLOTTING FUNCTIONS
###############################################################################

#' Save plot with error handling and consistent formatting
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution in DPI
#' @param output_dir Output directory
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300, output_dir = PLOTS_DIR) {
  
  if (is.null(plot)) {
    log_warn("No plot to save for", filename)
    return(FALSE)
  }
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  output_path <- file.path(output_dir, filename)
  
  success <- tryCatch({
    ggplot2::ggsave(
      filename = output_path,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    log_info("Saved plot:", filename)
    TRUE
  }, error = function(e) {
    log_error("Failed to save plot", filename, ":", e$message)
    FALSE
  })
  
  return(success)
}

###############################################################################
# SYSTEM UTILITIES
###############################################################################

#' Get system information for reproducibility
get_system_info <- function() {
  info <- list(
    r_version = R.version.string,
    platform = R.version$platform,
    os = Sys.info()["sysname"],
    user = Sys.info()["user"],
    date = Sys.Date(),
    time = Sys.time(),
    working_directory = getwd(),
    project_root = PROJECT_ROOT
  )
  
  return(info)
}

#' Create analysis summary report
create_analysis_summary <- function(script_name, start_time = NULL, additional_info = NULL) {
  
  if (is.null(start_time)) start_time <- Sys.time()
  
  end_time <- Sys.time()
  duration <- end_time - start_time
  
  summary_info <- list(
    script = script_name,
    project = PROJECT_NAME,
    version = PROJECT_VERSION,
    start_time = start_time,
    end_time = end_time,
    duration = duration,
    system_info = get_system_info()
  )
  
  if (!is.null(additional_info)) {
    summary_info <- c(summary_info, additional_info)
  }
  
  log_info("=== Analysis Summary ===")
  log_info("Script:", script_name)
  log_info("Duration:", round(as.numeric(duration, units = "mins"), 2), "minutes")
  log_info("Completed at:", format(end_time, "%Y-%m-%d %H:%M:%S"))
  
  return(summary_info)
}

###############################################################################
# DATA PREPARATION SPECIFIC FUNCTIONS
###############################################################################

#' Set up logging specifically for data preparation
#' @param script_name Name of the script for logging
setup_logging <- function(script_name) {
  init_logging(script_name, log_level = "INFO")
}

#' Validate dependencies for data preparation
#' @param required_packages Vector of required package names
validate_dependencies <- function(required_packages) {
  log_info("Validating dependencies...")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      log_error("Required package not found: {pkg}")
      stop("Missing required package: ", pkg)
    }
  }
  
  log_info("All dependencies validated")
}

#' Validate file paths for data preparation
#' @param required_paths List of required paths
validate_file_paths <- function(required_paths) {
  log_info("Validating file paths...")
  
  missing_paths <- c()
  for (path_name in names(required_paths)) {
    path <- required_paths[[path_name]]
    if (!file.exists(path) && !dir.exists(path)) {
      missing_paths <- c(missing_paths, paste0(path_name, ": ", path))
      log_error("Missing path: {path_name} -> {path}")
    }
  }
  
  if (length(missing_paths) > 0) {
    log_warning("Some paths missing - this may be expected if running for first time")
  } else {
    log_info("All file paths validated")
  }
}

#' Validate data directories exist
validate_data_directories <- function() {
  log_info("Validating data directories...")
  
  required_dirs <- c(
    PATHS$data,
    PATHS$results,
    file.path(PATHS$data, "OUTCOME_MR-PREG"),
    file.path(PATHS$data, "OUTCOME_PPH"),
    file.path(PATHS$data, "OUTCOME_FINNGEN")
  )
  
  for (dir_path in required_dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      log_info("Created directory: {dir_path}")
    }
  }
  
  log_info("Data directories validated")
}

#' Load SNP instruments from available sources
load_snp_list <- function() {
  log_info("Loading SNP instruments...")
  
  # Try multiple sources for SNP list
  snps_file <- file.path(PATHS$results, "endo_instruments_41.txt")
  clumped_file <- file.path(PATHS$results, "endometriosis_clumped_snps.tsv")
  
  if (file.exists(snps_file)) {
    log_info("Loading SNPs from: {basename(snps_file)}")
    snps <- scan(snps_file, what = character(), sep = "\n", quiet = TRUE, strip.white = TRUE)
    snps <- unique(snps[nzchar(snps)])
    snps <- snps[!tolower(snps) %in% c("snp", "rsid")]
  } else if (file.exists(clumped_file)) {
    log_info("Loading SNPs from: {basename(clumped_file)}")
    tmp <- data.table::fread(clumped_file, select = "SNP", showProgress = FALSE)
    snps <- unique(tmp$SNP)
  } else {
    stop("Could not find SNP instruments. Run script 1 first or place 'endo_instruments_41.txt' in results directory")
  }
  
  # Validate SNP format (basic check for rsID format)
  valid_snps <- snps[grepl("^rs[0-9]+$", snps)]
  if (length(valid_snps) < length(snps)) {
    log_warning("{length(snps) - length(valid_snps)} SNPs do not match rsID format")
  }
  
  return(valid_snps)
}

#' Process MR-PREG outcomes
#' @param snps Vector of SNP rsIDs to filter for
process_mr_preg_outcomes <- function(snps) {
  log_info("Processing MR-PREG adverse pregnancy outcomes...")
  
  mr_preg_file <- file.path(PATHS$data, "OUTCOME_MR-PREG", "ma_out_dat.txt")
  
  if (!file.exists(mr_preg_file)) {
    log_warning("MR-PREG file not found: {mr_preg_file}")
    return(data.frame())
  }
  
  # Read and clean the file
  raw_lines <- readLines(mr_preg_file)
  header <- raw_lines[1]
  data_clean <- raw_lines[-1][!grepl("\\.png$", raw_lines[-1])]
  txt_clean <- c(header, data_clean)
  
  log_info("Cleaned {length(data_clean)} data lines (removed {length(raw_lines) - 1 - length(data_clean)} non-data lines)")
  
  mr_preg_df <- read.table(
    text = txt_clean,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Diagnostic information
  log_info("MR-PREG data dimensions: {nrow(mr_preg_df)} rows x {ncol(mr_preg_df)} columns")
  log_info("Available outcomes: {length(unique(mr_preg_df$Phenotype))}")
  log_info("SNPs in data: {length(unique(mr_preg_df$SNP))}")
  log_info("Mean beta (log OR): {sprintf('%.4f', mean(mr_preg_df$beta, na.rm = TRUE))}")
  log_info("Mean SE: {sprintf('%.4f', mean(mr_preg_df$se, na.rm = TRUE))}")
  
  # Format for TwoSampleMR
  mr_preg_dat <- TwoSampleMR::format_data(
    dat = mr_preg_df,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "samplesize",
    chr_col = "chr",
    pos_col = "pos_b38",
    phenotype_col = "Phenotype",
    id_col = "StudyID"
  )
  
  return(mr_preg_dat)
}

#' Helper function for calculating SE from p-value
.se_from_p <- function(beta, pval) {
  p2 <- pmin(pmax(pval/2, .Machine$double.xmin), 1 - 1e-16)
  abs(beta) / qnorm(p2, lower.tail = FALSE)
}

#' Helper function to auto-detect OR vs beta
.to_beta <- function(effect_vec, outcome_name = "") {
  med <- suppressWarnings(median(effect_vec, na.rm = TRUE))
  if (is.finite(med) && med > 0.5 && med < 1.5) {
    log_info("  {outcome_name}: Auto-detected OR format (median = {sprintf('%.3f', med)}) → converting to log OR")
    log(effect_vec)  # OR → beta
  } else {
    log_info("  {outcome_name}: Auto-detected beta format (median = {sprintf('%.3f', med)}) → using as-is")
    effect_vec      # already beta
  }
}

#' Process single PPH outcome file
#' @param base Outcome name
#' @param snps Vector of SNPs to filter for
process_single_pph <- function(base, snps) {
  path <- file.path(PATHS$data, "OUTCOME_PPH", paste0(base, ".txt"))
  
  if (!file.exists(path)) {
    log_error("PPH file not found: {path}")
    return(list(raw = data.frame(), mr = data.frame()))
  }
  
  # Expected columns in PPH files
  common_cols_pph <- c(
    "Chr", "PosB38", "Marker", "rsName", "OA", "EA",
    "EAfrq", "Effect", "P", "Phet", "I2", "nCoh"
  )
  
  # Read file
  raw <- data.table::fread(
    path,
    col.names = common_cols_pph,
    data.table = FALSE,
    showProgress = FALSE
  )
  
  # Drop any repeated header row
  if (nrow(raw) > 0 && all(tolower(as.character(raw[1, ])) == tolower(common_cols_pph))) {
    raw <- raw[-1, , drop = FALSE]
  }
  
  # Filter for selected SNPs
  df <- raw %>%
    dplyr::filter(rsName %in% snps) %>%
    dplyr::transmute(
      outcome = base,
      rsid = rsName,
      chr = gsub("^chr", "", as.character(Chr)),
      pos = as.integer(PosB38),
      a1 = EA,              # effect allele (EA)
      a2 = OA,              # other allele (OA)
      eaf = suppressWarnings(as.numeric(EAfrq)),
      beta = .to_beta(suppressWarnings(as.numeric(Effect)), base),
      pval = suppressWarnings(as.numeric(P))
    ) %>%
    dplyr::mutate(
      SE = .se_from_p(beta, pval),
      Units = "logOR",
      Gene = NA_character_,
      n = NA_integer_
    ) %>%
    as.data.frame()
  
  # Format for MR
  mr_pph <- TwoSampleMR::format_data(
    dat = df,
    type = "outcome",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "a1",
    other_allele_col = "a2",
    eaf_col = "eaf",
    pval_col = "pval",
    phenotype_col = "outcome",
    samplesize_col = "n",
    units_col = "Units",
    gene_col = "Gene",
    chr_col = "chr",
    pos_col = "pos"
  )
  
  # Ensure consistent naming
  mr_pph$outcome <- base
  mr_pph$id.outcome <- base
  
  return(list(raw = df, mr = mr_pph))
}

#' Process PPH outcomes with progress tracking
#' @param snps Vector of SNPs to filter for
#' @param pph_files Vector of PPH outcome files to process
process_pph_outcomes <- function(snps, pph_files) {
  log_info("Processing {length(pph_files)} PPH outcomes...")
  
  pb_pph <- progress::progress_bar$new(
    format = "PPH [:bar] :percent :current/:total ETA: :eta",
    total = length(pph_files),
    clear = FALSE
  )
  
  res_pph <- list()
  for (i in seq_along(pph_files)) {
    pb_pph$tick()
    res_pph[[i]] <- process_single_pph(pph_files[i], snps)
  }
  names(res_pph) <- pph_files
  
  return(res_pph)
}

#' Process single FinnGen outcome file
#' @param base Outcome name
#' @param snps Vector of SNPs to filter for
process_single_finngen <- function(base, snps) {
  path <- file.path(PATHS$data, "OUTCOME_FINNGEN", base)
  
  if (!file.exists(path)) {
    log_error("FinnGen file not found: {path}")
    return(list(raw = data.frame(), mr = data.frame()))
  }
  
  # Expected columns in FinnGen files
  common_cols_fg <- c(
    "#chrom", "pos", "ref", "alt", "rsids",
    "pval", "beta", "sebeta",
    "af_alt", "af_alt_cases", "af_alt_controls"
  )
  
  raw <- data.table::fread(
    path,
    select = common_cols_fg
  )
  
  # Filter for selected SNPs
  df <- raw %>%
    dplyr::filter(rsids %in% snps) %>%
    dplyr::transmute(
      outcome = base,
      rsid = rsids,
      chr = as.character(`#chrom`),
      pos = pos,
      a1 = alt,
      a2 = ref,
      eaf = as.numeric(af_alt),
      beta = as.numeric(beta),
      SE = as.numeric(sebeta),
      pval = as.numeric(pval),
      Units = "logOR",
      Gene = NA_character_,
      n = NA_integer_
    ) %>%
    as.data.frame()
  
  # Format for MR
  mr_fg <- TwoSampleMR::format_data(
    dat = df,
    type = "outcome",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "a1",
    other_allele_col = "a2",
    eaf_col = "eaf",
    pval_col = "pval",
    phenotype_col = "outcome",
    samplesize_col = "n",
    units_col = "Units",
    gene_col = "Gene"
  )
  
  # Ensure consistent naming
  mr_fg$outcome <- base
  mr_fg$id.outcome <- base
  
  return(list(raw = df, mr = mr_fg))
}

#' Process FinnGen outcomes with progress tracking
#' @param snps Vector of SNPs to filter for
#' @param finngen_files Vector of FinnGen outcome files to process
process_finngen_outcomes <- function(snps, finngen_files) {
  log_info("Processing {length(finngen_files)} FinnGen outcomes...")
  
  pb_fg <- progress::progress_bar$new(
    format = "FinnGen [:bar] :percent :current/:total ETA: :eta",
    total = length(finngen_files),
    clear = FALSE
  )
  
  res_fg <- list()
  for (i in seq_along(finngen_files)) {
    pb_fg$tick()
    res_fg[[i]] <- process_single_finngen(finngen_files[i], snps)
  }
  names(res_fg) <- finngen_files
  
  return(res_fg)
}

#' Create QC summary for processed outcomes
#' @param results_list List of processing results
#' @param source_name Name of the data source
create_qc_summary <- function(results_list, source_name) {
  qc_summary <- dplyr::bind_rows(lapply(results_list, `[[`, "raw")) %>%
    dplyr::group_by(outcome) %>%
    dplyr::summarize(
      nsnp = dplyr::n(),
      mean_beta = mean(beta, na.rm = TRUE),
      sd_beta = sd(beta, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(source = source_name)
  
  return(qc_summary)
}

#' Combine all outcome datasets
#' @param mr_preg_dat MR-PREG formatted data
#' @param res_pph PPH results list
#' @param res_fg FinnGen results list
combine_outcome_datasets <- function(mr_preg_dat, res_pph, res_fg) {
  log_info("Combining outcome datasets...")
  
  datasets <- list()
  
  # Add MR-PREG data if available
  if (!is.null(mr_preg_dat) && nrow(mr_preg_dat) > 0) {
    # Remove problematic columns
    mr_preg_clean <- mr_preg_dat %>%
      dplyr::select(-dplyr::any_of(c("chr.outcome", "pos.outcome")))
    datasets$mr_preg <- mr_preg_clean
    log_info("Added MR-PREG data: {nrow(mr_preg_clean)} observations")
  }
  
  # Add PPH data if available
  if (length(res_pph) > 0) {
    pph_data <- dplyr::bind_rows(lapply(res_pph, `[[`, "mr"))
    if (nrow(pph_data) > 0) {
      datasets$pph <- pph_data
      log_info("Added PPH data: {nrow(pph_data)} observations")
    }
  }
  
  # Add FinnGen data if available
  if (length(res_fg) > 0) {
    fg_data <- dplyr::bind_rows(lapply(res_fg, `[[`, "mr"))
    if (nrow(fg_data) > 0) {
      datasets$finngen <- fg_data
      log_info("Added FinnGen data: {nrow(fg_data)} observations")
    }
  }
  
  if (length(datasets) == 0) {
    log_error("No datasets available for combination")
    return(data.frame())
  }
  
  # Combine all datasets
  combined_data <- dplyr::bind_rows(datasets)
  
  log_info("Combined {length(datasets)} datasets into {nrow(combined_data)} total observations")
  
  return(combined_data)
}

#' Validate and export data with comprehensive checks
#' @param data Data frame to export
#' @param file_path Full path for export
#' @param expected_min_rows Minimum expected number of rows
#' @param required_columns Required column names
#' @param description Description for logging
validate_and_export_data <- function(data, file_path, expected_min_rows = 1, 
                                   required_columns = NULL, description = "data") {
  
  log_info("Validating {description} before export...")
  
  # Check data exists and has content
  if (is.null(data) || nrow(data) == 0) {
    log_error("No data to export for {description}")
    stop("Cannot export empty dataset")
  }
  
  # Check minimum rows
  if (nrow(data) < expected_min_rows) {
    log_warning("{description} has only {nrow(data)} rows, expected at least {expected_min_rows}")
  }
  
  # Check required columns
  if (!is.null(required_columns)) {
    missing_cols <- setdiff(required_columns, names(data))
    if (length(missing_cols) > 0) {
      log_error("Missing required columns in {description}: {paste(missing_cols, collapse = ', ')}")
      stop("Missing required columns for export")
    }
  }
  
  # Check for infinite or invalid values
  numeric_cols <- sapply(data, is.numeric)
  if (any(numeric_cols)) {
    for (col in names(data)[numeric_cols]) {
      n_infinite <- sum(is.infinite(data[[col]]))
      n_na <- sum(is.na(data[[col]]))
      if (n_infinite > 0) {
        log_warning("Column {col} has {n_infinite} infinite values")
      }
      if (n_na > 0) {
        log_info("Column {col} has {n_na} NA values")
      }
    }
  }
  
  # Export the data
  write.table(
    data,
    file = file_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  log_info("Successfully exported {description}: {basename(file_path)}")
  log_info("Exported {nrow(data)} rows and {ncol(data)} columns")
  
  return(TRUE)
}

#' Generate comprehensive coverage report
#' @param mr_preg_dat MR-PREG data
#' @param res_pph PPH results
#' @param res_fg FinnGen results
generate_coverage_report <- function(mr_preg_dat = NULL, res_pph = list(), res_fg = list()) {
  
  log_info("Generating SNP coverage report...")
  
  # Helper function to count SNPs per outcome
  count_snps <- function(res_list, label) {
    if (length(res_list) == 0) return(data.frame())
    
    dplyr::bind_rows(lapply(names(res_list), function(nm) {
      data.frame(
        outcome = nm,
        nsnp = length(unique(res_list[[nm]]$raw$rsid)),
        source = label,
        stringsAsFactors = FALSE
      )
    }))
  }
  
  # Count SNPs for each source
  coverage_detailed <- list()
  
  # MR-PREG
  if (!is.null(mr_preg_dat) && nrow(mr_preg_dat) > 0) {
    mr_preg_counts <- mr_preg_dat %>%
      dplyr::group_by(outcome) %>%
      dplyr::summarise(nsnp = dplyr::n_distinct(SNP), .groups = "drop") %>%
      dplyr::mutate(source = "MR-PREG")
    coverage_detailed$mr_preg <- mr_preg_counts
  }
  
  # PPH
  if (length(res_pph) > 0) {
    pph_counts <- count_snps(res_pph, "PPH")
    if (nrow(pph_counts) > 0) {
      coverage_detailed$pph <- pph_counts
    }
  }
  
  # FinnGen
  if (length(res_fg) > 0) {
    fg_counts <- count_snps(res_fg, "FinnGen")
    if (nrow(fg_counts) > 0) {
      coverage_detailed$finngen <- fg_counts
    }
  }
  
  # Combine all coverage data
  if (length(coverage_detailed) > 0) {
    coverage_summary <- dplyr::bind_rows(coverage_detailed) %>%
      dplyr::arrange(source, dplyr::desc(nsnp))
    
    # Summary by source
    source_summary <- coverage_summary %>%
      dplyr::group_by(source) %>%
      dplyr::summarise(
        outcomes = dplyr::n(),
        mean_nsnp = mean(nsnp),
        median_nsnp = median(nsnp),
        total_nsnp = sum(nsnp),
        .groups = "drop"
      )
  } else {
    coverage_summary <- data.frame()
    source_summary <- data.frame()
  }
  
  return(list(
    detailed = coverage_summary,
    summary = coverage_summary,
    by_source = source_summary
  ))
}

#' Generate forest plots if MR results are available
#' @param mr_results_file Path to MR results file
generate_forest_plots <- function(mr_results_file) {
  
  if (!file.exists(mr_results_file)) {
    log_warning("MR results file not found: {mr_results_file}")
    return(FALSE)
  }
  
  log_info("Loading MR results for forest plots...")
  
  # This is a placeholder - the actual forest plot generation
  # should be handled by the plotting script
  log_info("Forest plots would be generated here")
  log_info("Recommend running the dedicated plotting script instead")
  
  return(TRUE)
}

#' Validate final output quality
#' @param data Final combined dataset
validate_final_output <- function(data) {
  
  log_info("Performing final validation...")
  
  issues <- c()
  
  # Check basic structure
  if (is.null(data) || nrow(data) == 0) {
    issues <- c(issues, "Dataset is empty")
  }
  
  # Check required columns
  required_cols <- c("SNP", "beta.outcome", "se.outcome", "outcome")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    issues <- c(issues, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check for reasonable values
  if ("beta.outcome" %in% names(data)) {
    extreme_betas <- sum(abs(data$beta.outcome) > 5, na.rm = TRUE)
    if (extreme_betas > 0) {
      issues <- c(issues, paste("Extreme beta values detected:", extreme_betas))
    }
  }
  
  # Check outcome coverage
  if ("outcome" %in% names(data)) {
    n_outcomes <- length(unique(data$outcome))
    if (n_outcomes < 10) {
      issues <- c(issues, paste("Low outcome coverage:", n_outcomes, "outcomes"))
    }
  }
  
  # Check SNP coverage
  if ("SNP" %in% names(data)) {
    n_snps <- length(unique(data$SNP))
    if (n_snps < 30) {
      issues <- c(issues, paste("Low SNP coverage:", n_snps, "SNPs"))
    }
  }
  
  is_valid <- length(issues) == 0
  
  if (is_valid) {
    log_info("Final validation passed successfully")
  } else {
    log_warning("Final validation detected {length(issues)} issues")
  }
  
  return(list(
    valid = is_valid,
    issues = issues
  ))
}

log_info("Enhanced utility functions loaded successfully")