#!/usr/bin/env Rscript
###############################################################################
# run_analysis.R - endoMR-PREG Master Analysis Pipeline
# Executes complete analysis workflow with comprehensive error handling
# 
# Usage:
#   Rscript run_analysis.R                    # Run all steps
#   Rscript run_analysis.R --steps 1:5       # Run steps 1-5 only
#   Rscript run_analysis.R --continue        # Continue on errors
###############################################################################

# Load configuration and utilities
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}

source(here::here("config", "config.R"))
source(here::here("config", "utils.R"))

# Initialize pipeline logging
init_logging("run_analysis_pipeline", log_level = "INFO")

### 1. PIPELINE CONFIGURATION ###############################################

# Define the analysis workflow
PIPELINE_SCRIPTS <- list(
  "01" = list(
    file = "01_select_instruments.R",
    name = "SNP Instrument Selection", 
    description = "Select and validate genetic instruments for endometriosis",
    outputs = c("endometriosis_clumped_snps.tsv", "endo_instruments_41.txt")
  ),
  "02" = list(
    file = "02_prepare_outcomes.R", 
    name = "Outcome Data Preparation",
    description = "Process pregnancy outcome data from multiple sources",
    outputs = c("stu_out_dat.txt", "combined_outcomes_data.csv")
  ),
  "03" = list(
    file = "03_harmonise_data.R",
    name = "Data Harmonisation", 
    description = "Harmonise exposure and outcome data for MR analysis",
    outputs = c("harmonised_rahmioglu_bpo.csv")
  ),
  "04" = list(
    file = "04_mendelian_randomization.R",
    name = "Mendelian Randomization Analysis",
    description = "Main MR analyses with sensitivity tests",
    outputs = c("ivw_results.csv", "all_mr_methods.csv", "heterogeneity_results.csv")
  ),
  "05" = list(
    file = "05_generate_tables.R", 
    name = "Results Tables Generation",
    description = "Generate publication-ready tables",
    outputs = c("main_results_table.csv", "supplementary_tables.csv")
  ),
  "06" = list(
    file = "06_create_plots.R",
    name = "Visualization Creation", 
    description = "Create scatter plots and diagnostic plots",
    outputs = c("scatter_plots/", "diagnostic_plots/")
  ),
  "07" = list(
    file = "07_sensitivity_analysis.R",
    name = "Sensitivity Analysis",
    description = "Leave-one-cohort-out and additional robustness tests", 
    outputs = c("leave_one_cohort_results.csv", "sensitivity_summary.csv")
  ),
  "08" = list(
    file = "08_forest_plots.R",
    name = "Forest Plot Generation",
    description = "Create publication-quality forest plots",
    outputs = c("forest_endoMR-PREG_IVW.png", "forest_endoMR-PREG_multiMethod.png")
  ),
  "09" = list(
    file = "09_supplementary_tables.R", 
    name = "Supplementary Materials",
    description = "Generate supplementary tables and materials",
    outputs = c("Supplementary_Table_S2.csv", "supplementary_materials/")
  )
)

# Pipeline configuration
PIPELINE_CONFIG <- list(
  validate_outputs = TRUE,
  create_backups = TRUE,
  generate_report = TRUE,
  parallel_compatible_steps = c("06", "08"),  # Steps that can run in parallel
  critical_outputs = c("harmonised_rahmioglu_bpo.csv", "ivw_results.csv")
)

### 2. COMMAND LINE ARGUMENT PARSING ########################################
parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  config <- list(
    steps = 1:9,
    continue_on_error = FALSE,
    skip_validation = FALSE,
    force_rerun = FALSE,
    verbose = TRUE
  )
  
  if ("--help" %in% args || "-h" %in% args) {
    cat("endoMR-PREG Analysis Pipeline\n\n")
    cat("Usage: Rscript run_analysis.R [options]\n\n")
    cat("Options:\n")
    cat("  --steps X:Y         Run steps X through Y (default: 1:9)\n")
    cat("  --continue          Continue on errors (default: stop on first error)\n")
    cat("  --skip-validation   Skip output validation\n") 
    cat("  --force             Force rerun even if outputs exist\n")
    cat("  --quiet             Reduce output verbosity\n")
    cat("  --help, -h          Show this help\n\n")
    cat("Available steps:\n")
    for (step in names(PIPELINE_SCRIPTS)) {
      script_info <- PIPELINE_SCRIPTS[[step]]
      cat(sprintf("  %s: %s\n", step, script_info$name))
    }
    quit(status = 0)
  }
  
  if ("--continue" %in% args) config$continue_on_error <- TRUE
  if ("--skip-validation" %in% args) config$skip_validation <- TRUE
  if ("--force" %in% args) config$force_rerun <- TRUE
  if ("--quiet" %in% args) config$verbose <- FALSE
  
  # Parse steps argument
  steps_arg <- grep("--steps", args, value = TRUE)
  if (length(steps_arg) > 0) {
    steps_value <- sub("--steps\\s*", "", steps_arg)
    if (grepl(":", steps_value)) {
      step_range <- as.numeric(strsplit(steps_value, ":")[[1]])
      config$steps <- step_range[1]:step_range[2]
    } else {
      config$steps <- as.numeric(strsplit(steps_value, ",")[[1]])
    }
  }
  
  return(config)
}

### 3. UTILITY FUNCTIONS ####################################################
validate_step_outputs <- function(step_id, expected_outputs) {
  if (!PIPELINE_CONFIG$validate_outputs) return(TRUE)
  
  log_info("Validating outputs for step", step_id)
  missing_outputs <- c()
  
  for (output in expected_outputs) {
    # Check in multiple possible locations
    possible_paths <- c(
      file.path(RESULTS_DIR, output),
      file.path(PLOTS_DIR, output), 
      file.path(DATA_DIR, output),
      file.path(PROJECT_ROOT, output)
    )
    
    found <- any(sapply(possible_paths, function(p) file.exists(p) || dir.exists(p)))
    
    if (!found) {
      missing_outputs <- c(missing_outputs, output)
    }
  }
  
  if (length(missing_outputs) > 0) {
    log_warn("Missing expected outputs:", paste(missing_outputs, collapse = ", "))
    return(FALSE)
  }
  
  log_info("All expected outputs validated for step", step_id)
  return(TRUE)
}

create_step_backup <- function(step_id) {
  if (!PIPELINE_CONFIG$create_backups) return(TRUE)
  
  backup_dir <- file.path(PROJECT_ROOT, "backups", Sys.Date(), paste0("step_", step_id))
  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Backup key results
  key_dirs <- c(RESULTS_DIR, PLOTS_DIR)
  for (dir_path in key_dirs) {
    if (dir.exists(dir_path)) {
      backup_target <- file.path(backup_dir, basename(dir_path))
      file.copy(dir_path, backup_target, recursive = TRUE)
    }
  }
  
  log_info("Created backup for step", step_id, "in", backup_dir)
  return(TRUE)
}

check_step_prerequisites <- function(step_id) {
  # Define prerequisites for each step
  prerequisites <- list(
    "02" = "01",  # Step 2 needs step 1 outputs
    "03" = c("01", "02"),  # Step 3 needs steps 1 and 2
    "04" = "03",  # MR analysis needs harmonised data
    "05" = "04",  # Tables need MR results
    "06" = "04",  # Plots need MR results  
    "07" = "04",  # Sensitivity needs MR results
    "08" = "04",  # Forest plots need MR results
    "09" = c("04", "05")   # Supplementary needs MR and tables
  )
  
  if (!step_id %in% names(prerequisites)) return(TRUE)
  
  required_steps <- prerequisites[[step_id]]
  for (req_step in required_steps) {
    script_info <- PIPELINE_SCRIPTS[[req_step]]
    outputs_exist <- validate_step_outputs(req_step, script_info$outputs)
    
    if (!outputs_exist) {
      log_error("Step", step_id, "requires outputs from step", req_step, "which are missing")
      return(FALSE)
    }
  }
  
  return(TRUE)
}

### 4. MAIN PIPELINE FUNCTION ###############################################
run_pipeline <- function(config) {
  
  log_info("=== Starting endoMR-PREG Analysis Pipeline ===")
  log_info("Configuration:")
  log_info("- Steps to run:", paste(config$steps, collapse = ", "))
  log_info("- Continue on error:", config$continue_on_error)
  log_info("- Force rerun:", config$force_rerun)
  
  # Initialize pipeline tracking
  pipeline_start_time <- Sys.time()
  step_results <- list()
  
  # Validate environment
  log_info("Validating analysis environment...")
  validate_data_directories()
  
  # Run each requested step
  for (step_num in config$steps) {
    step_id <- sprintf("%02d", step_num)
    
    if (!step_id %in% names(PIPELINE_SCRIPTS)) {
      log_warn("Step", step_id, "not found in pipeline, skipping")
      next
    }
    
    script_info <- PIPELINE_SCRIPTS[[step_id]]
    script_path <- file.path(SCRIPTS_DIR, script_info$file)
    
    log_info("=== Step", step_id, ":", script_info$name, "===")
    log_info("Description:", script_info$description)
    
    # Check if step can be skipped (outputs exist and not forcing rerun)
    if (!config$force_rerun) {
      outputs_exist <- validate_step_outputs(step_id, script_info$outputs)
      if (outputs_exist) {
        log_info("Step", step_id, "outputs already exist, skipping (use --force to rerun)")
        step_results[[step_id]] <- list(status = "skipped", duration = 0)
        next
      }
    }
    
    # Check prerequisites
    if (!check_step_prerequisites(step_id)) {
      log_error("Prerequisites not met for step", step_id)
      if (!config$continue_on_error) {
        stop("Pipeline halted due to missing prerequisites")
      }
      step_results[[step_id]] <- list(status = "failed_prerequisites", duration = 0)
      next
    }
    
    # Create backup before running
    create_step_backup(step_id)
    
    # Run the step
    step_start_time <- Sys.time()
    step_status <- tryCatch({
      
      if (!file.exists(script_path)) {
        stop("Script file not found: ", script_path)
      }
      
      log_info("Executing:", script_info$file)
      source(script_path)
      
      # Validate outputs were created
      if (!config$skip_validation) {
        outputs_created <- validate_step_outputs(step_id, script_info$outputs)
        if (!outputs_created) {
          warning("Some expected outputs were not created")
        }
      }
      
      "completed"
      
    }, error = function(e) {
      log_error("Step", step_id, "failed:", e$message)
      if (!config$continue_on_error) {
        stop("Pipeline halted at step ", step_id, ": ", e$message)
      }
      "failed"
    })
    
    step_end_time <- Sys.time()
    step_duration <- as.numeric(step_end_time - step_start_time, units = "mins")
    
    step_results[[step_id]] <- list(
      status = step_status,
      duration = step_duration,
      start_time = step_start_time,
      end_time = step_end_time
    )
    
    log_info("Step", step_id, "completed with status:", step_status)
    log_info("Duration:", round(step_duration, 2), "minutes")
  }
  
  # Generate pipeline report
  if (PIPELINE_CONFIG$generate_report) {
    generate_pipeline_report(step_results, pipeline_start_time)
  }
  
  log_info("=== Pipeline execution completed ===")
  return(step_results)
}

### 5. PIPELINE REPORT GENERATION ###########################################
generate_pipeline_report <- function(step_results, start_time) {
  end_time <- Sys.time()
  total_duration <- as.numeric(end_time - start_time, units = "mins")
  
  # Create report data
  report_data <- data.frame(
    step = names(step_results),
    status = sapply(step_results, `[[`, "status"),
    duration_mins = round(sapply(step_results, `[[`, "duration"), 2),
    stringsAsFactors = FALSE
  )
  
  # Summary statistics
  n_completed <- sum(report_data$status == "completed")
  n_failed <- sum(report_data$status == "failed")
  n_skipped <- sum(report_data$status == "skipped")
  
  # Generate report
  report_path <- file.path(RESULTS_DIR, paste0("pipeline_report_", Sys.Date(), ".txt"))
  
  cat("endoMR-PREG Analysis Pipeline Report\n", file = report_path)
  cat("=====================================\n\n", file = report_path, append = TRUE)
  cat("Execution Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = report_path, append = TRUE)
  cat("Total Duration:", round(total_duration, 2), "minutes\n\n", file = report_path, append = TRUE)
  
  cat("Summary:\n", file = report_path, append = TRUE)
  cat("- Steps completed:", n_completed, "\n", file = report_path, append = TRUE)
  cat("- Steps failed:", n_failed, "\n", file = report_path, append = TRUE)
  cat("- Steps skipped:", n_skipped, "\n\n", file = report_path, append = TRUE)
  
  cat("Step Details:\n", file = report_path, append = TRUE)
  write.table(report_data, file = report_path, append = TRUE, sep = "\t", 
              row.names = FALSE, quote = FALSE)
  
  log_info("Pipeline report saved:", report_path)
}

### 6. MAIN EXECUTION ########################################################
if (!interactive()) {
  # Parse command line arguments and run pipeline
  config <- parse_arguments()
  
  tryCatch({
    results <- run_pipeline(config)
    
    # Exit with appropriate code
    failed_steps <- sapply(results, function(x) x$status == "failed")
    if (any(failed_steps)) {
      log_warn("Pipeline completed with failures")
      quit(status = 1)
    } else {
      log_info("Pipeline completed successfully")
      quit(status = 0)
    }
    
  }, error = function(e) {
    log_error("Pipeline execution failed:", e$message)
    quit(status = 1)
  })
}
