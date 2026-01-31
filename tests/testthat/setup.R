# Test Setup for fmri.pipeline
# =============================
# This file is loaded before tests run. It should contain only:
# - Package imports needed for testing
# - Global test configuration
# - References to helper files (loaded automatically by testthat)
#
# Heavy object builders and integration test utilities have been moved to:
# - helpers-mock_factories.R: Lightweight mock object creators for unit tests
# - helpers-integration.R: Real data loaders for integration tests (optional)

library(data.table)
library(dplyr)
library(fmri.pipeline)

# ==============================================================================
# Global Test Configuration
# ==============================================================================

# Path to test data directory (for integration tests that need real data)
# This should only be used by integration tests, not unit tests
TEST_DATA_DIR <- Sys.getenv("FMRI_PIPELINE_TEST_DATA", 
                            unset = "/proj/mnhallqlab/projects/fmri.pipeline_test_data")

# Default settings for tests
TEST_TR <- 1.0
TEST_DROP_VOLUMES <- 2L

# ==============================================================================
# Simple Data Accessors (for tests that need real data files)
# ==============================================================================

#' Get the test data directory path
#' @return Path to test data directory
get_test_data_dir <- function() {
  if (!dir.exists(TEST_DATA_DIR)) {
    skip("Test data directory not available")
  }
  TEST_DATA_DIR
}

#' Check if integration test data is available
#' @return TRUE if test data directory exists
has_integration_test_data <- function() {
  dir.exists(TEST_DATA_DIR)
}

# ==============================================================================
# Legacy Compatibility Functions
# ==============================================================================
# These functions maintain backward compatibility with existing tests
# but should be migrated to use mock factories for true unit tests

#' Provide a minimal instantiation of the gpa list
#' @return a minimal gpa list
#' @note Prefer create_mock_gpa_minimal() from helpers-mock_factories.R
get_gpa_minimal <- function() {
  # Use the mock factory for consistency
  create_mock_gpa_minimal()
}

# ==============================================================================
# Integration Test Utilities
# ==============================================================================
# These functions are for integration tests that need to work with real data
# They should NOT be used in unit tests

#' Build a GPA object from real test data files (integration tests only)
#'
#' @param analysis_name Name for the analysis
#' @param test_data_base_dir Base directory containing test data
#' @param trial_data_file Filename of trial data
#' @param run_data_file Filename of run data
#' @param subject_data_file Filename of subject data
#' @param ... Additional arguments passed to setup_glm_pipeline
#' @return A gpa object
build_integration_gpa <- function(
    analysis_name = "integration_test",
    test_data_base_dir = TEST_DATA_DIR,
    trial_data_file = "sample_trial_data.csv.gz",
    run_data_file = "sample_run_data.csv",
    subject_data_file = "sample_subject_data.csv",
    scheduler = "local",
    drop_volumes = TEST_DROP_VOLUMES,
    tr = TEST_TR,
    ...
) {
  if (!has_integration_test_data()) {
    stop("Integration test data not available at: ", test_data_base_dir)
  }
  
  trial_df <- read.csv(file.path(test_data_base_dir, trial_data_file))
  run_df <- read.csv(file.path(test_data_base_dir, run_data_file))
  subj_df <- read.csv(file.path(test_data_base_dir, subject_data_file))
  
  gpa <- setup_glm_pipeline(
    analysis_name = analysis_name,
    scheduler = scheduler,
    trial_data = trial_df,
    run_data = run_df,
    subject_data = subj_df,
    n_expected_runs = 8,
    tr = tr,
    drop_volumes = drop_volumes,
    l1_models = NULL, 
    l2_models = NULL, 
    l3_models = NULL,
    confound_settings = list(
      motion_params_file = "motion.par",
      confound_input_file = "nuisance_regressors.txt",
      confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
      l1_confound_regressors = c("csf", "dcsf", "wm", "dwm")
    ),
    ...
  )
  
  return(gpa)
}

# ==============================================================================
# Skip Helpers
# ==============================================================================

#' Skip test if running in CI without test data
skip_without_test_data <- function() {
  if (!has_integration_test_data()) {
    skip("Integration test data not available")
  }
}

#' Skip test if not running interactively (for manual/exploratory tests)
skip_if_not_interactive <- function() {
  if (!interactive()) {
    skip("Test requires interactive session")
  }
}

#' Skip test if on SLURM cluster but not in a job (for cluster-specific tests)
skip_without_slurm <- function() {
  if (Sys.getenv("SLURM_JOB_ID") == "") {
    skip("Test requires SLURM job environment")
  }
}
