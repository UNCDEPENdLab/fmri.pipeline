# Unit tests for glm_pipeline_arguments (GPA) object
# ===================================================
# Tests for the core gpa object structure and methods

library(testthat)

# ==============================================================================
# Tests for GPA object structure
# ==============================================================================

test_that("create_mock_gpa returns valid gpa object", {
  gpa <- create_mock_gpa()
  
  expect_s3_class(gpa, "glm_pipeline_arguments")
  expect_s3_class(gpa, "list")
})

test_that("mock gpa has required fields", {
  gpa <- create_mock_gpa()
  
  required_fields <- c("analysis_name", "tr", "subject_data", "run_data", 
                       "trial_data", "output_directory", "l1_models", 
                       "l2_models", "l3_models")
  
  for (field in required_fields) {
    expect_true(field %in% names(gpa), 
                info = sprintf("Missing required field: %s", field))
  }
})

test_that("mock gpa subject_data has required columns", {
  gpa <- create_mock_gpa(n_subjects = 5)
  
  expect_s3_class(gpa$subject_data, "data.frame")
  expect_true("id" %in% names(gpa$subject_data))
  expect_equal(nrow(gpa$subject_data), 5)
})

test_that("mock gpa run_data has required columns", {
  gpa <- create_mock_gpa(n_subjects = 3, n_runs = 4)
  
  expect_s3_class(gpa$run_data, "data.frame")
  expect_true(all(c("id", "run_number", "run_volumes") %in% names(gpa$run_data)))
  expect_equal(nrow(gpa$run_data), 3 * 4)  # n_subjects * n_runs
})

test_that("mock gpa trial_data has required columns", {
  gpa <- create_mock_gpa(n_subjects = 2, n_runs = 2, n_trials = 10)
  
  expect_s3_class(gpa$trial_data, "data.frame")
  expect_true(all(c("id", "run_number", "trial", "onset") %in% names(gpa$trial_data)))
  expect_equal(nrow(gpa$trial_data), 2 * 2 * 10)  # n_subjects * n_runs * n_trials
})

# ==============================================================================
# Tests for GPA summary method
# ==============================================================================

test_that("summary.glm_pipeline_arguments prints without error", {
  gpa <- create_mock_gpa()
  
  # Should not throw an error
  expect_output(summary(gpa))
})

test_that("summary works with minimal gpa", {
  gpa <- create_mock_gpa_minimal()
  
  expect_output(summary(gpa))
})

test_that("summary works with l1_models included", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  
  expect_output(summary(gpa))
})

test_that("summary works with all model levels included", {
  gpa <- create_mock_gpa(
    include_l1_models = TRUE,
    include_l2_models = TRUE,
    include_l3_models = TRUE
  )
  
  expect_output(summary(gpa))
})

# ==============================================================================
# Tests for GPA parameter validation
# ==============================================================================

test_that("gpa preserves tr parameter", {
  gpa <- create_mock_gpa(tr = 2.5)
  
  expect_equal(gpa$tr, 2.5)
})

test_that("gpa preserves analysis_name", {
  gpa <- create_mock_gpa(analysis_name = "my_test_analysis")
  
  expect_equal(gpa$analysis_name, "my_test_analysis")
})

test_that("gpa preserves output_directory", {
  test_dir <- tempdir()
  gpa <- create_mock_gpa(output_directory = test_dir)
  
  expect_equal(gpa$output_directory, test_dir)
})

# ==============================================================================
# Tests for GPA with different configurations
# ==============================================================================

test_that("gpa handles single subject", {
  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2)
  
  expect_equal(nrow(gpa$subject_data), 1)
  expect_equal(nrow(gpa$run_data), 2)
})

test_that("gpa handles many subjects", {
  gpa <- create_mock_gpa(n_subjects = 50, n_runs = 1)
  
  expect_equal(nrow(gpa$subject_data), 50)
  expect_equal(nrow(gpa$run_data), 50)
})

test_that("gpa handles many runs", {
  gpa <- create_mock_gpa(n_subjects = 2, n_runs = 10)
  
  expect_equal(nrow(gpa$run_data), 20)
})

# ==============================================================================
# Tests for L1 models structure
# ==============================================================================

test_that("l1_models has correct structure when included", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  
  expect_true(!is.null(gpa$l1_models))
  expect_true("signals" %in% names(gpa$l1_models))
  expect_true("models" %in% names(gpa$l1_models))
})

test_that("l1_models signals have required fields", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  
  for (sig_name in names(gpa$l1_models$signals)) {
    sig <- gpa$l1_models$signals[[sig_name]]
    expect_true("name" %in% names(sig), 
                info = sprintf("Signal %s missing 'name' field", sig_name))
    expect_true("event" %in% names(sig),
                info = sprintf("Signal %s missing 'event' field", sig_name))
  }
})

# ==============================================================================
# Tests for L2/L3 models structure
# ==============================================================================

test_that("l2_models has correct structure when included", {
  gpa <- create_mock_gpa(include_l2_models = TRUE)
  
  expect_true(!is.null(gpa$l2_models))
  expect_true("models" %in% names(gpa$l2_models))
})

test_that("l3_models has correct structure when included", {
  gpa <- create_mock_gpa(include_l3_models = TRUE)
  
  expect_true(!is.null(gpa$l3_models))
  expect_true("models" %in% names(gpa$l3_models))
})

# ==============================================================================
# Tests for NULL handling
# ==============================================================================

test_that("gpa handles NULL l1_models", {
  gpa <- create_mock_gpa(include_l1_models = FALSE)
  
  expect_null(gpa$l1_models)
})

test_that("gpa handles NULL l2_models", {
  gpa <- create_mock_gpa(include_l2_models = FALSE)
  
  expect_null(gpa$l2_models)
})

test_that("gpa handles NULL l3_models", {
  gpa <- create_mock_gpa(include_l3_models = FALSE)
  
  expect_null(gpa$l3_models)
})
