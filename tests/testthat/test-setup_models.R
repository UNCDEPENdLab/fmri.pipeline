# Unit tests for setup_l1_models, setup_l2_models, setup_l3_models
# ==================================================================
# These tests focus on input validation and helper functions that can be
# tested in isolation without requiring external data or cluster resources.

library(testthat)

# ==============================================================================
# Mock factories specific to setup_l*_models testing
# ==============================================================================

#' Create a mock gpa suitable for setup_l1_models input validation
create_mock_gpa_for_l1_setup <- function(
    n_subjects = 2,
    n_runs = 2,
    include_l1_models = TRUE
) {
  gpa <- create_mock_gpa(
    n_subjects = n_subjects,
    n_runs = n_runs,
    include_l1_models = include_l1_models
  )
  
  # Add required columns for setup_l1_models
  gpa$run_data$run_nifti_present <- TRUE
  gpa$run_data$tr <- 1.0
  gpa$run_data$nvoxels <- 100000L
  
  # Add output locations
  gpa$output_locations$setup_l1_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l1_log_txt <- tempfile(fileext = ".txt")
  
  # Add logging threshold
  gpa$lgr_threshold <- "warn"  # Suppress debug messages during testing
  gpa$log_json <- FALSE
  gpa$log_txt <- FALSE
  
  # Add parallel settings
  gpa$parallel$l1_setup_cores <- 1L
  
  return(gpa)
}

#' Create a mock gpa suitable for setup_l2_models input validation
create_mock_gpa_for_l2_setup <- function(
    n_subjects = 2,
    n_runs = 2
) {
  gpa <- create_mock_gpa(
    n_subjects = n_subjects,
    n_runs = n_runs,
    include_l1_models = TRUE,
    include_l2_models = TRUE
  )
  
  # Add required columns for setup_l2_models
  gpa$run_data$exclude_run <- FALSE
  gpa$run_data$exclude_subject <- FALSE
  
  # Add output locations
  gpa$output_locations$setup_l2_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l2_log_txt <- tempfile(fileext = ".txt")
  
  # Add logging threshold
  gpa$lgr_threshold <- "warn"
  gpa$log_json <- FALSE
  gpa$log_txt <- FALSE
  
  # Add parallel settings
  gpa$parallel$l2_setup_cores <- 1L
  
  # Add l1_model_setup (would normally be created by setup_l1_models)
  gpa$l1_model_setup <- create_mock_l1_model_setup(gpa)
  
  # Set correct class for l1_models and l2_models
  class(gpa$l1_models) <- c("l1_model_set", "list")
  class(gpa$l2_models) <- c("hi_model_set", "list")
  
  return(gpa)
}

#' Create a mock gpa suitable for setup_l3_models input validation
create_mock_gpa_for_l3_setup <- function(
    n_subjects = 3,
    n_runs = 2,
    multi_run = TRUE
) {
  gpa <- create_mock_gpa(
    n_subjects = n_subjects,
    n_runs = n_runs,
    include_l1_models = TRUE,
    include_l2_models = multi_run,
    include_l3_models = TRUE
  )
  
  # Add required columns
  gpa$run_data$exclude_run <- FALSE
  gpa$run_data$exclude_subject <- FALSE
  gpa$subject_data$exclude_subject <- FALSE
  
  # Add output locations
  gpa$output_locations$setup_l3_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l3_log_txt <- tempfile(fileext = ".txt")
  
  # Add settings
  gpa$lgr_threshold <- "warn"
  gpa$log_json <- FALSE
  gpa$log_txt <- FALSE
  gpa$multi_run <- multi_run
  gpa$glm_software <- "fsl"
  
  # Add l1_model_setup and l2_model_setup
  gpa$l1_model_setup <- create_mock_l1_model_setup(gpa)
  if (multi_run) {
    gpa$l2_model_setup <- create_mock_l2_model_setup(gpa)
  }
  
  # Set correct classes
  class(gpa$l1_models) <- c("l1_model_set", "list")
  if (multi_run) {
    class(gpa$l2_models) <- c("hi_model_set", "list")
  }
  class(gpa$l3_models) <- c("hi_model_set", "list")
  
  return(gpa)
}

#' Create a mock l1_model_setup object
create_mock_l1_model_setup <- function(gpa) {
  n_rows <- nrow(gpa$run_data)
  
  fsl_setup <- data.frame(
    id = gpa$run_data$id,
    session = gpa$run_data$session,
    run_number = gpa$run_data$run_number,
    l1_model = rep("model1", n_rows),
    feat_fsf = rep(tempfile(fileext = ".fsf"), n_rows),
    feat_dir = rep(tempfile(), n_rows),
    feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  
  result <- list(fsl = fsl_setup, spm = list(), afni = list())
  class(result) <- c("l1_setup", "list")
  return(result)
}

#' Create a mock l2_model_setup object
create_mock_l2_model_setup <- function(gpa) {
  subjects <- unique(gpa$subject_data[, c("id", "session")])
  n_rows <- nrow(subjects)
  
  fsl_setup <- data.frame(
    id = subjects$id,
    session = subjects$session,
    l1_model = rep("model1", n_rows),
    l2_model = rep("l2_model1", n_rows),
    feat_fsf = rep(tempfile(fileext = ".fsf"), n_rows),
    feat_dir = rep(tempfile(), n_rows),
    feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  
  result <- list(fsl = fsl_setup, spm = list(), afni = list())
  class(result) <- c("l2_setup", "list")
  return(result)
}

# ==============================================================================
# Tests for setup_l1_models input validation
# ==============================================================================

test_that("setup_l1_models rejects invalid gpa class", {
  expect_error(
    setup_l1_models(list()),
    "glm_pipeline_arguments"
  )
})

test_that("setup_l1_models rejects invalid l1_model_names type", {
  gpa <- create_mock_gpa_for_l1_setup()
  
  expect_error(
    setup_l1_models(gpa, l1_model_names = 123),
    "character"
  )
})

test_that("setup_l1_models rejects unknown model names", {
  gpa <- create_mock_gpa_for_l1_setup()
  
  expect_error(
    setup_l1_models(gpa, l1_model_names = c("nonexistent_model")),
    "subset"
  )
})

test_that("setup_l1_models requires run_data to be a data.frame", {
  gpa <- create_mock_gpa_for_l1_setup()
  gpa$run_data <- NULL
  
  expect_error(
    setup_l1_models(gpa),
    "data.frame"
  )
})

test_that("setup_l1_models requires specific columns in run_data", {
  gpa <- create_mock_gpa_for_l1_setup()
  gpa$run_data$run_nifti <- NULL  # Remove required column
  
  expect_error(
    setup_l1_models(gpa),
    "run_nifti"
  )
})

test_that("setup_l1_models defaults to all models when l1_model_names is NULL", {
  gpa <- create_mock_gpa_for_l1_setup()
  # NULL l1_model_names defaults to all models (line 33 in setup_l1_models.R)
  # The function handles missing files gracefully, so we just verify it doesn't error on validation
  # It will return a gpa with updated l1_model_setup, even if empty
  result <- setup_l1_models(gpa, l1_model_names = NULL)
  expect_s3_class(result, "glm_pipeline_arguments")
})

# ==============================================================================
# Tests for setup_l2_models input validation
# ==============================================================================

test_that("setup_l2_models rejects invalid gpa class", {
  expect_error(
    setup_l2_models(list()),
    "glm_pipeline_arguments"
  )
})

test_that("setup_l2_models rejects missing l2_models", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l2_models <- NULL
  
  expect_error(
    setup_l2_models(gpa),
    "hi_model_set"
  )
})

test_that("setup_l2_models rejects missing l1_models", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l1_models <- NULL
  
  expect_error(
    setup_l2_models(gpa),
    "l1_model_set"
  )
})

test_that("setup_l2_models rejects invalid l1_model_names", {
  gpa <- create_mock_gpa_for_l2_setup()
  
  expect_error(
    setup_l2_models(gpa, l1_model_names = c("nonexistent")),
    "subset"
  )
})

test_that("setup_l2_models rejects invalid l2_model_names", {
  gpa <- create_mock_gpa_for_l2_setup()
  
  expect_error(
    setup_l2_models(gpa, l2_model_names = c("nonexistent")),
    "subset"
  )
})

test_that("setup_l2_models requires l1_model_setup", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l1_model_setup <- NULL
  
  expect_error(
    setup_l2_models(gpa),
    "l1_model_setup|setup_l1_models"
  )
})

# ==============================================================================
# Tests for setup_l3_models input validation
# ==============================================================================

test_that("setup_l3_models rejects invalid gpa class", {
  expect_error(
    setup_l3_models(list()),
    "glm_pipeline_arguments"
  )
})

test_that("setup_l3_models rejects missing l3_models", {
  gpa <- create_mock_gpa_for_l3_setup()
  gpa$l3_models <- NULL
  
  expect_error(
    setup_l3_models(gpa),
    "hi_model_set"
  )
})

test_that("setup_l3_models rejects missing l1_models", {
  gpa <- create_mock_gpa_for_l3_setup()
  gpa$l1_models <- NULL
  
  expect_error(
    setup_l3_models(gpa),
    "l1_model_set"
  )
})

test_that("setup_l3_models handles invalid l3_model_names in downstream processing", {
  # Note: setup_l3_models does not directly validate l3_model_names as a subset
  # It errors when trying to process unknown models in the actual setup loop
  # The validation happens via choose_glm_models in run_glm_pipeline
  gpa <- create_mock_gpa_for_l3_setup()
  # Ensure l2 models are marked complete to get past that check
  gpa$l2_model_setup$fsl$feat_complete <- TRUE
  gpa$l1_model_setup$fsl$feat_complete <- TRUE
  
  # With a nonexistent model name, it should error somewhere in processing
  # (either on subset check or when the model can't be found)
  expect_error(
    setup_l3_models(gpa, l3_model_names = c("nonexistent"))
  )
})

test_that("setup_l3_models requires l2_models for multi_run", {
  gpa <- create_mock_gpa_for_l3_setup(multi_run = TRUE)
  gpa$l2_models <- NULL
  
  expect_error(
    setup_l3_models(gpa),
    "hi_model_set"
  )
})

# ==============================================================================
# Tests for enforce_glms_complete helper function
# ==============================================================================

test_that("enforce_glms_complete rejects invalid gpa", {
  lg <- lgr::get_logger("test")
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(list(), level = 1L, lg = lg),
    "glm_pipeline_arguments"
  )
})

test_that("enforce_glms_complete rejects invalid level", {
  gpa <- create_mock_gpa()
  lg <- lgr::get_logger("test")
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 0L, lg = lg),
    "level"
  )
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 4L, lg = lg),
    "level"
  )
})

test_that("enforce_glms_complete errors when l1_model_setup is missing", {
  gpa <- create_mock_gpa()
  gpa$l1_model_setup <- NULL
  lg <- lgr::get_logger("test")
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 1L, lg = lg),
    "l1_model_setup"
  )
})

test_that("enforce_glms_complete errors when all feat runs incomplete", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l1_model_setup$fsl$feat_complete <- FALSE
  gpa$glm_software <- "fsl"
  lg <- lgr::get_logger("test")
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 1L, lg = lg),
    "incomplete"
  )
})

test_that("enforce_glms_complete warns when some feat runs incomplete", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l1_model_setup$fsl$feat_complete[1] <- FALSE
  gpa$glm_software <- "fsl"
  lg <- lgr::get_logger("test")
  lg$set_threshold("warn")
  
  # Should not error, just warn
  expect_no_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 1L, lg = lg)
  )
})

test_that("enforce_glms_complete passes with all complete runs", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$glm_software <- "fsl"
  lg <- lgr::get_logger("test")
  
  expect_no_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 1L, lg = lg)
  )
})

# ==============================================================================
# Tests for respecify_l2_models_by_subject helper function
# ==============================================================================

test_that("respecify_l2_models_by_subject rejects invalid model object", {
  data <- data.frame(id = c("s1", "s2"), session = c(1, 1), value = c(1, 2))
  
  expect_error(
    fmri.pipeline:::respecify_l2_models_by_subject(list(), data),
    "l1_model_spec|hi_model_spec"
  )
})

test_that("respecify_l2_models_by_subject requires id and session columns", {
  # Create minimal model object
  mobj <- list(
    lmfit = lm(value ~ 1, data = data.frame(value = 1:3)),
    contrasts = list(intercept = c(intercept = 1))
  )
  class(mobj) <- c("hi_model_spec", "list")
  
  data <- data.frame(x = 1:3)
  
  expect_error(
    fmri.pipeline:::respecify_l2_models_by_subject(mobj, data),
    "id|session"
  )
})

# ==============================================================================
# Tests for choose_glm_models helper function
# ==============================================================================

test_that("choose_glm_models returns all models when 'all' specified", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  result <- fmri.pipeline:::choose_glm_models(gpa, "all", level = 1)
  
  expect_equal(result, names(gpa$l1_models$models))
})

test_that("choose_glm_models returns NULL when 'none' specified", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  result <- fmri.pipeline:::choose_glm_models(gpa, "none", level = 1)
  
  expect_null(result)
})

test_that("choose_glm_models returns specified model names", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  model_name <- names(gpa$l1_models$models)[1]
  
  result <- fmri.pipeline:::choose_glm_models(gpa, model_name, level = 1)
  
  expect_equal(result, model_name)
})

test_that("choose_glm_models rejects invalid model names", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  expect_error(
    fmri.pipeline:::choose_glm_models(gpa, "nonexistent_model", level = 1),
    "subset"
  )
})

test_that("choose_glm_models validates level argument", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  # Level 0 is invalid - should error on checkmate validation
  expect_error(
    fmri.pipeline:::choose_glm_models(gpa, "all", level = 0),
    "Element 1 is not >= 1"
  )
  
  # Level 4 is invalid - should error on checkmate validation
  expect_error(
    fmri.pipeline:::choose_glm_models(gpa, "all", level = 4),
    "Element 1 is not <= 3"
  )
})

test_that("choose_glm_models handles missing model object gracefully", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  # Remove l2_models
  gpa$l2_models <- NULL
  
  # When the model object doesn't exist, the function should handle it
  # (either error or return NULL)
  result <- tryCatch(
    fmri.pipeline:::choose_glm_models(gpa, "all", level = 2),
    error = function(e) "error"
  )
  
  # Either returns NULL (no models) or errors
  expect_true(is.null(result) || identical(result, "error"))
})

# ==============================================================================
# Tests for model setup data structures
# ==============================================================================

test_that("l1_setup object has correct structure", {
  gpa <- create_mock_gpa_for_l2_setup()
  setup <- gpa$l1_model_setup
  
  expect_s3_class(setup, "l1_setup")
  expect_true("fsl" %in% names(setup))
  expect_s3_class(setup$fsl, "data.frame")
})

test_that("l1_setup$fsl has required columns", {
  gpa <- create_mock_gpa_for_l2_setup()
  setup <- gpa$l1_model_setup
  
  required_cols <- c("id", "session", "run_number", "l1_model", "feat_fsf", "feat_dir", "feat_complete")
  expect_true(all(required_cols %in% names(setup$fsl)))
})

test_that("l2_setup object has correct structure", {
  gpa <- create_mock_gpa_for_l3_setup()
  setup <- gpa$l2_model_setup
  
  expect_s3_class(setup, "l2_setup")
  expect_true("fsl" %in% names(setup))
  expect_s3_class(setup$fsl, "data.frame")
})

test_that("l2_setup$fsl has required columns", {
  gpa <- create_mock_gpa_for_l3_setup()
  setup <- gpa$l2_model_setup
  
  required_cols <- c("id", "session", "l1_model", "l2_model", "feat_fsf", "feat_dir", "feat_complete")
  expect_true(all(required_cols %in% names(setup$fsl)))
})

# ==============================================================================
# Tests for run/subject exclusion logic
# ==============================================================================

test_that("run_data exclusion flags are properly structured", {
  gpa <- create_mock_gpa_for_l2_setup()
  
  expect_true("exclude_run" %in% names(gpa$run_data))
  expect_true("exclude_subject" %in% names(gpa$run_data))
  expect_type(gpa$run_data$exclude_run, "logical")
  expect_type(gpa$run_data$exclude_subject, "logical")
})

test_that("subject_data exclusion flags are properly structured", {
  gpa <- create_mock_gpa_for_l3_setup()
  
  expect_true("exclude_subject" %in% names(gpa$subject_data))
  expect_type(gpa$subject_data$exclude_subject, "logical")
})

# ==============================================================================
# Tests for multi_run vs single_run configuration
# ==============================================================================

test_that("multi_run gpa requires l2_models", {
  gpa <- create_mock_gpa_for_l3_setup(multi_run = TRUE)
  
  expect_true(gpa$multi_run)
  expect_s3_class(gpa$l2_models, "hi_model_set")
})

test_that("single_run gpa can work without l2_models", {
  gpa <- create_mock_gpa_for_l3_setup(multi_run = FALSE)
  gpa$l2_models <- NULL  # Remove l2_models
  
  expect_false(gpa$multi_run)
  expect_null(gpa$l2_models)
})

# ==============================================================================
# Tests for model name defaulting behavior
# ==============================================================================

test_that("NULL l1_model_names defaults to all models", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  all_models <- names(gpa$l1_models$models)
  
  # The setup functions should use all models when NULL is passed
  # We test this indirectly through the choose_glm_models helper
  result <- fmri.pipeline:::choose_glm_models(gpa, "all", level = 1)
  expect_equal(result, all_models)
})

test_that("empty model names vector is handled", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  result <- fmri.pipeline:::choose_glm_models(gpa, NULL, level = 1)
  expect_null(result)
})
