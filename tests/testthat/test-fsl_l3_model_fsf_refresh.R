test_that("fsl_l3_model rewrites existing FSF when model is queued to run", {
  outdir <- tempfile("l3_fsf_refresh_")
  dir.create(outdir, recursive = TRUE)

  l3_df <- data.frame(
    id = "sub1",
    session = 1L,
    l1_model = "l1m",
    l2_model = "l2m",
    l3_model = "l3m",
    l1_cope_name = "EV_choice",
    l2_cope_name = "EV_Intercept",
    cope_file = "dummy.nii.gz",
    stringsAsFactors = FALSE
  )

  gpa <- list(
    lgr_threshold = "warn",
    multi_run = TRUE,
    l3_models = list(models = list(l3m = list())),
    output_locations = list(feat_l3_fsf = "FEAT_l1c-{l1_contrast}.fsf"),
    additional = list(feat_l3_args = NULL),
    glm_settings = list(fsl = list(force_l3_creation = FALSE))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  l3_fsf <- file.path(outdir, "FEAT_l1c-EV_choice.fsf")
  writeLines(c("# SENTINEL_OLD_FILE", "set fmri(robust_yn) 1"), con = l3_fsf)

  mock_feat_status <- function(feat_dir, feat_fsf, feat_complete) {
    data.frame(
      feat_fsf = feat_fsf,
      feat_fsf_modified_date = Sys.time(),
      feat_fsf_exists = file.exists(feat_fsf),
      feat_dir = feat_dir,
      feat_dir_exists = dir.exists(feat_dir),
      feat_execution_start = as.POSIXct(NA),
      feat_execution_end = as.POSIXct(NA),
      feat_execution_min = NA_real_,
      feat_complete = feat_complete,
      feat_failed = !feat_complete,
      stringsAsFactors = FALSE
    )
  }

  testthat::local_mocked_bindings(
    respecify_l3_model = function(mobj, new_data) {
      list(
        metadata = new_data[, c("id", "session"), drop = FALSE],
        model_matrix = matrix(1, nrow = nrow(new_data), ncol = 1),
        contrasts = matrix(1, nrow = 1, ncol = 1),
        fsl_outlier_deweighting = FALSE
      )
    },
    fsl_generate_fsf_ev_syntax = function(inputs, dmat) character(0),
    fsl_generate_fsf_contrast_syntax = function(contrasts) character(0),
    get_output_directory = function(...) outdir,
    get_feat_status = function(feat_dir, feat_fsf, lg = NULL) {
      mock_feat_status(feat_dir, feat_fsf, feat_complete = FALSE)
    },
    .package = "fmri.pipeline"
  )

  result <- fmri.pipeline:::fsl_l3_model(l3_df = l3_df, gpa = gpa)
  expect_true(isTRUE(result$to_run[1L]))

  rewritten <- readLines(l3_fsf)
  expect_true(any(grepl("^set fmri\\(robust_yn\\) 0$", rewritten)))
  expect_false(any(grepl("SENTINEL_OLD_FILE", rewritten, fixed = TRUE)))
})


test_that("fsl_l3_model keeps existing FSF when model is already complete", {
  outdir <- tempfile("l3_fsf_skip_")
  dir.create(outdir, recursive = TRUE)

  l3_df <- data.frame(
    id = "sub1",
    session = 1L,
    l1_model = "l1m",
    l2_model = "l2m",
    l3_model = "l3m",
    l1_cope_name = "EV_choice",
    l2_cope_name = "EV_Intercept",
    cope_file = "dummy.nii.gz",
    stringsAsFactors = FALSE
  )

  gpa <- list(
    lgr_threshold = "warn",
    multi_run = TRUE,
    l3_models = list(models = list(l3m = list())),
    output_locations = list(feat_l3_fsf = "FEAT_l1c-{l1_contrast}.fsf"),
    additional = list(feat_l3_args = NULL),
    glm_settings = list(fsl = list(force_l3_creation = FALSE))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  l3_fsf <- file.path(outdir, "FEAT_l1c-EV_choice.fsf")
  writeLines(c("# SENTINEL_OLD_FILE", "set fmri(robust_yn) 1"), con = l3_fsf)

  mock_feat_status <- function(feat_dir, feat_fsf, feat_complete) {
    data.frame(
      feat_fsf = feat_fsf,
      feat_fsf_modified_date = Sys.time(),
      feat_fsf_exists = file.exists(feat_fsf),
      feat_dir = feat_dir,
      feat_dir_exists = dir.exists(feat_dir),
      feat_execution_start = as.POSIXct(NA),
      feat_execution_end = as.POSIXct(NA),
      feat_execution_min = NA_real_,
      feat_complete = feat_complete,
      feat_failed = !feat_complete,
      stringsAsFactors = FALSE
    )
  }

  testthat::local_mocked_bindings(
    respecify_l3_model = function(mobj, new_data) {
      list(
        metadata = new_data[, c("id", "session"), drop = FALSE],
        model_matrix = matrix(1, nrow = nrow(new_data), ncol = 1),
        contrasts = matrix(1, nrow = 1, ncol = 1),
        fsl_outlier_deweighting = FALSE
      )
    },
    fsl_generate_fsf_ev_syntax = function(inputs, dmat) character(0),
    fsl_generate_fsf_contrast_syntax = function(contrasts) character(0),
    get_output_directory = function(...) outdir,
    get_feat_status = function(feat_dir, feat_fsf, lg = NULL) {
      mock_feat_status(feat_dir, feat_fsf, feat_complete = TRUE)
    },
    .package = "fmri.pipeline"
  )

  result <- fmri.pipeline:::fsl_l3_model(l3_df = l3_df, gpa = gpa)
  expect_false(isTRUE(result$to_run[1L]))

  existing <- readLines(l3_fsf)
  expect_true(any(grepl("SENTINEL_OLD_FILE", existing, fixed = TRUE)))
  expect_true(any(grepl("^set fmri\\(robust_yn\\) 1$", existing)))
})
