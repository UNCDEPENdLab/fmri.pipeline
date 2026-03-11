library(testthat)

test_that("get_l1_confounds initializes defaults when cache is missing and no confounds are requested", {
  tdir <- tempfile("confound_cache_fallback_")
  dir.create(tdir, recursive = TRUE)

  sqlite_db <- file.path(tdir, "cache.sqlite")
  run_nifti <- file.path(tdir, "run.nii.gz")
  file.create(run_nifti)

  run_df <- data.frame(
    id = "sub01",
    session = 1L,
    run_number = 1L,
    run_nifti = run_nifti,
    run_volumes = 120L,
    drop_volumes = 0L,
    stringsAsFactors = FALSE
  )

  gpa <- list(
    analysis_name = "cache_fallback_test",
    output_directory = tdir,
    output_locations = list(
      sqlite_db = sqlite_db,
      feat_sub_directory = file.path(tdir, "feat_l1", "sub-{id}", "ses-{session}")
    ),
    run_data = run_df,
    confound_settings = list(
      motion_params_file = NULL,
      motion_params_colnames = NULL,
      confound_input_file = NULL,
      confound_input_colnames = NULL,
      l1_confound_regressors = NULL,
      exclude_run = NULL,
      truncate_run = NULL,
      exclude_subject = NULL,
      spike_volumes = NULL
    ),
    drop_volumes = 0L,
    lgr_threshold = "warn"
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- get_l1_confounds(run_df = run_df, gpa = gpa)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_false(out$exclude_run[1L])
  expect_equal(out$orig_volumes[1L], 120L)
  expect_equal(out$first_volume[1L], 1L)
  expect_equal(out$last_volume[1L], 120L)
  expect_equal(out$truncate_volumes[1L], 0L)
  expect_equal(out$run_volumes[1L], 120L)
  expect_true(is.na(out$l1_confound_file[1L]))

  cached <- fmri.pipeline:::read_df_sqlite(
    gpa = gpa, id = "sub01", session = 1L, run_number = 1L,
    table = "l1_run_calculations"
  )
  expect_true(is.data.frame(cached))
  expect_equal(nrow(cached), 1L)
  expect_equal(cached$first_volume[1L], 1L)
  expect_equal(cached$last_volume[1L], 120L)
})

test_that("test_generate_run_exclusion always returns exclude_run", {
  lg <- lgr::get_logger("test/get_l1_confounds")
  lg$set_threshold("warn")
  gpa <- list(
    confound_settings = list(exclude_run = NULL),
    output_locations = list(sqlite_db = tempfile(fileext = ".sqlite"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- fmri.pipeline:::test_generate_run_exclusion(
    gpa = gpa,
    id = "sub01",
    session = 1L,
    run_number = 1L,
    l1_cached_df = NULL,
    lg = lg
  )

  expect_true("exclude_run" %in% names(out))
  expect_identical(out$exclude_run, FALSE)
})
