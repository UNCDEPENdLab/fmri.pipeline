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

  out <- fmri.pipeline:::get_l1_confounds(run_df = run_df, gpa = gpa)
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

test_that("confound defaults do not assume motion, confounds, or spike regressors", {
  defaults <- eval(formals(fmri.pipeline::setup_glm_pipeline)$confound_settings)

  expect_null(defaults$motion_params_file)
  expect_null(defaults$motion_params_colnames)
  expect_null(defaults$confound_input_file)
  expect_null(defaults$spike_volumes)
})

test_that("spike_volumes is disabled when no motion files are present", {
  lg <- lgr::get_logger("test/disable_unavailable_spike_volumes")
  lg$set_threshold("fatal")

  gpa <- list(
    confound_settings = list(spike_volumes = "framewise_displacement > 0.9"),
    run_data = data.frame(
      id = "sub01",
      session = 1L,
      run_number = 1L,
      motion_params_present = FALSE,
      stringsAsFactors = FALSE
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- fmri.pipeline:::disable_unavailable_spike_volumes(gpa, lg)

  expect_null(out$confound_settings$spike_volumes)
})

test_that("read_confounds_motion_parameters returns an empty confound frame when no inputs are present", {
  tdir <- tempfile("missing_input_fallback_")
  dir.create(tdir, recursive = TRUE)

  lg <- lgr::get_logger("test/read_confounds_motion_parameters")
  lg$set_threshold("warn")

  run_df <- data.frame(
    id = "sub01",
    session = 1L,
    run_number = 1L,
    run_nifti = file.path(tdir, "run.nii.gz"),
    run_volumes = 3L,
    tr = 2,
    motion_params_file = NA_character_,
    motion_params_present = FALSE,
    confound_input_file = NA_character_,
    confound_input_file_present = FALSE,
    stringsAsFactors = FALSE
  )

  gpa <- list(
    output_locations = list(sqlite_db = file.path(tdir, "cache.sqlite")),
    confound_settings = list(
      motion_params_colnames = NULL,
      na_strings = "NA"
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- fmri.pipeline:::read_confounds_motion_parameters(
    gpa = gpa,
    id = "sub01",
    session = 1L,
    run_number = 1L,
    run_df = run_df,
    expected_l1_confound_file = NA_character_,
    lg = lg
  )

  expect_true(is.data.frame(out$all_confounds))
  expect_null(out$motion_df)
  expect_equal(names(out$all_confounds), c("volume", "time"))
  expect_equal(out$all_confounds$volume, 1:3)
  expect_equal(out$all_confounds$time, c(0, 2, 4))
})

test_that("read_confounds_motion_parameters uses gpa TR when run_df has no tr column", {
  tdir <- tempfile("missing_input_gpa_tr_")
  dir.create(tdir, recursive = TRUE)

  lg <- lgr::get_logger("test/read_confounds_motion_parameters_gpa_tr")
  lg$set_threshold("warn")

  run_df <- data.frame(
    id = "sub01",
    session = 1L,
    run_number = 1L,
    run_nifti = file.path(tdir, "run.nii.gz"),
    run_volumes = 4L,
    motion_params_file = NA_character_,
    motion_params_present = FALSE,
    confound_input_file = NA_character_,
    confound_input_file_present = FALSE,
    stringsAsFactors = FALSE
  )

  gpa <- list(
    tr = 1.5,
    output_locations = list(sqlite_db = file.path(tdir, "cache.sqlite")),
    confound_settings = list(
      motion_params_colnames = NULL,
      na_strings = "NA"
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- fmri.pipeline:::read_confounds_motion_parameters(
    gpa = gpa,
    id = "sub01",
    session = 1L,
    run_number = 1L,
    run_df = run_df,
    expected_l1_confound_file = NA_character_,
    lg = lg
  )

  expect_equal(out$all_confounds$time, c(0, 1.5, 3, 4.5))
})

test_that("get_l1_confounds handles spike request without motion by returning no confound file", {
  tdir <- tempfile("spike_without_motion_")
  dir.create(tdir, recursive = TRUE)

  sqlite_db <- file.path(tdir, "cache.sqlite")
  run_nifti <- file.path(tdir, "run.nii.gz")
  file.create(run_nifti)

  run_df <- data.frame(
    id = "sub01",
    session = 1L,
    run_number = 1L,
    run_nifti = run_nifti,
    run_volumes = 5L,
    drop_volumes = 0L,
    tr = 1,
    motion_params_file = NA_character_,
    motion_params_present = FALSE,
    confound_input_file = NA_character_,
    confound_input_file_present = FALSE,
    stringsAsFactors = FALSE
  )

  gpa <- list(
    analysis_name = "spike_without_motion_test",
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
      spike_volumes = "framewise_displacement > 0.9",
      na_strings = "NA"
    ),
    drop_volumes = 0L,
    lgr_threshold = "fatal"
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- fmri.pipeline:::get_l1_confounds(run_df = run_df, gpa = gpa)

  expect_true(is.na(out$l1_confound_file[1L]))
  expect_false(out$exclude_run[1L])
  expect_equal(out$first_volume[1L], 1L)
  expect_equal(out$last_volume[1L], 5L)

  cached <- fmri.pipeline:::read_df_sqlite(
    gpa = gpa, id = "sub01", session = 1L, run_number = 1L,
    table = "l1_run_calculations"
  )
  expect_equal(nrow(cached), 1L)
  expect_true(is.na(cached$l1_confound_file[1L]))
})

test_that("test_generate_l1_confounds distinguishes no-confound and requested-confound cases", {
  tdir <- tempfile("generate_l1_confounds_")
  dir.create(tdir, recursive = TRUE)

  lg <- lgr::get_logger("test/test_generate_l1_confounds")
  lg$set_threshold("warn")

  gpa <- list(
    confound_settings = list(
      l1_confound_regressors = NULL,
      spike_volumes = NULL
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- fmri.pipeline:::test_generate_l1_confounds(gpa, tdir, 2L, lg)
  expect_false(out$generate_l1_confounds)
  expect_true(is.na(out$expected_l1_confound_file))

  gpa$confound_settings$spike_volumes <- "framewise_displacement > 0.9"
  out <- fmri.pipeline:::test_generate_l1_confounds(gpa, tdir, 2L, lg)
  expect_true(out$generate_l1_confounds)
  expect_equal(out$expected_l1_confound_file, file.path(tdir, "run2_l1_confounds.txt"))

  file.create(out$expected_l1_confound_file)
  out <- fmri.pipeline:::test_generate_l1_confounds(gpa, tdir, 2L, lg)
  expect_false(out$generate_l1_confounds)
  expect_equal(out$expected_l1_confound_file, file.path(tdir, "run2_l1_confounds.txt"))
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
