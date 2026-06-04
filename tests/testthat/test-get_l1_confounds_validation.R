library(testthat)

test_that("get_l1_confounds fails clearly when requested confound columns are missing", {
  tdir <- tempfile("confound_validation_")
  dir.create(tdir, recursive = TRUE)

  sqlite_db <- file.path(tdir, "cache.sqlite")
  run_nifti <- file.path(tdir, "run.nii.gz")
  confounds_tsv <- file.path(tdir, "confounds.tsv")

  file.create(run_nifti)
  write.table(
    data.frame(V1 = 1:5, V2 = 6:10, V3 = 11:15, V4 = 16:20, V5 = 0L),
    file = confounds_tsv,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  run_df <- data.frame(
    id = "sub01",
    session = 1L,
    run_number = 1L,
    run_nifti = run_nifti,
    run_volumes = 5L,
    drop_volumes = 0L,
    tr = 1.0,
    confound_input_file = confounds_tsv,
    confound_input_file_present = TRUE,
    motion_params_present = FALSE,
    stringsAsFactors = FALSE
  )

  gpa <- list(
    analysis_name = "confound_validation_test",
    output_directory = tdir,
    output_locations = list(
      sqlite_db = sqlite_db,
      feat_sub_directory = file.path(tdir, "feat_l1", "sub-{id}")
    ),
    run_data = run_df,
    confound_settings = list(
      motion_params_file = NULL,
      motion_params_colnames = NULL,
      confound_input_file = NULL,
      confound_input_colnames = NULL,
      l1_confound_regressors = c("csf", "csf_derivative1", "white_matter", "white_matter_derivative1"),
      exclude_run = NULL,
      truncate_run = NULL,
      exclude_subject = NULL,
      spike_volumes = NULL
    ),
    drop_volumes = 0L,
    lgr_threshold = "warn"
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  err <- expect_error(
    fmri.pipeline:::get_l1_confounds(run_df = run_df, gpa = gpa),
    regexp = "requested confound columns are missing"
  )

  expect_match(conditionMessage(err), "confound_input_colnames")
  expect_false(grepl("replacement has 1 row, data has 0", conditionMessage(err), fixed = TRUE))
})

test_that("read_confounds_motion_parameters fails clearly when motion and confound rows differ", {
  tdir <- tempfile("confound_motion_row_mismatch_")
  dir.create(tdir, recursive = TRUE)

  confound_file <- file.path(tdir, "confounds.tsv")
  motion_file <- file.path(tdir, "motion.par")

  write.table(
    data.frame(csf = c(1, 2, 3)),
    file = confound_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    matrix(
      c(
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 2, 0, 0,
        0, 0, 0, 3, 0, 0
      ),
      ncol = 6,
      byrow = TRUE
    ),
    file = motion_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  run_df <- data.frame(
    id = "sub01",
    session = 1L,
    run_number = 1L,
    run_nifti = file.path(tdir, "run.nii.gz"),
    run_volumes = 4L,
    tr = 1,
    confound_input_file = confound_file,
    confound_input_file_present = TRUE,
    motion_params_file = motion_file,
    motion_params_present = TRUE,
    stringsAsFactors = FALSE
  )

  gpa <- list(
    output_locations = list(sqlite_db = file.path(tdir, "cache.sqlite")),
    confound_settings = list(
      motion_params_colnames = c("rx", "ry", "rz", "tx", "ty", "tz"),
      confound_input_colnames = NULL,
      na_strings = "NA"
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  lg <- lgr::get_logger("test/read_confounds_motion_parameters_row_mismatch")
  lg$set_threshold("fatal")

  err <- expect_error(
    fmri.pipeline:::read_confounds_motion_parameters(
      gpa = gpa,
      id = "sub01",
      session = 1L,
      run_number = 1L,
      run_df = run_df,
      expected_l1_confound_file = file.path(tdir, "run1_l1_confounds.txt"),
      lg = lg
    ),
    regexp = "Number of rows in motion_df is 4 and in confound_df is 3"
  )

  expect_false(grepl("empty_set", conditionMessage(err), fixed = TRUE))
})
