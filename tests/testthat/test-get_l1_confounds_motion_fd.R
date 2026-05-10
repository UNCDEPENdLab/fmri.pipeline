test_that("get_l1_confounds exposes motion-derived framewise_displacement to exclude_run", {
  skip_if_not_installed("data.table")

  tmp_dir <- tempfile("get_l1_confounds_motion_fd_")
  dir.create(tmp_dir, recursive = TRUE)

  confound_file <- file.path(tmp_dir, "confounds.tsv")
  motion_file <- file.path(tmp_dir, "motion.par")
  sqlite_db <- file.path(tmp_dir, "test.sqlite")

  write.table(
    data.frame(
      csf = c(10, 11, 12, 13),
      wm = c(20, 21, 22, 23)
    ),
    file = confound_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  motion_mat <- matrix(
    c(
      0, 0, 0, 0.0, 0, 0,
      0, 0, 0, 0.1, 0, 0,
      0, 0, 0, 0.2, 0, 0,
      0, 0, 0, 0.3, 0, 0
    ),
    ncol = 6,
    byrow = TRUE
  )
  write.table(
    motion_mat,
    file = motion_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  gpa <- list(
    output_directory = tmp_dir,
    output_locations = list(
      sqlite_db = sqlite_db,
      feat_sub_directory = "{gpa$output_directory}/feat_l1/sub-{id}"
    ),
    run_data = data.frame(
      id = "sub1",
      session = 1L,
      run_number = 1L,
      run_nifti = file.path(tmp_dir, "bold.nii.gz"),
      run_volumes = 4L,
      drop_volumes = 0L,
      tr = 1,
      confound_input_file = confound_file,
      confound_input_file_present = TRUE,
      motion_params_file = motion_file,
      motion_params_present = TRUE,
      stringsAsFactors = FALSE
    ),
    drop_volumes = 0L,
    confound_settings = list(
      motion_params_colnames = c("rx", "ry", "rz", "tx", "ty", "tz"),
      l1_confound_regressors = c("csf", "wm"),
      exclude_run = "max(framewise_displacement) > 0.09",
      run_exclusion_columns = "framewise_displacement",
      truncate_run = NULL,
      spike_volumes = NULL,
      na_strings = "NA"
    ),
    lgr_threshold = "fatal"
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- fmri.pipeline:::get_l1_confounds(
    id = "sub1",
    session = 1L,
    run_number = 1L,
    gpa = gpa
  )

  expect_true(out$exclude_run)
  expect_true("framewise_displacement" %in% names(out$exclude_data))
  expect_equal(
    out$exclude_data$framewise_displacement,
    c(0, 0.1, 0.1, 0.1),
    tolerance = 1e-8
  )

  motion_df <- fmri.pipeline:::read_df_sqlite(
    gpa = gpa,
    id = "sub1",
    session = 1L,
    run_number = 1L,
    table = "l1_motion_parameters"
  )
  expect_true("framewise_displacement" %in% names(motion_df))
  expect_equal(
    motion_df$framewise_displacement,
    c(0, 0.1, 0.1, 0.1),
    tolerance = 1e-8
  )
})

test_that("get_l1_confounds handles a single retained confound column", {
  skip_if_not_installed("data.table")

  tmp_dir <- tempfile("get_l1_confounds_single_col_")
  dir.create(tmp_dir, recursive = TRUE)

  confound_file <- file.path(tmp_dir, "confounds.tsv")
  sqlite_db <- file.path(tmp_dir, "test.sqlite")

  write.table(
    data.frame(csf = c(10, 11, 12, 13)),
    file = confound_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  gpa <- list(
    output_directory = tmp_dir,
    output_locations = list(
      sqlite_db = sqlite_db,
      feat_sub_directory = "{gpa$output_directory}/feat_l1/sub-{id}"
    ),
    run_data = data.frame(
      id = "sub1",
      session = 1L,
      run_number = 1L,
      run_nifti = file.path(tmp_dir, "bold.nii.gz"),
      run_volumes = 4L,
      drop_volumes = 0L,
      tr = 1,
      confound_input_file = confound_file,
      confound_input_file_present = TRUE,
      motion_params_file = NA_character_,
      motion_params_present = FALSE,
      stringsAsFactors = FALSE
    ),
    drop_volumes = 0L,
    confound_settings = list(
      motion_params_colnames = c("rx", "ry", "rz", "tx", "ty", "tz"),
      l1_confound_regressors = "csf",
      exclude_run = NULL,
      run_exclusion_columns = NULL,
      truncate_run = NULL,
      spike_volumes = NULL,
      na_strings = "NA"
    ),
    lgr_threshold = "fatal"
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  out <- fmri.pipeline:::get_l1_confounds(
    id = "sub1",
    session = 1L,
    run_number = 1L,
    gpa = gpa
  )

  expect_false(out$exclude_run)
  expect_equal(names(out$l1_confounds_df), "csf")
  expect_equal(ncol(out$l1_confounds_df), 1L)
  expect_equal(out$l1_confounds_df$csf, c(-1.5, -0.5, 0.5, 1.5))
})
