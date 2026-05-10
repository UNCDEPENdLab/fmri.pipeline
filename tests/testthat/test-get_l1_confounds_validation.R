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
    get_l1_confounds(run_df = run_df, gpa = gpa),
    regexp = "requested confound columns are missing"
  )

  expect_match(conditionMessage(err), "confound_input_colnames")
  expect_false(grepl("replacement has 1 row, data has 0", conditionMessage(err), fixed = TRUE))
})
