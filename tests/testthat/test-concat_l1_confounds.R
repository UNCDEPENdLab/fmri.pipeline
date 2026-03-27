test_that("concat_l1_confounds unions spikes and pads zeros", {
  skip_if_not_installed("data.table")

  tmp_dir <- tempfile("confound_concat_")
  dir.create(tmp_dir, recursive = TRUE)
  sqlite_db <- file.path(tmp_dir, "test.sqlite")

  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = FALSE, output_directory = tmp_dir)
  gpa$output_locations$sqlite_db <- sqlite_db
  gpa$glm_software <- "spm"
  gpa$parallel$compute_environment <- list(global = character(0), fsl = character(0), afni = character(0), spm = character(0), r = character(0))

  run1 <- data.frame(
    csf = c(1, 2, 3),
    wm = c(4, 5, 6),
    expr1_spike_1 = c(0, 1, 0),
    stringsAsFactors = FALSE
  )
  run2 <- data.frame(
    csf = c(7, 8, 9),
    wm = c(10, 11, 12),
    expr1_spike_1 = c(0, 0, 1),
    expr1_spike_2 = c(1, 0, 0),
    stringsAsFactors = FALSE
  )

  f1 <- file.path(tmp_dir, "run1_l1_confounds.txt")
  f2 <- file.path(tmp_dir, "run2_l1_confounds.txt")
  write.table(run1, f1, row.names = FALSE, col.names = FALSE)
  write.table(run2, f2, row.names = FALSE, col.names = FALSE)

  cols1 <- data.frame(
    col_index = 1:3,
    col_name = c("csf", "wm", "expr1_spike_1"),
    is_spike = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  cols2 <- data.frame(
    col_index = 1:4,
    col_name = c("csf", "wm", "expr1_spike_1", "expr1_spike_2"),
    is_spike = c(FALSE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  fmri.pipeline:::insert_df_sqlite(gpa, id = "sub1", session = 1L, run_number = 1L, data = cols1, table = "l1_confound_columns", immediate = TRUE)
  fmri.pipeline:::insert_df_sqlite(gpa, id = "sub1", session = 1L, run_number = 2L, data = cols2, table = "l1_confound_columns", immediate = TRUE)

  concat_file <- fmri.pipeline:::concat_l1_confounds(
    gpa = gpa, id = "sub1", session = 1L,
    run_numbers = c(1L, 2L),
    confound_files = c(f1, f2),
    output_dir = tmp_dir
  )

  expect_true(file.exists(concat_file))
  concat_df <- data.table::fread(concat_file, header = FALSE, data.table = FALSE)
  expect_equal(nrow(concat_df), 6)
  expect_equal(ncol(concat_df), 4)

  # run1 rows should have zeros in the final spike column
  expect_true(all(concat_df[1:3, 4] == 0))
})
