library(testthat)

make_setup_l1_test_gpa <- function(n_subjects = 1, n_runs = 2) {
  gpa <- create_mock_gpa(
    n_subjects = n_subjects,
    n_runs = n_runs,
    include_l1_models = TRUE
  )

  gpa$run_data$run_nifti_present <- TRUE
  gpa$run_data$tr <- 1.0
  gpa$run_data$nvoxels <- 100000L
  gpa$output_locations$setup_l1_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l1_log_txt <- tempfile(fileext = ".txt")
  gpa$lgr_threshold <- "warn"
  gpa$log_json <- FALSE
  gpa$log_txt <- FALSE
  gpa$parallel$l1_setup_cores <- 1L

  gpa
}

test_that("prepare_subject_run_context drops runs without matching event rows", {
  gpa <- make_setup_l1_test_gpa(n_subjects = 1, n_runs = 2)
  tmp_dir <- tempfile("prepare_subject_ctx_")
  dir.create(tmp_dir, recursive = TRUE)

  gpa$output_directory <- tmp_dir
  gpa$subject_data$mr_dir <- file.path(tmp_dir, gpa$subject_data$id)
  gpa$run_data$mr_dir <- file.path(tmp_dir, gpa$run_data$id, paste0("run", gpa$run_data$run_number))
  gpa$run_data$run_nifti <- file.path(gpa$run_data$mr_dir, "bold.nii.gz")
  for (d in unique(gpa$run_data$mr_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (nii in gpa$run_data$run_nifti) file.create(nii)

  gpa$glm_software <- character(0)
  gpa$level_backends <- list(l1 = character(0), l2 = character(0), l3 = character(0))
  gpa$l1_models$events <- list(
    stimulus = list(
      name = "stimulus",
      data = subset(
        gpa$trial_data,
        run_number == 1,
        select = c("id", "session", "run_number", "trial", "onset", "duration", "event")
      )
    )
  )

  ctx <- fmri.pipeline:::init_l1_setup_context(gpa, "model1")
  subject_ctx <- fmri.pipeline:::prepare_subject_run_context(
    subj_df = gpa$subject_data[1, , drop = FALSE],
    gpa = ctx$gpa,
    ctx = ctx
  )

  expect_equal(subject_ctx$run_df$run_number, 1)
  expect_equal(subject_ctx$mr_df$run_number, 1)
  expect_equal(sort(unique(subject_ctx$m_events$run_number)), 1)
})


test_that("load_cached_l1_bdm normalizes legacy run_4d_files caches", {
  tmp_file <- tempfile(fileext = ".RData")
  d_obj <- list(run_4d_files = c("run1.nii.gz", "run2.nii.gz"))
  class(d_obj) <- c("list", "bdm")
  save(d_obj, file = tmp_file)

  lg <- lgr::get_logger("tests/setup_l1_cache")
  lg$set_threshold("error")
  cached <- fmri.pipeline:::load_cached_l1_bdm(tmp_file, lg)

  expect_false(cached$run_bdm)
  expect_equal(cached$d_obj$run_nifti, c("run1.nii.gz", "run2.nii.gz"))
})


test_that("resolve_spm_confound_inputs concatenates valid confounds when configured", {
  tmp_dir <- tempfile("resolve_spm_confounds_")
  dir.create(tmp_dir, recursive = TRUE)
  conf1 <- file.path(tmp_dir, "conf1.txt")
  conf2 <- file.path(tmp_dir, "conf2.txt")
  writeLines("1 2 3", conf1)
  writeLines("4 5 6", conf2)

  gpa <- make_setup_l1_test_gpa(n_subjects = 1, n_runs = 2)
  gpa$glm_settings <- list(spm = list(concatenate_runs = TRUE))

  lg <- lgr::get_logger("tests/setup_l1_confounds")
  lg$set_threshold("error")
  concat_path <- file.path(tmp_dir, "concat.txt")

  local_mocked_bindings(
    concat_l1_confounds = function(...) concat_path,
    .package = "fmri.pipeline"
  )

  out <- fmri.pipeline:::resolve_spm_confound_inputs(
    gpa = gpa,
    id = "sub1",
    session = 1L,
    run_numbers = c(1L, 2L),
    confound_files = c(conf1, conf2),
    output_dir = tmp_dir,
    lg = lg
  )

  expect_equal(out, concat_path)
})


test_that("resolve_spm_confound_inputs returns NULL when concatenation is requested with no confounds", {
  gpa <- make_setup_l1_test_gpa(n_subjects = 1, n_runs = 2)
  gpa$glm_settings <- list(spm = list(concatenate_runs = TRUE))

  lg <- lgr::get_logger("tests/setup_l1_confounds_none")
  lg$set_threshold("error")
  out <- fmri.pipeline:::resolve_spm_confound_inputs(
    gpa = gpa,
    id = "sub1",
    session = 1L,
    run_numbers = c(1L, 2L),
    confound_files = c("", ""),
    output_dir = tempdir(),
    lg = lg
  )

  expect_null(out)
})


test_that("prepare_pooled_spm_context remaps session-specific run numbers into pooled order", {
  tmp_dir <- tempfile("prepare_pooled_spm_")
  dir.create(tmp_dir, recursive = TRUE)

  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 1, include_l1_models = TRUE, output_directory = tmp_dir)
  gpa$subject_data <- data.frame(
    id = c("sub1", "sub1"),
    session = c(1L, 2L),
    mr_dir = file.path(tmp_dir, "sub1"),
    stringsAsFactors = FALSE
  )
  gpa$run_data <- data.frame(
    id = c("sub1", "sub1"),
    session = c(1L, 2L),
    run_number = c(1L, 1L),
    run_volumes = c(100L, 120L),
    drop_volumes = c(2L, 2L),
    mr_dir = c(file.path(tmp_dir, "sub1", "ses-1"), file.path(tmp_dir, "sub1", "ses-2")),
    run_nifti = c("bold.nii.gz", "bold.nii.gz"),
    exclude_run = c(FALSE, FALSE),
    has_l1_setup = c(FALSE, FALSE),
    has_l1_complete = c(FALSE, FALSE),
    run_nifti_present = c(TRUE, TRUE),
    tr = c(1.0, 1.0),
    nvoxels = c(100000L, 100000L),
    stringsAsFactors = FALSE
  )
  for (d in unique(gpa$run_data$mr_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(d, "bold.nii.gz"))
  }

  gpa$glm_software <- "spm"
  gpa$level_backends <- list(l1 = "spm", l2 = character(0), l3 = character(0))
  gpa$glm_settings <- list(spm = list(concatenate_runs = FALSE))
  gpa$additional <- list(bdm_args = list(plot = FALSE))
  gpa$use_preconvolve <- FALSE
  gpa$output_locations$spm_l1_directory <- file.path(tmp_dir, "spm_l1", "{id}", "{l1_model}")

  signal_values <- expand.grid(
    id = "sub1",
    session = c(1L, 2L),
    run_number = 1L,
    trial = 1:2,
    stringsAsFactors = FALSE
  )
  signal_values$value <- seq_len(nrow(signal_values))
  gpa$l1_models$signals$parametric_signal$value <- signal_values
  gpa$l1_models$models$model1$signals <- "parametric_signal"
  gpa$l1_models$events <- list(
    stimulus = list(
      name = "stimulus",
      data = data.frame(
        id = c("sub1", "sub1"),
        session = c(1L, 2L),
        run_number = c(1L, 1L),
        trial = c(1L, 1L),
        onset = c(0.5, 0.75),
        duration = c(1, 1),
        event = c("stimulus", "stimulus"),
        stringsAsFactors = FALSE
      )
    )
  )

  lg <- lgr::get_logger("tests/setup_l1_pooled")
  lg$set_threshold("error")
  pooled_ctx <- fmri.pipeline:::prepare_pooled_spm_context(
    subject_ctx = list(subj_id = "sub1"),
    gpa = gpa,
    model_name = "model1",
    lg = lg
  )

  expect_equal(pooled_ctx$pooled_run_df$run_number, c(1L, 2L))
  expect_equal(pooled_ctx$pooled_run_df$source_run_number, c(1L, 1L))
  expect_equal(pooled_ctx$m_events_spm$run_number, c(1L, 2L))
  expect_equal(pooled_ctx$mr_df_spm$session, c(1L, 2L))
})

test_that("drop_runs_without_events drops pooled SPM runs with missing events", {
  lg <- lgr::get_logger("tests/drop_runs_pooled")
  lg$set_threshold("error")

  run_df <- data.frame(
    run_number = 1:3,
    run_nifti = paste0("bold", 1:3, ".nii.gz"),
    stringsAsFactors = FALSE
  )
  events <- data.frame(
    run_number = c(1L, 1L, 3L),
    trial = c(1L, 2L, 1L),
    onset = c(0.0, 2.0, 0.0),
    duration = c(1, 1, 1),
    stringsAsFactors = FALSE
  )

  # Pooled path: run 2 has imaging but no events
  result <- fmri.pipeline:::drop_runs_without_events(
    run_df = run_df,
    events = events,
    subj_id = "sub1",
    model_name = "model1",
    lg = lg,
    pooled = TRUE
  )

  expect_equal(result$run_df$run_number, c(1L, 3L))
  expect_equal(sort(unique(result$events$run_number)), c(1L, 3L))
  expect_equal(nrow(result$events), 3L)
})

test_that("drop_runs_without_events is a no-op when all runs have events", {
  lg <- lgr::get_logger("tests/drop_runs_noop")
  lg$set_threshold("error")

  run_df <- data.frame(
    run_number = 1:2,
    run_nifti = paste0("bold", 1:2, ".nii.gz"),
    stringsAsFactors = FALSE
  )
  events <- data.frame(
    run_number = c(1L, 2L),
    trial = c(1L, 1L),
    onset = c(0.0, 0.0),
    duration = c(1, 1),
    stringsAsFactors = FALSE
  )

  result <- fmri.pipeline:::drop_runs_without_events(
    run_df = run_df,
    events = events,
    subj_id = "sub1",
    subj_session = 1L,
    lg = lg,
    pooled = FALSE
  )

  expect_equal(result$run_df$run_number, 1:2)
  expect_equal(nrow(result$events), 2L)
})
