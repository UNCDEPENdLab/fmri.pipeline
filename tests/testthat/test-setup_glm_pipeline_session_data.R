library(testthat)

test_that("setup_glm_pipeline derives session_data internally from run and subject data", {
  subject_data <- data.frame(
    id = c("sub1", "sub1", "sub2", "sub2"),
    session = c(1L, 2L, 1L, 2L),
    session_label = c("baseline", "followup", "baseline", "followup"),
    group = c("control", "control", "patient", "patient"),
    exclude_subject = FALSE,
    stringsAsFactors = FALSE
  )

  run_data <- expand.grid(
    id = c("sub1", "sub2"),
    session = c(1L, 2L),
    run_number = c(1L, 2L),
    stringsAsFactors = FALSE
  )
  run_data <- run_data[order(run_data$id, run_data$session, run_data$run_number), ]
  run_data$tr <- 1.5
  run_data$session_label <- subject_data$session_label[
    match(paste(run_data$id, run_data$session), paste(subject_data$id, subject_data$session))
  ]
  run_data$scanner_day <- ifelse(run_data$session == 1L, 0L, 14L)
  run_data$run_condition <- paste0("run", run_data$run_number)
  run_data$exclude_run <- FALSE

  trial_data <- do.call(
    rbind,
    lapply(seq_len(nrow(run_data)), function(ii) {
      data.frame(
        id = run_data$id[ii],
        session = run_data$session[ii],
        run_number = run_data$run_number[ii],
        trial = seq_len(2L),
        onset = c(0, 4),
        duration = 1,
        condition = c("face", "house"),
        stringsAsFactors = FALSE
      )
    })
  )

  gpa <- setup_glm_pipeline(
    analysis_name = "session_data_internal_test",
    scheduler = "local",
    output_directory = tempdir(),
    subject_data = subject_data,
    run_data = run_data,
    trial_data = trial_data,
    l1_models = NULL,
    l2_models = NULL,
    l3_models = NULL,
    glm_software = "fsl",
    n_expected_runs = 2L,
    lgr_threshold = "off"
  )

  expect_false("session_data" %in% names(formals(setup_glm_pipeline)))
  expect_s3_class(gpa$session_data, "bg_session_data")
  expect_equal(nrow(gpa$session_data), 4L)
  expect_equal(nrow(unique(gpa$session_data[, c("id", "session")])), 4L)
  expect_true(all(c("session_label", "group", "scanner_day", "tr") %in% names(gpa$session_data)))
  expect_false("run_condition" %in% names(gpa$session_data))
  expect_equal(gpa$session_data$scanner_day, c(0L, 14L, 0L, 14L))
})

test_that("setup_glm_pipeline copies incoming data.tables to data.frames", {
  trial_data <- data.table::data.table(
    id = c("sub1", "sub1"),
    trial = 1:2,
    onset = c(0, 4),
    duration = 1,
    condition = c("face", "house")
  )
  output_directory <- tempfile("setup_glm_dt_input_")

  gpa <- setup_glm_pipeline(
    analysis_name = basename(output_directory),
    scheduler = "local",
    output_directory = output_directory,
    trial_data = trial_data,
    tr = 1.5,
    l1_models = NULL,
    l2_models = NULL,
    l3_models = NULL,
    glm_software = "fsl",
    n_expected_runs = 1L,
    lgr_threshold = "off"
  )

  expect_true(data.table::is.data.table(trial_data))
  expect_false("session" %in% names(trial_data))
  expect_false("run_number" %in% names(trial_data))
  expect_false(data.table::is.data.table(gpa$trial_data))
  expect_false(data.table::is.data.table(gpa$run_data))
  expect_false(data.table::is.data.table(gpa$subject_data))
  expect_equal(gpa$trial_data$session, rep(1L, 2L))
  expect_equal(gpa$trial_data$run_number, rep(1L, 2L))
  expect_s3_class(gpa$trial_data, "bg_trial_data")
  expect_s3_class(gpa$run_data, "bg_run_data")
  expect_s3_class(gpa$subject_data, "bg_subject_data")
})
