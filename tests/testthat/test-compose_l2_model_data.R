library(testthat)

test_that("compose_l2_model_data skips identical session-level duplicates already present in run_data", {
  gpa <- create_mock_gpa()
  gpa$run_data$tr <- 1
  gpa$run_data$session_label <- paste0("ses", gpa$run_data$session)

  gpa$session_data <- unique(gpa$run_data[, c("id", "session", "tr", "session_label")])
  gpa$session_data$session_covariate <- seq_len(nrow(gpa$session_data))

  result <- fmri.pipeline:::compose_l2_model_data(gpa)

  expect_false("session_tr" %in% names(result))
  expect_false("session_session_label" %in% names(result))
  expect_true("session_covariate" %in% names(result))
  expect_identical(attr(result, "l2_var_origin")[["tr"]], "run")
  expect_identical(attr(result, "l2_var_origin")[["session_label"]], "run")
  expect_identical(attr(result, "l2_var_origin")[["session_covariate"]], "session")
  expect_null(attr(result, "session_rename_map"))
})

test_that("compose_l2_model_data renames only conflicting session-level columns with different values", {
  gpa <- create_mock_gpa()
  gpa$run_data$tr <- 1
  gpa$run_data$session_label <- paste0("ses", gpa$run_data$session)

  gpa$session_data <- unique(gpa$run_data[, c("id", "session", "tr", "session_label")])
  gpa$session_data$tr <- gpa$session_data$tr + 0.5

  result <- fmri.pipeline:::compose_l2_model_data(gpa)

  expect_true("session_tr" %in% names(result))
  expect_false("session_session_label" %in% names(result))
  expect_identical(unname(attr(result, "session_rename_map")), "session_tr")
  expect_identical(names(attr(result, "session_rename_map")), "tr")
  expect_equal(unique(result$tr), 1)
  expect_equal(unique(result$session_tr), 1.5)
})
