library(testthat)

test_that("get_3dlmer_status returns a refresh-compatible data frame", {
  out_dir <- tempfile("afni_3dlmer_status_")
  dir.create(out_dir, recursive = TRUE)
  prefix <- file.path(out_dir, "lmer_result")
  head_file <- paste0(prefix, "+tlrc.HEAD")
  file.create(head_file)

  status <- fmri.pipeline:::get_3dlmer_status(prefix)

  expect_true(is.data.frame(status))
  expect_identical(status$output_file, prefix)
  expect_true(status$output_file_exists)
  expect_identical(status$output_file_resolved, head_file)
  expect_true(status$feat_complete)
  expect_false(status$feat_failed)
})

test_that("refresh_glm_status updates AFNI L3 tables via get_3dlmer_status", {
  out_dir <- tempfile("afni_3dlmer_refresh_")
  dir.create(out_dir, recursive = TRUE)
  prefix <- file.path(out_dir, "lmer_result")
  file.create(paste0(prefix, "+tlrc.HEAD"))

  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- "afni"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  gpa$l3_model_setup <- structure(
    list(
      afni = data.frame(
        output_file = prefix,
        feat_complete = FALSE,
        stringsAsFactors = FALSE
      )
    ),
    class = c("l3_setup", "list")
  )

  res <- fmri.pipeline:::refresh_glm_status(gpa, level = 3L, glm_software = "afni")

  expect_true(res$l3_model_setup$afni$feat_complete[1])
  expect_true(res$l3_model_setup$afni$output_file_exists[1])
  expect_identical(res$l3_model_setup$afni$output_file_resolved[1], paste0(prefix, "+tlrc.HEAD"))
})
