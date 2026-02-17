test_that("run_glm_pipeline tolerates l1 none selection", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$multi_run <- TRUE
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))
  gpa$glm_backend_specs <- list(
    fsl = list(l1_run = "run_feat_sepjobs", l3_run = "run_feat_sepjobs")
  )

  mock_job <- function(...) {
    obj <- list(...)
    obj$copy <- function(...) {
      args <- list(...)
      new <- modifyList(obj, args)
      new$copy <- obj$copy
      new
    }
    obj
  }

  mock_sequence <- function(joblist, sequence_id) {
    list(submit = function(...) NULL)
  }

  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = "none",
      l2_model_names = "l2_model1",
      l3_model_names = "l3_model1"
    ),
    R_batch_job = list(new = mock_job),
    R_batch_sequence = list(new = mock_sequence),
    .package = "fmri.pipeline"
  )
  expect_null(result)
})
