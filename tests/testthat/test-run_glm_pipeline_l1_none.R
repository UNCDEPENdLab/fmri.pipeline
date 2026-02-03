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
    fsl = list(l1_run = "fake_l1", l3_run = "fake_l3")
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

  choose <- function(...) {
    list(
      l1_model_names = NULL,
      l2_model_names = "l2int",
      l3_model_names = "l3int",
      l1_model_name = NULL,
      l2_model_name = "l2int",
      l3_model_name = "l3int"
    )
  }

  get_backends <- function(...) {
    list(fsl = list(l1_run = function() {}, l3_run = function() {}))
  }

  default_specs <- function() gpa$glm_backend_specs

  result <- testthat::with_mocked_bindings(
    run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = "l2int", l3_model_names = "l3int"),
    choose_glm_set = choose,
    get_glm_backends = get_backends,
    default_glm_backend_specs = default_specs,
    R_batch_job = list(new = mock_job),
    R_batch_sequence = list(new = mock_sequence)
  )
  expect_null(result)
})
