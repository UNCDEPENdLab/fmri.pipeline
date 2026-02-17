make_mock_logger <- function() {
  list(
    set_threshold = function(...) NULL,
    info = function(...) NULL,
    warn = function(...) NULL,
    debug = function(...) NULL,
    error = function(...) NULL
  )
}

make_mock_job_new <- function() {
  function(...) {
    obj <- list(...)
    obj$depends_on_parents <- NULL
    obj$wait_for_children <- FALSE
    obj$input_rdata_file <- NULL
    obj$output_rdata_file <- NULL
    obj$copy <- function(...) {
      args <- list(...)
      new <- modifyList(obj, args)
      new$copy <- obj$copy
      new
    }
    obj
  }
}

make_mock_sequence_new <- function(store_env) {
  function(joblist, sequence_id) {
    store_env$joblist <- joblist
    store_env$sequence_id <- sequence_id
    list(submit = function(...) NULL)
  }
}

flatten_jobs <- function(joblist) {
  out <- list()
  for (item in joblist) {
    if (is.null(item)) next
    if (is.list(item) && !is.null(item$job_name)) {
      out <- c(out, list(item))
    } else if (is.list(item)) {
      out <- c(out, flatten_jobs(item))
    }
  }
  out
}

make_model_list <- function(l1 = NULL, l2 = NULL, l3 = NULL) {
  list(
    l1_model_names = l1,
    l2_model_names = l2,
    l3_model_names = l3,
    l1_model_name = l1,
    l2_model_name = l2,
    l3_model_name = l3
  )
}

test_that("L1 none + L3 requested waits for split_backend_caches", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- FALSE
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

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = NULL, l3_model_names = "l3_model1"),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  expect_true(length(jobs) > 0L, info = "No jobs were captured; check backend resolution and mock binding setup.")
  split_job <- jobs[vapply(jobs, function(j) j$job_name == "split_backend_caches", logical(1))]
  l3_job <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l3_fsl", logical(1))]

  expect_null(result)
  expect_length(split_job, 1L)
  expect_length(l3_job, 1L)
  expect_equal(split_job[[1]]$depends_on_parents, "finalize_configuration")
  expect_equal(l3_job[[1]]$depends_on_parents, "split_backend_caches")
})

test_that("L1 none + L2 requested waits for split_backend_caches", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
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

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = "l2_model1", l3_model_names = "none"),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  expect_true(length(jobs) > 0L, info = "No jobs were captured; check backend resolution and mock binding setup.")
  split_job <- jobs[vapply(jobs, function(j) j$job_name == "split_backend_caches", logical(1))]
  l2_job <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l2", logical(1))]

  expect_null(result)
  expect_length(split_job, 1L)
  expect_length(l2_job, 1L)
  expect_equal(l2_job[[1]]$depends_on_parents, "split_backend_caches")
})

test_that("L2 requested but cannot run completes without L2 job", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- FALSE
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

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = "l2_model1", l3_model_names = "none"),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  l2_jobs <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l2", logical(1))]

  expect_null(result)
  expect_length(l2_jobs, 0L)
})

test_that("No L3 requested completes without L3 job", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
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

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = "l2_model1", l3_model_names = "none"),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  l3_jobs <- jobs[vapply(jobs, function(j) grepl("^setup_run_l3_", j$job_name), logical(1))]

  expect_null(result)
  expect_length(l3_jobs, 0L)
})
