test_that("setup_l3_models refreshes and enforces the resolved producer stage", {
  refresh_calls <- list()
  enforce_calls <- list()

  l1_models <- structure(
    list(models = list(model1 = list())),
    class = c("l1_model_set", "list")
  )
  l2_models <- structure(
    list(models = list(l2_model1 = list(l1_model_names = "model1"))),
    class = c("hi_model_set", "list")
  )
  l3_models <- structure(
    list(models = list(
      l3_model1 = list(
        execution_backend = "afni",
        producer_backend = "fsl",
        l3_input_mode = "3dlmer"
      )
    )),
    class = c("hi_model_set", "list")
  )

  gpa <- list(
    glm_software = "afni",
    multi_run = TRUE,
    lgr_threshold = "info",
    log_txt = FALSE,
    log_json = FALSE,
    output_locations = list(
      setup_l3_log_txt = file.path(tempdir(), "setup_l3_models.txt"),
      setup_l3_log_json = file.path(tempdir(), "setup_l3_models.json")
    ),
    subject_data = data.frame(
      id = paste0("sub", 1:4),
      session = 1L,
      exclude_subject = FALSE,
      stringsAsFactors = FALSE
    ),
    run_data = data.frame(
      id = rep(paste0("sub", 1:4), each = 2L),
      session = 1L,
      run_number = rep(1:2, times = 4L),
      exclude_run = FALSE,
      exclude_subject = FALSE,
      stringsAsFactors = FALSE
    ),
    l1_models = l1_models,
    l2_models = l2_models,
    l3_models = l3_models
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  testthat::local_mocked_bindings(
    normalize_longitudinal_model_signatures = function(gpa, lg) gpa,
    refresh_l1_cope_names = function(gpa, lg = NULL) gpa,
    resolve_l2_l3_compatible_pairs = function(...) {
      data.frame(l2_model = "l2_model1", l3_model = "l3_model1", stringsAsFactors = FALSE)
    },
    refresh_glm_status = function(gpa, level = 1L, lg = NULL, glm_software = NULL) {
      refresh_calls[[length(refresh_calls) + 1L]] <<- list(
        level = level,
        glm_software = glm_software
      )
      gpa
    },
    enforce_glms_complete = function(gpa, level = 1L, lg = NULL, glm_software = NULL) {
      enforce_calls[[length(enforce_calls) + 1L]] <<- list(
        level = level,
        glm_software = glm_software
      )
      invisible(NULL)
    },
    afni_3dlmer_setup = function(...) {
      list(
        metadata = NULL,
        data = data.frame(),
        id_cols = c("l1_model", "l2_model", "l3_model", "l1_cope_name", "l2_cope_name")
      )
    },
    .package = "fmri.pipeline"
  )

  res <- fmri.pipeline:::setup_l3_models(
    gpa,
    l1_model_names = "model1",
    l2_model_names = "l2_model1",
    l3_model_names = "l3_model1",
    backend = "afni"
  )

  expect_true(inherits(res$l3_model_setup, "l3_setup"))
  expect_identical(refresh_calls[[1L]]$level, 2L)
  expect_identical(refresh_calls[[1L]]$glm_software, "fsl")
  expect_identical(enforce_calls[[1L]]$level, 2L)
  expect_identical(enforce_calls[[1L]]$glm_software, "fsl")
  expect_true(any(vapply(refresh_calls, function(x) identical(x$level, 3L), logical(1))))
})
