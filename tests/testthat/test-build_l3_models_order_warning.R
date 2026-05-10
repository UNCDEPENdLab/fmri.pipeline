test_that("build_l3_models warns when L2 models are not yet defined", {
  gpa <- create_mock_gpa(
    n_subjects = 3,
    n_runs = 2,
    include_l1_models = TRUE,
    include_l2_models = FALSE,
    include_l3_models = FALSE
  )

  spec <- list(
    l3_models = list(
      l3_age = list(
        name = "l3_age",
        model_formula = "~ age",
        num2fac = character(0),
        covariate_transform = list(),
        reference_level = list(),
        contrasts = list(diagonal = TRUE),
        l3_input_mode = "per_session"
      )
    )
  )

  spec_file <- tempfile(fileext = ".json")
  jsonlite::write_json(spec, spec_file, auto_unbox = TRUE, pretty = TRUE)

  expect_warning(
    out <- build_l3_models(gpa, from_spec_file = spec_file),
    "built before any L2 models were defined"
  )

  expect_s3_class(out$l3_models, "hi_model_set")
  expect_true("l3_age" %in% names(out$l3_models$models))
})

test_that("build_l3_models warning includes backend-specific SPM guidance", {
  gpa <- create_mock_gpa(
    n_subjects = 3,
    n_runs = 2,
    include_l1_models = TRUE,
    include_l2_models = FALSE,
    include_l3_models = FALSE
  )
  gpa$glm_software <- "spm"

  spec <- list(
    l3_models = list(
      l3_age = list(
        name = "l3_age",
        model_formula = "~ age",
        num2fac = character(0),
        covariate_transform = list(),
        reference_level = list(),
        contrasts = list(diagonal = TRUE),
        l3_input_mode = "per_session"
      )
    )
  )

  spec_file <- tempfile(fileext = ".json")
  jsonlite::write_json(spec, spec_file, auto_unbox = TRUE, pretty = TRUE)

  expect_warning(
    build_l3_models(gpa, from_spec_file = spec_file),
    "SPM, L1 setup also skips projected L2 run/session regressors"
  )
})
