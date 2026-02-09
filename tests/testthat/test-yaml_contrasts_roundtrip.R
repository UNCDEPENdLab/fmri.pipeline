test_that("YAML contrasts are nested and roundtrip into contrast_spec", {
  skip_if_not_installed("yaml")

  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 1, n_trials = 5, output_directory = tempdir())

  spec <- list(
    onsets = list("onset"),
    durations = list("duration"),
    events = list(
      stim = list(onset = "onset", duration = "duration")
    ),
    signals = list(
      stim = list(event = "stim", normalization = "none", value_fixed = 1)
    ),
    l1_models = list(
      basic = list(
        signals = list("stim"),
        contrasts = list(diagonal = TRUE)
      )
    )
  )

  spec_file <- tempfile(fileext = ".yaml")
  writeLines(yaml::as.yaml(spec), spec_file)

  gpa <- build_l1_models(gpa, from_spec_file = spec_file)
  mm <- gpa$l1_models$models$basic
  expect_true(isTRUE(mm$contrast_spec$diagonal))

  out_file <- tempfile(fileext = ".yaml")
  l1_cfg <- fmri.pipeline:::get_l1_config(gpa)
  writeLines(yaml::as.yaml(l1_cfg), out_file)
  out_list <- yaml::read_yaml(out_file)

  expect_true("contrasts" %in% names(out_list$l1_models$basic))
  expect_true(isTRUE(out_list$l1_models$basic$contrasts$diagonal))
  expect_false("diagonal" %in% names(out_list$l1_models$basic))

  gpa2 <- create_mock_gpa(n_subjects = 1, n_runs = 1, n_trials = 5, output_directory = tempdir())
  gpa2 <- build_l1_models(gpa2, from_spec_file = out_file)
  mm2 <- gpa2$l1_models$models$basic
  expect_true(isTRUE(mm2$contrast_spec$diagonal))
})
