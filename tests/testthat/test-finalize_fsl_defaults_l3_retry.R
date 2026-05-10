test_that("finalize_pipeline_configuration defaults L3 outlier retry toggle to TRUE", {
  finalize_path <- normalizePath(file.path(pkg_root, "R", "finalize_pipeline_configuration.R"), mustWork = TRUE)
  lines <- readLines(finalize_path, warn = FALSE)

  expect_true(any(grepl("auto_retry_l3_excessive_outliers\\s*=\\s*TRUE", lines)))
})

test_that("finalize_pipeline_configuration defaults SPM concatenate_runs to FALSE", {
  finalize_path <- normalizePath(file.path(pkg_root, "R", "finalize_pipeline_configuration.R"), mustWork = TRUE)
  lines <- readLines(finalize_path, warn = FALSE)

  expect_true(any(grepl("concatenate_runs\\s*=\\s*FALSE", lines)))
  expect_true(any(grepl("l1_session_mode\\s*=\\s*\"separate\"", lines)))
  expect_true(any(grepl("l2_projection_interaction_contrast_modes\\s*=\\s*NULL", lines)))
})
