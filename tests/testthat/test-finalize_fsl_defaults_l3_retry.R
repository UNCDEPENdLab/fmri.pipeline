test_that("finalize_pipeline_configuration defaults L3 outlier retry toggle to TRUE", {
  finalize_path <- source_tree_file("R", "finalize_pipeline_configuration.R")
  lines <- readLines(finalize_path, warn = FALSE)

  expect_true(any(grepl("auto_retry_l3_excessive_outliers\\s*=\\s*TRUE", lines)))
})

test_that("finalize_pipeline_configuration defaults SPM concatenate_runs to FALSE", {
  finalize_path <- source_tree_file("R", "finalize_pipeline_configuration.R")
  lines <- readLines(finalize_path, warn = FALSE)

  expect_true(any(grepl("concatenate_runs\\s*=\\s*FALSE", lines)))
  expect_true(any(grepl("l1_session_mode\\s*=\\s*\"separate\"", lines)))
  expect_true(any(grepl("l2_projection_interaction_contrast_modes\\s*=\\s*NULL", lines)))
})

test_that("multi-run FEAT L3 AFNI default groups L2 copes as sub-briks", {
  gpa <- create_mock_gpa()
  gpa$multi_run <- TRUE

  lg <- lgr::get_logger("test_finalize_defaults")
  lg$set_threshold("fatal")
  gpa <- fmri.pipeline:::setup_output_locations(gpa, lg)

  expect_match(gpa$output_locations$feat_l3_combined_filename, "L2m-\\{l2_model\\}")
  expect_match(gpa$output_locations$feat_l3_combined_filename, "l1c-\\{l1_cope_name\\}")
  expect_no_match(gpa$output_locations$feat_l3_combined_filename, "l2_cope_name")
  expect_match(gpa$output_locations$feat_l3_combined_briknames, "\\{l2_cope_name\\}")
  expect_match(gpa$output_locations$feat_l3_combined_briknames, "\\{l3_cope_name\\}")
})
