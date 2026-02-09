# test output directory selection for sub/ses paths

test_that("get_output_directory uses software-specific sub/ses directories", {
  gpa <- create_mock_gpa()
  gpa$output_locations$feat_sub_directory <- file.path("{gpa$output_directory}", "feat_l1", "sub-{id}")
  gpa$output_locations$feat_ses_directory <- gpa$output_locations$feat_sub_directory
  gpa$output_locations$spm_sub_directory <- file.path("{gpa$output_directory}", "spm_l1", "sub-{id}")
  gpa$output_locations$afni_sub_directory <- file.path("{gpa$output_directory}", "afni_l1", "sub-{id}")

  fsl_sub <- get_output_directory(id = "sub1", session = 1L, run_number = 1L, gpa = gpa, what = "sub", glm_software = "fsl")
  spm_sub <- get_output_directory(id = "sub1", session = 1L, run_number = 1L, gpa = gpa, what = "sub", glm_software = "spm")
  afni_sub <- get_output_directory(id = "sub1", session = 1L, run_number = 1L, gpa = gpa, what = "sub", glm_software = "afni")

  expect_match(fsl_sub, file.path(gpa$output_directory, "feat_l1", "sub-sub1"), fixed = TRUE)
  expect_match(spm_sub, file.path(gpa$output_directory, "spm_l1", "sub-sub1"), fixed = TRUE)
  expect_match(afni_sub, file.path(gpa$output_directory, "afni_l1", "sub-sub1"), fixed = TRUE)

  fsl_ses <- get_output_directory(id = "sub1", session = 1L, run_number = 1L, gpa = gpa, what = "ses", glm_software = "fsl")
  spm_ses <- get_output_directory(id = "sub1", session = 1L, run_number = 1L, gpa = gpa, what = "ses", glm_software = "spm")
  afni_ses <- get_output_directory(id = "sub1", session = 1L, run_number = 1L, gpa = gpa, what = "ses", glm_software = "afni")

  expect_match(fsl_ses, file.path(gpa$output_directory, "feat_l1", "sub-sub1"), fixed = TRUE)
  expect_match(spm_ses, file.path(gpa$output_directory, "spm_l1", "sub-sub1"), fixed = TRUE)
  expect_match(afni_ses, file.path(gpa$output_directory, "afni_l1", "sub-sub1"), fixed = TRUE)
})
