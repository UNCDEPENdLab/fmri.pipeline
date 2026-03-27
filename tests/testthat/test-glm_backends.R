# Tests for GLM backend registry helpers

test_that("get_glm_backends returns resolved backends from gpa", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- c("fsl", "spm")
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  backends <- fmri.pipeline:::get_glm_backends(gpa)

  expect_type(backends, "list")
  expect_true(all(c("fsl", "spm") %in% names(backends)))
  expect_true(is.function(backends$fsl$l1_setup))
  expect_true(is.function(backends$spm$l1_setup))
})
