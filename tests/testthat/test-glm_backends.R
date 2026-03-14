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

test_that("backend dependency helpers capture AFNI->FSL L2 requirements", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- c("afni", "spm")
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  backends <- fmri.pipeline:::get_glm_backends(gpa)

  expect_true(fmri.pipeline:::backend_l3_requires_l2(backends$afni))
  expect_identical(fmri.pipeline:::backend_l3_l2_source_backend(backends$afni), "fsl")
  expect_false(fmri.pipeline:::backend_l3_requires_l2(backends$spm))
  expect_true(fmri.pipeline:::any_l3_backend_requires_l2(backends))
  expect_identical(fmri.pipeline:::get_l3_dependency_backends(backends), "fsl")
})
