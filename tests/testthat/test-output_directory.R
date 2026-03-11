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

test_that("get_output_directory routes pooled L2 models to subject-root path", {
  gpa <- create_mock_gpa()
  gpa$output_locations$feat_l2_directory <- file.path("{gpa$output_directory}", "feat_l2", "sub-{id}", "ses-{session}", "{l1_model}")
  gpa$output_locations$feat_l2_id_scope_directory <- file.path("{gpa$output_directory}", "feat_l2", "sub-{id}", "{l1_model}")

  gpa$l1_models <- list(models = list(model1 = list(name = "model1")))
  class(gpa$l1_models) <- c("l1_model_set", "list")
  gpa$l2_models <- list(models = list(
    l2_pooled = list(name = "l2_pooled", l2_scope = "id"),
    l2_session = list(name = "l2_session", l2_scope = "id_session")
  ))
  class(gpa$l2_models) <- c("hi_model_set", "list")

  pooled_dir <- get_output_directory(
    id = "sub1", session = 0L, l1_model = "model1", l2_model = "l2_pooled",
    what = "l2", gpa = gpa, glm_software = "fsl"
  )
  session_dir <- get_output_directory(
    id = "sub1", session = 2L, l1_model = "model1", l2_model = "l2_session",
    what = "l2", gpa = gpa, glm_software = "fsl"
  )

  expect_match(pooled_dir, file.path(gpa$output_directory, "feat_l2", "sub-sub1", "model1"), fixed = TRUE)
  expect_false(grepl("ses-0", pooled_dir, fixed = TRUE))
  expect_match(session_dir, file.path(gpa$output_directory, "feat_l2", "sub-sub1", "ses-2", "model1"), fixed = TRUE)
})

test_that("get_output_directory pooled L2 fallback removes ses token when id path is unset", {
  gpa <- create_mock_gpa()
  gpa$output_locations$feat_l2_directory <- file.path("{gpa$output_directory}", "feat_l2", "sub-{id}", "ses-{session}", "{l1_model}")
  gpa$output_locations$feat_l2_id_scope_directory <- NULL

  gpa$l1_models <- list(models = list(model1 = list(name = "model1")))
  class(gpa$l1_models) <- c("l1_model_set", "list")
  gpa$l2_models <- list(models = list(l2_pooled = list(name = "l2_pooled", l2_scope = "id")))
  class(gpa$l2_models) <- c("hi_model_set", "list")

  pooled_dir <- get_output_directory(
    id = "sub1", session = 0L, l1_model = "model1", l2_model = "l2_pooled",
    what = "l2", gpa = gpa, glm_software = "fsl"
  )

  expect_match(pooled_dir, file.path(gpa$output_directory, "feat_l2", "sub-sub1", "model1"), fixed = TRUE)
  expect_false(grepl("ses-0", pooled_dir, fixed = TRUE))
})

test_that("get_output_directory routes pooled SPM L1 models to subject-root path", {
  gpa <- create_mock_gpa()
  gpa$glm_settings <- list(spm = list(l1_session_mode = "pooled"))
  gpa$output_locations$spm_l1_directory <- file.path("{gpa$output_directory}", "spm_l1", "sub-{id}", "ses-{session}", "{l1_model}")
  gpa$output_locations$spm_l1_id_scope_directory <- file.path("{gpa$output_directory}", "spm_l1", "sub-{id}", "{l1_model}")
  gpa$l1_models <- list(models = list(model1 = list(name = "model1")))
  class(gpa$l1_models) <- c("l1_model_set", "list")

  pooled_dir <- get_output_directory(
    id = "sub1", session = 0L, l1_model = "model1",
    what = "l1", gpa = gpa, glm_software = "spm"
  )

  expect_match(pooled_dir, file.path(gpa$output_directory, "spm_l1", "sub-sub1", "model1"), fixed = TRUE)
  expect_false(grepl("ses-0", pooled_dir, fixed = TRUE))
})
