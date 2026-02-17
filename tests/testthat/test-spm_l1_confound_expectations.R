test_that("spm_l1_model enforces per-run vs concatenated confound expectations", {
  tmp_dir <- tempfile("spm_l1_expect_")
  dir.create(tmp_dir, recursive = TRUE)
  nifti1 <- file.path(tmp_dir, "run1.nii.gz")
  nifti2 <- file.path(tmp_dir, "run2.nii.gz")
  file.create(nifti1)
  file.create(nifti2)

  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = TRUE, output_directory = tmp_dir)
  gpa$glm_software <- "spm"
  gpa$glm_settings <- list(spm = list(concatenate_runs = TRUE))
  gpa$output_locations$spm_l1_directory <- file.path(tmp_dir, "{id}", "ses-{session}", "{l1_model}")
  gpa$parallel$compute_environment <- list(global = character(0), fsl = character(0), afni = character(0), spm = character(0), r = character(0))
  gpa$l1_models$models <- list(model1 = list(name = "model1", signals = character(0)))

  d_obj <- list(run_volumes = c(10, 10), run_niftis = c(nifti1, nifti2))
  class(d_obj) <- c("bdm", "list")

  dummy_status <- function(spm_dir, lg = NULL, prefix = NULL) {
    data.frame(
      spm_dir = spm_dir,
      spm_dir_exists = TRUE,
      spm_mat_exists = FALSE,
      spm_complete = FALSE,
      spm_failed = NA,
      stringsAsFactors = FALSE
    )
  }

  testthat::local_mocked_bindings(
    generate_spm_mat = function(...) {
      list(gunzip_cmds = character(0), contrast_cmds = character(0))
    },
    get_spm_status = dummy_status,
    .env = asNamespace("fmri.pipeline")
  )

  expect_error(
    spm_l1_model(
      id = "sub1", session = 1L,
      l1_confound_files = c("a.txt", "b.txt"),
      d_obj = d_obj, gpa = gpa, model_name = "model1", run_nifti = c(nifti1, nifti2)
    ),
    "concatenated runs"
  )

  gpa$glm_settings$spm$concatenate_runs <- FALSE
  expect_error(
    spm_l1_model(
      id = "sub1", session = 1L,
      l1_confound_files = "single.txt",
      d_obj = d_obj, gpa = gpa, model_name = "model1", run_nifti = c(nifti1, nifti2)
    ),
    "must match length"
  )
})
