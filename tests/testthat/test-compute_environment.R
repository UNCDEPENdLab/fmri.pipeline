make_fake_executable <- function(bin_dir, name, body = "echo fake") {
  path <- file.path(bin_dir, name)
  writeLines(c("#!/bin/sh", body), path)
  Sys.chmod(path, "0755")
  path
}

make_compute_env_gpa <- function(bin_dir, l3_input_mode = NULL) {
  gpa <- create_mock_gpa(include_l3_models = TRUE)
  gpa$glm_software <- "afni"
  gpa$glm_settings <- list(afni = list())
  gpa$parallel$compute_environment <- list(
    global = paste("export PATH=", shQuote(bin_dir), sep = ""),
    afni = character(0),
    fsl = character(0),
    spm = character(0),
    r = character(0)
  )
  gpa$l3_models$models$l3_model1$l3_input_mode <- l3_input_mode
  gpa
}

test_that("test_compute_environment does not require 3dLMEr for generic AFNI models", {
  bin_dir <- tempfile("fake_afni_bin_")
  dir.create(bin_dir)
  make_fake_executable(bin_dir, "afni", "echo afni")

  gpa <- make_compute_env_gpa(bin_dir, l3_input_mode = "per_session")

  expect_no_error(
    fmri.pipeline:::test_compute_environment(gpa, what = "afni", stop_on_fail = TRUE)
  )
})

test_that("test_compute_environment requires 3dLMEr only for 3dlmer model sets", {
  bin_dir <- tempfile("fake_3dlmer_bin_")
  dir.create(bin_dir)
  make_fake_executable(bin_dir, "afni", "echo afni")
  make_fake_executable(bin_dir, "R", "echo R")
  make_fake_executable(bin_dir, "Rscript", "exit 0")

  gpa <- make_compute_env_gpa(bin_dir, l3_input_mode = "3dlmer")

  expect_error(
    fmri.pipeline:::test_compute_environment(gpa, what = "afni", stop_on_fail = TRUE),
    "compute environment"
  )

  make_fake_executable(bin_dir, "3dLMEr", "echo 3dLMEr")
  expect_no_error(
    fmri.pipeline:::test_compute_environment(gpa, what = "afni", stop_on_fail = TRUE)
  )
})
