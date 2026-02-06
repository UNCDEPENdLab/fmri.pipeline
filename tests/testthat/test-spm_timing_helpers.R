test_that("infer_spm_microtime reads SliceTiming metadata", {
  tmp_dir <- tempfile("spm_timing_")
  dir.create(tmp_dir, recursive = TRUE)

  run_nifti <- file.path(tmp_dir, "sub-01_task-test_run-01_bold.nii.gz")
  json_path <- sub("\\.nii\\.gz$", ".json", run_nifti)

  jsonlite::write_json(
    list(SliceTiming = c(0, 0.5, 1.0, 1.5), ReferenceSlice = 3),
    json_path,
    auto_unbox = TRUE
  )

  res <- fmri.pipeline:::infer_spm_microtime(run_nifti, fmri_t = NULL, fmri_t0 = NULL, lg = NULL)
  expect_equal(res$fmri_t, 4L)
  expect_equal(res$fmri_t0, 3L)
})

test_that("infer_spm_microtime defaults fmri_t0 when no reference slice found", {
  tmp_dir <- tempfile("spm_timing_noref_")
  dir.create(tmp_dir, recursive = TRUE)

  run_nifti <- file.path(tmp_dir, "sub-01_task-test_run-01_bold.nii.gz")
  json_path <- sub("\\.nii\\.gz$", ".json", run_nifti)

  jsonlite::write_json(
    list(SliceTiming = c(0, 0.5, 1.0, 1.5)),
    json_path,
    auto_unbox = TRUE
  )

  res <- fmri.pipeline:::infer_spm_microtime(run_nifti, fmri_t = NULL, fmri_t0 = NULL, lg = NULL)
  expect_equal(res$fmri_t, 4L)
  expect_equal(res$fmri_t0, 1L)
})

test_that("infer_spm_microtime reads .stimes files when JSON is absent", {
  tmp_dir <- tempfile("spm_timing_stimes_")
  dir.create(tmp_dir, recursive = TRUE)

  run_nifti <- file.path(tmp_dir, "sub-01_task-test_run-01_bold.nii.gz")
  stimes_path <- file.path(tmp_dir, ".stimes")
  writeLines("0,0.5,1.0,1.5", stimes_path)

  res <- fmri.pipeline:::infer_spm_microtime(run_nifti, fmri_t = NULL, fmri_t0 = NULL, lg = NULL)
  expect_equal(res$fmri_t, 4L)
  expect_equal(res$fmri_t0, 1L)
})
