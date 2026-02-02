test_that("truncate_runs writes truncated NIfTI via subset_nifti_volumes", {
  skip_if_not_installed("RNifti")

  dims <- c(2, 2, 2, 15)
  img_data <- array(rnorm(prod(dims)), dim = dims)
  img <- RNifti::asNifti(img_data)
  infile <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(img, infile)

  gpa <- list(
    confound_settings = list(
      truncate_run = "volume > 12"
    ),
    drop_volumes = 0L
  )
  class(gpa) <- c("list", "glm_pipeline_arguments")

  mr_df <- data.frame(
    id = "sub1",
    session = 1L,
    run_number = 1L,
    run_nifti = infile,
    run_volumes = dims[4],
    drop_volumes = 0L,
    stringsAsFactors = FALSE
  )

  lg <- lgr::get_logger("test-truncate-runs")
  lg$set_threshold("fatal")

  truncation_data <- data.frame(volume = seq_len(dims[4]))

  out <- fmri.pipeline:::truncate_runs(
    mr_df,
    gpa = gpa,
    subj_outdir = tempdir(),
    truncation_data = truncation_data,
    lg = lg
  )

  expect_true(file.exists(out$run_nifti))
  expect_lt(out$run_volumes, mr_df$run_volumes)
  dims_out <- fmri.pipeline:::get_nifti_dim(out$run_nifti)
  expect_equal(dims_out[4], 12)
})
