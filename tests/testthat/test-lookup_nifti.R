test_that("lookup_* functions work correctly on a valid 4D NIfTI file", {
  skip_if_not_installed("RNifti")
  
  # Create synthetic 4D image with known dimensions
  dims <- c(5, 6, 7, 8)
  img_data <- array(rnorm(prod(dims)), dim = dims)
  img <- RNifti::asNifti(img_data)
  tmpfile <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(img, tmpfile)
  
  # Test lookup_run_volumes
  expect_equal(lookup_run_volumes(tmpfile), dims[4])
  
  # Test lookup_nvoxels
  expect_equal(lookup_nvoxels(tmpfile), prod(dims))
  
  # Test lookup_dim
  expected_dims <- c(dim_x = dims[1], dim_y = dims[2], dim_z = dims[3], dim_t = dims[4])
  expect_equal(lookup_dim(tmpfile), expected_dims)
})
