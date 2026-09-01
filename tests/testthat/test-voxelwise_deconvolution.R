test_that("the retired external Bush 2011 backend is rejected", {
  expect_error(
    voxelwise_deconvolution(
      niftis = character(0),
      TR = 1,
      algorithm = "bush2011_external"
    ),
    "one of"
  )
})
