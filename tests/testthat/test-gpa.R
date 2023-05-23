test_that("Basic GPA setup works", {
  gpa <- get_gpa_minimal()
  expect_equal(1, 1)
})

test_that("Basic GPA list can be dumped to yaml with utility function", {
  gpa <- get_gpa_minimal()

  dump_to_yaml(gpa)

  expect_true(dir.exists(file.path(gpa.output_directory, "test.yaml")))
})