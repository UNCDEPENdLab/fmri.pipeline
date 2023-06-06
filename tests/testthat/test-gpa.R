test_that("Basic GPA setup works", {
  gpa <- get_gpa_minimal()
  expect_equal(1, 1)
})

test_that("Summarize GPA setup", {
  gpa <- get_gpa_minimal()
  summarize_l1_models(gpa)
  expect_equal(1, 1)
})

test_that("Summarize GPA", {
  gpa <- get_gpa_minimal()
  summarize_pipeline(gpa)
  expect_equal(1, 1)
})
