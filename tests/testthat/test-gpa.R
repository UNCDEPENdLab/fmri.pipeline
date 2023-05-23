test_that("Basic GPA setup works", {
  gpa <- get_gpa_minimal()
  expect_equal(1, 1)
})

test_that("Summarize GPA setup", {
  gpa <- get_gpa()
  summarize_l1_models(gpa)
})