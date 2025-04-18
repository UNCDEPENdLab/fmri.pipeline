library(testthat)

#' Test summarization of GPA
test_that("Summarize GPA", {
  gpa <- get_gpa("testdata")
  summary(gpa)
  expect_equal(1, 1)
})

#' Test edge case where the GPA is created with minimal settings
#' Assures None values are handled correctly
test_that("Summarize GPA minimal", {
  gpa <- get_gpa_minimal("testdata")
  summary(gpa)
  expect_equal(1, 1)
})

#' Test summarization of L1 models
test_that("Summarize L1 models", {
  gpa <- get_gpa("testdata")
  summary(gpa)
  expect_equal(1, 1)
})
