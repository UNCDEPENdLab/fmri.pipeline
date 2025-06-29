library(testthat)
library(fmri)

# Define reusable variables
TR <- 2
n_vols <- 100
onsets <- seq(1, 80, by = 5)

# Unit tests for fmri.stimulus

# 1. Basic functionality - volume units
test_that("fmri.stimulus returns correct length for volume units", {
  stim <- fmri.stimulus(n_vols = n_vols, onsets = onsets, durations = 1,
                        units = "volumes", tr = TR, convolve = TRUE)
  expect_equal(length(stim), n_vols)
})

# 2. Basic functionality - time units
test_that("fmri.stimulus returns correct length for time units", {
  stim <- fmri.stimulus(n_vols = n_vols, onsets = TR * (onsets - 1), durations = TR,
                        units = "time", tr = TR, convolve = TRUE)
  expect_equal(length(stim), n_vols)
})

# 3. Equivalence of volume- and time-based onsets
test_that("volume- and time-based onsets yield similar results", {
  stim_vol <- fmri.stimulus(n_vols = n_vols, onsets = onsets, durations = 1,
                            units = "volumes", tr = TR, convolve = TRUE)
  stim_time <- fmri.stimulus(n_vols = n_vols, onsets = TR * (onsets - 1), durations = TR,
                             units = "time", tr = TR, convolve = TRUE)
  expect_gt(cor(stim_vol, stim_time), 0.99)
})

# 4. PPI functionality with ts_multiplier
test_that("fmri.stimulus handles ts_multiplier with convolution", {
  tsvec <- (1:n_vols) / 90
  stim <- fmri.stimulus(n_vols = n_vols, onsets = TR * (onsets - 1), durations = TR,
                        units = "time", tr = TR, ts_multiplier = tsvec, convolve = TRUE)
  expect_equal(length(stim), n_vols)
  expect_true(any(stim != 0))
})

# 5. Unconvolved stimulus
test_that("fmri.stimulus returns correct unconvolved stimulus", {
  stim <- fmri.stimulus(n_vols = n_vols, onsets = TR * (onsets - 1), durations = TR,
                        units = "time", tr = TR, convolve = FALSE)
  expect_equal(length(stim), n_vols)
  expect_equal(length(which(stim != 0)), length(onsets) * TR)
})

# 6. ts_multiplier only
test_that("fmri.stimulus works with ts_multiplier only", {
  tsvec <- (1:n_vols) / 90
  stim <- fmri.stimulus(n_vols = n_vols, tr = TR, ts_multiplier = tsvec, convolve = TRUE)
  expect_equal(length(stim), n_vols)
  expect_true(any(stim != 0))
})