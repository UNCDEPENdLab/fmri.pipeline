# Tests to verify that the optimized evtmax_1 implementation
# produces correct results compared to the original per-event convolution approach.
# The new approach uses a small padded window for peak computation, which:
#   - Is faster (especially with many unique durations like RTs)
#   - Is more correct for end-of-run events (no dependence on run length)
#   - May produce tiny floating-point differences for non-boundary events due to
#     different FFT sizes (small window vs full run)

# Helper: create a reg matrix as expected by convolve_regressor
make_reg <- function(trials, onsets, durations, values, reg_name = "test_reg") {
  reg <- cbind(trial = trials, onset = onsets, duration = durations, value = values)
  attr(reg, "event") <- "cue"
  attr(reg, "reg_name") <- reg_name
  attr(reg, "physio_only") <- FALSE
  reg
}

# Reference implementation: the OLD per-event convolution algorithm
# This is an exact copy of the logic from the original convolve_regressor
# for the evtmax_1 (normeach=TRUE) path.
evtmax1_per_event_reference <- function(n_vols, times, durations, values, tr, hrf_params, convolve = TRUE) {
  normed_events <- sapply(seq_along(times), function(i) {
    stim_conv <- fmri.stimulus(
      n_vols = n_vols, values = 1.0, onsets = times[i], durations = durations[i],
      tr = tr, demean = FALSE, center_values = FALSE, convolve = convolve,
      a1 = hrf_params["a1"], a2 = hrf_params["a2"],
      b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"]
    )

    if (times[i] + durations[i] > (n_vols * tr - 20)) {
      mid_vol <- n_vols * tr / 2
      stim_at_center <- fmri.stimulus(
        n_vols = n_vols, values = 1.0, onsets = mid_vol, durations = durations[i],
        tr = tr, demean = FALSE, center_values = FALSE, convolve = convolve,
        a1 = hrf_params["a1"], a2 = hrf_params["a2"],
        b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"]
      )
      stim_conv <- stim_conv / max(stim_at_center)
    } else {
      stim_conv <- stim_conv / max(stim_conv)
    }

    stim_conv * values[i]
  })

  if (is.matrix(normed_events)) {
    rowSums(normed_events)
  } else {
    normed_events
  }
}

# ============================================================================
# Test: Small-window peak matches full-run peak for non-boundary events
# This validates that the small window doesn't introduce meaningful error
# ============================================================================

test_that("small-window HRF peak matches full-run peak", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  set.seed(42)
  test_durations <- c(0, 0.3, 0.5, 1.0, 1.5, 2.0, 3.5, 5.0, 10.0, runif(10, 0.2, 3.0))

  for (dur in test_durations) {
    # Full-run peak (old approach: event at center of run)
    peak_full <- max(fmri.stimulus(
      n_vols = n_vols, values = 1.0, onsets = n_vols * tr / 2, durations = dur,
      tr = tr, demean = FALSE, center_values = FALSE, convolve = TRUE,
      a1 = hrf_params["a1"], a2 = hrf_params["a2"],
      b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"]
    ))

    # Small-window peak (new approach)
    peak_onset <- 20
    peak_n_vols <- ceiling((peak_onset + dur + 40) / tr)
    peak_small <- max(fmri.stimulus(
      n_vols = peak_n_vols, values = 1.0, onsets = peak_onset, durations = dur,
      tr = tr, demean = FALSE, center_values = FALSE, convolve = TRUE,
      a1 = hrf_params["a1"], a2 = hrf_params["a2"],
      b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"]
    ))

    expect_equal(peak_small, peak_full, tolerance = 1e-10,
      label = sprintf("Small vs full-run peak for duration=%.4f", dur))
  }
})

# ============================================================================
# Test: HRF peak is position-independent (translation invariance)
# ============================================================================

test_that("HRF peak is position-independent within the run", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  test_durations <- c(0, 0.3, 0.5, 1.0, 1.5, 2.0, 3.5, 5.0, 10.0)
  test_onsets <- c(5, 20, 50, 80, 100, 140, 170) # avoid last 20s

  for (dur in test_durations) {
    peaks <- sapply(test_onsets, function(onset) {
      max(fmri.stimulus(
        n_vols = n_vols, values = 1.0, onsets = onset, durations = dur,
        tr = tr, demean = FALSE, center_values = FALSE, convolve = TRUE,
        a1 = hrf_params["a1"], a2 = hrf_params["a2"],
        b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"]
      ))
    })

    expect_equal(peaks, rep(peaks[1], length(peaks)), tolerance = 1e-12,
      label = sprintf("Peak invariance for duration=%.2f", dur))
  }
})

# ============================================================================
# Test: Optimized evtmax_1 vs per-event convolution (non-boundary events)
# ============================================================================

test_that("evtmax_1 optimization matches old approach: fixed duration, unit values", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  reg <- make_reg(
    trials = 1:15, onsets = seq(5, 170, length.out = 15),
    durations = rep(2.0, 15), values = rep(1, 15), reg_name = "fixed_unit"
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )
  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-8,
    label = "Fixed duration, unit values")
})

test_that("evtmax_1 optimization matches old approach: fixed duration, parametric values", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  set.seed(123)
  reg <- make_reg(
    trials = 1:15, onsets = seq(5, 170, length.out = 15),
    durations = rep(1.5, 15), values = rnorm(15), reg_name = "fixed_param"
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )
  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-8,
    label = "Fixed duration, parametric values")
})

test_that("evtmax_1 optimization matches old approach: variable (RT-like) durations", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  set.seed(456)
  n_evt <- 20
  reg <- make_reg(
    trials = 1:n_evt, onsets = seq(5, 170, length.out = n_evt),
    durations = runif(n_evt, 0.3, 2.5), values = rep(1, n_evt), reg_name = "var_unit"
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )
  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-8,
    label = "Variable RT-like durations, unit values")
})

test_that("evtmax_1 optimization matches old approach: many unique durations with parametric values", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  set.seed(789)
  n_evt <- 30
  reg <- make_reg(
    trials = 1:n_evt, onsets = sort(runif(n_evt, 5, 170)),
    durations = runif(n_evt, 0.2, 3.0), values = rnorm(n_evt), reg_name = "many_durs"
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )
  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-8,
    label = "Many unique durations, parametric values")
})

# ============================================================================
# Test: End-of-run events — new approach is more correct
# The old code used max(stim_conv) for non-boundary events but fell back to
# a center-of-run peak for events within 20s of the end. The new code uses
# a consistent small-window peak for ALL events, which is position-independent
# and gives the true peak height. Values will differ slightly for the tail
# events, but the new values are more correct.
# ============================================================================

test_that("evtmax_1 optimization handles end-of-run events correctly", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  # Include events in the last 20s of the run
  reg <- make_reg(
    trials = 1:8,
    onsets = c(10, 50, 100, 150, 175, 185, 190, 195),
    durations = c(1, 1.5, 2, 0.5, 1.2, 0.8, 2.5, 1),
    values = c(1, -0.5, 2, 0.3, -1.1, 0.8, 1.5, -0.2),
    reg_name = "end_run"
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )
  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  # Non-boundary events (onset + dur < n_vols*tr - 20 = 180s) should be very close.
  # Boundary events may differ because the new approach computes the true peak in a
  # padded context rather than using the old center-of-run fallback at full run length.
  # Use a relaxed tolerance to account for different FFT sizes in peak computation
  # while confirming the overall regressor shape is preserved.
  r <- cor(as.vector(new_result), old_result)
  expect_gt(r, 0.9999, label = "Correlation with old approach for mixed boundary/non-boundary events")

  # Also verify that a single event convolved with evtmax_1 has max <= 1.0,
  # regardless of position. This is the fundamental guarantee of evtmax_1.
  for (i in seq_len(nrow(reg))) {
    single_reg <- make_reg(
      trials = reg[i, "trial"], onsets = reg[i, "onset"],
      durations = reg[i, "duration"], values = 1.0, reg_name = "single"
    )
    single_result <- fmri.pipeline:::convolve_regressor(
      n_vols = n_vols, reg = single_reg, tr = tr, normalization = "evtmax_1",
      rm_zeros = FALSE, center_values = FALSE,
      convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
    )
    # Max should be very close to 1.0 for non-boundary events, and <= 1.0 for all
    expect_lte(max(single_result), 1.0 + 1e-10,
      label = sprintf("Peak <= 1.0 for event at onset=%.0f", reg[i, "onset"]))
  }
})

test_that("evtmax_1 optimization matches old approach: sub-second TR", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 0.5
  n_vols <- 400

  set.seed(111)
  n_evt <- 20
  reg <- make_reg(
    trials = 1:n_evt, onsets = sort(runif(n_evt, 5, 170)),
    durations = runif(n_evt, 0.2, 2.0), values = rnorm(n_evt), reg_name = "fast_tr"
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )
  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-8,
    label = "Sub-second TR with variable durations")
})

test_that("evtmax_1 optimization matches old approach: with center_values and rm_zeros", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  set.seed(222)
  n_evt <- 15
  vals <- rnorm(n_evt)
  vals[c(3, 7, 12)] <- 0

  reg <- make_reg(
    trials = 1:n_evt, onsets = seq(5, 170, length.out = n_evt),
    durations = runif(n_evt, 0.5, 2.0), values = vals, reg_name = "with_zeros"
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = TRUE, center_values = TRUE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = TRUE
  )
  ref_values <- cleaned$values
  if (length(ref_values) > 1L && sd(ref_values, na.rm = TRUE) > 1e-5) {
    ref_values <- ref_values - mean(ref_values, na.rm = TRUE)
  }
  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = ref_values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-8,
    label = "With center_values=TRUE and rm_zeros=TRUE")
})

test_that("evtmax_1 optimization matches old approach: single event", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  reg <- make_reg(
    trials = 1, onsets = 50, durations = 1.5, values = 2.3, reg_name = "single_event"
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )
  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-8,
    label = "Single event")
})

# ============================================================================
# Test: durmax_1 is unchanged by the optimization
# ============================================================================

test_that("durmax_1 normalization is unaffected by evtmax_1 optimization", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  set.seed(333)
  n_evt <- 15
  reg <- make_reg(
    trials = 1:n_evt, onsets = seq(5, 170, length.out = n_evt),
    durations = runif(n_evt, 0.5, 3.0), values = rnorm(n_evt), reg_name = "durmax_test"
  )

  hrf_boxcar <- fmri.stimulus(
    n_vols = 300 / tr, values = 1.0, onsets = 100, durations = 100, tr = tr,
    demean = FALSE,
    a1 = hrf_params["a1"], a2 = hrf_params["a2"],
    b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"]
  )
  hrf_max <- max(hrf_boxcar)

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )

  old_result <- sapply(seq_along(cleaned$times), function(i) {
    stim_conv <- fmri.stimulus(
      n_vols = n_vols, values = 1.0, onsets = cleaned$times[i], durations = cleaned$durations[i],
      tr = tr, demean = FALSE, center_values = FALSE, convolve = TRUE,
      a1 = hrf_params["a1"], a2 = hrf_params["a2"],
      b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"]
    )
    stim_conv / hrf_max * cleaned$values[i]
  })
  old_result <- rowSums(old_result)

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "durmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-10,
    label = "durmax_1 unchanged")
})

# ============================================================================
# Test: Performance — new approach is faster with many unique durations
# ============================================================================

test_that("evtmax_1 optimization is faster than per-event convolution with many unique durations", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 300

  set.seed(999)
  n_evt <- 50
  reg <- make_reg(
    trials = 1:n_evt, onsets = sort(runif(n_evt, 5, 270)),
    durations = runif(n_evt, 0.2, 3.0), # all unique RT-like durations
    values = rnorm(n_evt), reg_name = "perf_test"
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )

  t_new <- system.time({
    for (rep in 1:3) {
      fmri.pipeline:::convolve_regressor(
        n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
        rm_zeros = FALSE, center_values = FALSE,
        convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
      )
    }
  })["elapsed"]

  t_old <- system.time({
    for (rep in 1:3) {
      evtmax1_per_event_reference(
        n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
        values = cleaned$values, tr = tr, hrf_params = hrf_params
      )
    }
  })["elapsed"]

  # The new approach should be faster. The small-window peak lookups are much cheaper
  # than full-run convolutions, and the final single convolution replaces N separate ones.
  expect_lt(t_new, t_old,
    label = sprintf("New (%.3fs) should be faster than old (%.3fs)", t_new, t_old))
})

# ============================================================================
# Test: Overlapping events — additive construct_stimulus handles correctly
# ============================================================================

test_that("evtmax_1 handles overlapping events correctly via additive stimulus", {
  # Two events with long durations that overlap in time. With the old overwrite
  # semantics in construct_stimulus, the first event's contribution in the overlap
  # region would be clobbered. With additive accumulation, both contribute correctly.
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 200

  # Overlapping events: event 1 at t=20 with duration=5, event 2 at t=22 with duration=5
  # (events overlap from t=22 to t=25)
  reg <- make_reg(
    trials = 1:2, onsets = c(20, 22),
    durations = c(5, 5), values = c(1.5, 2.0)
  )

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )

  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-10,
    label = "overlapping events match per-event reference")

  # Also verify the combined result differs from a non-overlapping scenario
  # (i.e., the overlap is actually being handled, not just ignored)
  reg_nooverlap <- make_reg(
    trials = 1:2, onsets = c(20, 30),
    durations = c(5, 5), values = c(1.5, 2.0)
  )
  nooverlap_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg_nooverlap, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  expect_false(isTRUE(all.equal(as.vector(new_result), as.vector(nooverlap_result))),
    label = "overlapping vs non-overlapping should differ")
})

test_that("evtmax_1 handles many overlapping events with different durations and values", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
  tr <- 1.0
  n_vols <- 300

  # Dense events with RT-like durations — many will overlap
  set.seed(42)
  n_evt <- 20
  onsets <- sort(runif(n_evt, 10, 100))  # dense enough to cause overlaps
  durations <- runif(n_evt, 1.0, 6.0)
  values <- rnorm(n_evt, mean = 3, sd = 1)

  reg <- make_reg(trials = 1:n_evt, onsets = onsets, durations = durations, values = values)

  new_result <- fmri.pipeline:::convolve_regressor(
    n_vols = n_vols, reg = reg, tr = tr, normalization = "evtmax_1",
    rm_zeros = FALSE, center_values = FALSE,
    convmax_1 = FALSE, demean_convolved = FALSE, hrf_parameters = hrf_params
  )

  cleaned <- fmri.pipeline:::cleanup_regressor(
    reg[, "onset"], reg[, "duration"], reg[, "value"], rm_zeros = FALSE
  )

  old_result <- evtmax1_per_event_reference(
    n_vols = n_vols, times = cleaned$times, durations = cleaned$durations,
    values = cleaned$values, tr = tr, hrf_params = hrf_params
  )

  expect_equal(as.vector(new_result), old_result, tolerance = 1e-10,
    label = "many overlapping events match per-event reference")
})

# ============================================================================
# Test: R and C++ convolution paths produce identical results
# ============================================================================

test_that("R and C++ convolution paths produce identical results for impulse event", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)

  r_result <- fmri.stimulus(
    n_vols = 200, values = 1.0, onsets = 50, durations = 0,
    tr = 1.0, demean = FALSE, center_values = FALSE, convolve = TRUE,
    a1 = hrf_params["a1"], a2 = hrf_params["a2"],
    b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"],
    conv_method = "r"
  )

  cpp_result <- fmri.stimulus(
    n_vols = 200, values = 1.0, onsets = 50, durations = 0,
    tr = 1.0, demean = FALSE, center_values = FALSE, convolve = TRUE,
    a1 = hrf_params["a1"], a2 = hrf_params["a2"],
    b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"],
    conv_method = "cpp"
  )

  expect_equal(as.vector(r_result), as.vector(cpp_result), tolerance = 1e-12,
    label = "R and C++ paths match for impulse event")
})

test_that("R and C++ convolution paths match across durations and TRs", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)

  test_cases <- list(
    list(dur = 0.5, tr = 1.0, n_vols = 200, onset = 50),
    list(dur = 2.0, tr = 1.0, n_vols = 200, onset = 50),
    list(dur = 10.0, tr = 1.0, n_vols = 200, onset = 50),
    list(dur = 1.0, tr = 0.5, n_vols = 400, onset = 50),
    list(dur = 1.0, tr = 2.0, n_vols = 150, onset = 50),
    list(dur = 0, tr = 0.72, n_vols = 300, onset = 30)
  )

  for (tc in test_cases) {
    r_result <- fmri.stimulus(
      n_vols = tc$n_vols, values = 1.0, onsets = tc$onset, durations = tc$dur,
      tr = tc$tr, demean = FALSE, center_values = FALSE, convolve = TRUE,
      a1 = hrf_params["a1"], a2 = hrf_params["a2"],
      b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"],
      conv_method = "r"
    )

    cpp_result <- fmri.stimulus(
      n_vols = tc$n_vols, values = 1.0, onsets = tc$onset, durations = tc$dur,
      tr = tc$tr, demean = FALSE, center_values = FALSE, convolve = TRUE,
      a1 = hrf_params["a1"], a2 = hrf_params["a2"],
      b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"],
      conv_method = "cpp"
    )

    label <- sprintf("dur=%.1f, tr=%.2f, n_vols=%d", tc$dur, tc$tr, tc$n_vols)
    expect_equal(as.vector(r_result), as.vector(cpp_result), tolerance = 1e-12, label = label)
  }
})

test_that("R and C++ convolution paths match with multiple events", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)

  onsets <- c(10, 25, 50, 80, 120)
  durations <- c(1, 2, 0.5, 3, 1.5)
  values <- c(1.0, 0.5, 2.0, 0.8, 1.2)

  r_result <- fmri.stimulus(
    n_vols = 200, values = values, onsets = onsets, durations = durations,
    tr = 1.0, demean = FALSE, center_values = FALSE, convolve = TRUE,
    a1 = hrf_params["a1"], a2 = hrf_params["a2"],
    b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"],
    conv_method = "r"
  )

  cpp_result <- fmri.stimulus(
    n_vols = 200, values = values, onsets = onsets, durations = durations,
    tr = 1.0, demean = FALSE, center_values = FALSE, convolve = TRUE,
    a1 = hrf_params["a1"], a2 = hrf_params["a2"],
    b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"],
    conv_method = "cpp"
  )

  expect_equal(as.vector(r_result), as.vector(cpp_result), tolerance = 1e-12,
    label = "R and C++ match with multiple events of varying height/duration")
})

test_that("R and C++ convolution paths match for end-of-run events", {
  hrf_params <- c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)

  # Event near the end of the run
  r_result <- fmri.stimulus(
    n_vols = 200, values = 1.0, onsets = 195, durations = 2.0,
    tr = 1.0, demean = FALSE, center_values = FALSE, convolve = TRUE,
    a1 = hrf_params["a1"], a2 = hrf_params["a2"],
    b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"],
    conv_method = "r"
  )

  cpp_result <- fmri.stimulus(
    n_vols = 200, values = 1.0, onsets = 195, durations = 2.0,
    tr = 1.0, demean = FALSE, center_values = FALSE, convolve = TRUE,
    a1 = hrf_params["a1"], a2 = hrf_params["a2"],
    b1 = hrf_params["b1"], b2 = hrf_params["b2"], cc = hrf_params["cc"],
    conv_method = "cpp"
  )

  expect_equal(as.vector(r_result), as.vector(cpp_result), tolerance = 1e-12,
    label = "R and C++ match for end-of-run event")
})
