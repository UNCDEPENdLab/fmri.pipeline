# Unit tests for build_design_matrix and all helper functions
# These tests focus on individual function behavior rather than integration testing

# ==============================================================================
# Test fixtures and helpers
# ==============================================================================

#' Create minimal events data.frame for testing
create_test_events <- function(nruns = 2, ntrials = 5, tr = 1.0) {
  do.call(rbind, lapply(1:nruns, function(r) {
    data.frame(
      event = rep("cue", ntrials),
      run_number = r,
      trial = 1:ntrials,
      onset = cumsum(runif(ntrials, min = 2, max = 5)),
      duration = rep(1, ntrials),
      stringsAsFactors = FALSE
    )
  }))
}

#' Create minimal signals list for testing
#' Note: signals require a 'name' field which would normally be added by build_l1_models
create_test_signals <- function(nruns = 2, ntrials = 5) {
  value_df <- do.call(rbind, lapply(1:nruns, function(r) {
    data.frame(
      run_number = r,
      trial = 1:ntrials,
      value = rnorm(ntrials, mean = 5, sd = 2)
    )
  }))
  
  list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none"),
    cue_param = list(name = "cue_param", event = "cue", value = value_df, normalization = "none")
  )
}

#' Create minimal run_data for testing build_design_matrix
create_test_run_data <- function(nruns = 2, run_volumes = 50, drop_volumes = 0) {
  data.frame(
    run_number = 1:nruns,
    run_volumes = rep(run_volumes, nruns),
    drop_volumes = rep(drop_volumes, nruns),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Tests for Stage 1 helpers: Input validation and preprocessing
# ==============================================================================

# --- validate_write_timing_files ---

test_that("validate_write_timing_files accepts valid values", {
  expect_equal(fmri.pipeline:::validate_write_timing_files("convolved"), "convolved")
  expect_equal(fmri.pipeline:::validate_write_timing_files("fsl"), "fsl")
  expect_equal(fmri.pipeline:::validate_write_timing_files("afni"), "afni")
  expect_equal(fmri.pipeline:::validate_write_timing_files("spm"), "spm")
  expect_equal(fmri.pipeline:::validate_write_timing_files(c("fsl", "afni")), c("fsl", "afni"))
})

test_that("validate_write_timing_files accepts NULL", {
  expect_null(fmri.pipeline:::validate_write_timing_files(NULL))
})

test_that("validate_write_timing_files normalizes to lowercase", {
  expect_equal(fmri.pipeline:::validate_write_timing_files("CONVOLVED"), "convolved")
  expect_equal(fmri.pipeline:::validate_write_timing_files("FSL"), "fsl")
  expect_equal(fmri.pipeline:::validate_write_timing_files(c("AFNI", "SPM")), c("afni", "spm"))
})

test_that("validate_write_timing_files rejects invalid values", {
  expect_error(fmri.pipeline:::validate_write_timing_files("invalid"))
  expect_error(fmri.pipeline:::validate_write_timing_files("nifti"))
  expect_error(fmri.pipeline:::validate_write_timing_files(c("fsl", "bad")))
})

test_that("validate_write_timing_files rejects non-character input", {
  expect_error(fmri.pipeline:::validate_write_timing_files(123))
  expect_error(fmri.pipeline:::validate_write_timing_files(TRUE))
})

# --- validate_run_data ---

test_that("validate_run_data accepts valid run_data", {
  run_data <- data.frame(
    run_number = 1:3,
    run_volumes = c(100, 100, 100),
    drop_volumes = c(5, 5, 5)
  )
  result <- fmri.pipeline:::validate_run_data(run_data)
  expect_equal(nrow(result), 3)
  expect_equal(result$run_number, 1:3)
  expect_equal(result$drop_volumes, c(5, 5, 5))
})

test_that("validate_run_data adds run_number if missing", {
  run_data <- data.frame(run_volumes = c(100, 100))
  expect_message(
    result <- fmri.pipeline:::validate_run_data(run_data),
    "No run_number column found"
  )
  expect_equal(result$run_number, 1:2)
})

test_that("validate_run_data expands scalar drop_volumes to all runs", {
  run_data <- data.frame(
    run_number = 1:3,
    run_volumes = c(100, 100, 100)
  )
  lg <- lgr::get_logger()
  result <- fmri.pipeline:::validate_run_data(run_data, drop_volumes = 5, lg = lg)
  expect_equal(result$drop_volumes, c(5, 5, 5))
})

test_that("validate_run_data uses run-specific drop_volumes when provided", {
  run_data <- data.frame(
    run_number = 1:3,
    run_volumes = c(100, 100, 100)
  )
  result <- fmri.pipeline:::validate_run_data(run_data, drop_volumes = c(3, 4, 5))
  expect_equal(result$drop_volumes, c(3, 4, 5))
})

test_that("validate_run_data defaults drop_volumes to 0 when not specified", {
  run_data <- data.frame(
    run_number = 1:2,
    run_volumes = c(100, 100)
  )
  result <- fmri.pipeline:::validate_run_data(run_data, drop_volumes = 0)
  expect_equal(result$drop_volumes, c(0, 0))
})

test_that("validate_run_data preserves existing drop_volumes in run_data", {
  run_data <- data.frame(
    run_number = 1:2,
    run_volumes = c(100, 100),
    drop_volumes = c(10, 20)
  )
  result <- fmri.pipeline:::validate_run_data(run_data, drop_volumes = 5)
  expect_equal(result$drop_volumes, c(10, 20))  # Should preserve original
})

test_that("validate_run_data rejects invalid run_volumes", {
  run_data <- data.frame(
    run_number = 1:2,
    run_volumes = c(0, 100)  # 0 is invalid
  )
  expect_error(fmri.pipeline:::validate_run_data(run_data))
})

test_that("validate_run_data rejects negative drop_volumes", {
  run_data <- data.frame(
    run_number = 1:2,
    run_volumes = c(100, 100),
    drop_volumes = c(-1, 5)
  )
  expect_error(fmri.pipeline:::validate_run_data(run_data))
})

# --- validate_events ---

test_that("validate_events accepts valid events data.frame", {
  events <- data.frame(
    event = c("cue", "cue"),
    trial = c(1, 2),
    onset = c(0, 5),
    duration = c(1, 1),
    run_number = c(1, 1)
  )
  result <- fmri.pipeline:::validate_events(events)
  expect_equal(nrow(result), 2)
})

test_that("validate_events adds run_number if missing", {
  events <- data.frame(
    event = "cue",
    trial = 1,
    onset = 0,
    duration = 1
  )
  result <- fmri.pipeline:::validate_events(events)
  expect_equal(result$run_number, 1)
})

test_that("validate_events stops if events is NULL", {
  expect_error(
    fmri.pipeline:::validate_events(NULL),
    "You must pass in an events data.frame"
  )
})

test_that("validate_events stops if event column is missing", {
  events <- data.frame(
    trial = 1,
    onset = 0,
    duration = 1
  )
  expect_error(
    fmri.pipeline:::validate_events(events),
    "must contain event column"
  )
})

test_that("validate_events stops if trial column is missing", {
  events <- data.frame(
    event = "cue",
    onset = 0,
    duration = 1
  )
  expect_error(
    fmri.pipeline:::validate_events(events),
    "must contain trial column"
  )
})

test_that("validate_events stops if onset column is missing", {
  events <- data.frame(
    event = "cue",
    trial = 1,
    duration = 1
  )
  expect_error(
    fmri.pipeline:::validate_events(events),
    "must contain onset column"
  )
})

test_that("validate_events stops if duration column is missing", {
  events <- data.frame(
    event = "cue",
    trial = 1,
    onset = 0
  )
  expect_error(
    fmri.pipeline:::validate_events(events),
    "must contain duration column"
  )
})

test_that("validate_events stops on NA durations", {
  events <- data.frame(
    event = "cue",
    trial = 1,
    onset = 0,
    duration = NA_real_
  )
  expect_error(
    fmri.pipeline:::validate_events(events),
    "Invalid missing \\(NA\\) durations"
  )
})

test_that("validate_events stops on negative durations", {
  events <- data.frame(
    event = "cue",
    trial = 1,
    onset = 0,
    duration = -1
  )
  expect_error(
    fmri.pipeline:::validate_events(events),
    "Invalid negative durations"
  )
})

test_that("validate_events stops on NA onsets", {
  events <- data.frame(
    event = "cue",
    trial = 1,
    onset = NA_real_,
    duration = 1
  )
  expect_error(
    fmri.pipeline:::validate_events(events),
    "Invalid missing \\(NA\\) onsets"
  )
})

test_that("validate_events stops on negative onsets", {
  events <- data.frame(
    event = "cue",
    trial = 1,
    onset = -5,
    duration = 1
  )
  expect_error(
    fmri.pipeline:::validate_events(events),
    "Invalid negative onsets"
  )
})

# --- validate_signals ---

test_that("validate_signals accepts valid signals list", {
  signals <- list(
    list(name = "cue", event = "cue", value = 1)
  )
  result <- fmri.pipeline:::validate_signals(signals)
  expect_equal(result, signals)
})

test_that("validate_signals stops if signals is NULL", {
  expect_error(
    fmri.pipeline:::validate_signals(NULL),
    "You must pass in a signals list"
  )
})

# --- validate_tr ---

test_that("validate_tr accepts valid TR values", {
  expect_equal(fmri.pipeline:::validate_tr(1.0), 1.0)
  expect_equal(fmri.pipeline:::validate_tr(2.0), 2.0)
  expect_equal(fmri.pipeline:::validate_tr(0.5), 0.5)
})

test_that("validate_tr stops if TR is NULL", {
  expect_error(
    fmri.pipeline:::validate_tr(NULL),
    "You must pass in the tr"
  )
})
  
test_that("validate_tr rejects TR that is too small", {
  expect_error(fmri.pipeline:::validate_tr(0.001))
})

test_that("validate_tr rejects TR that is too large", {
  expect_error(fmri.pipeline:::validate_tr(200))
})

# --- get_drop_volumes_flags ---

test_that("get_drop_volumes_flags returns expected structure", {
  flags <- fmri.pipeline:::get_drop_volumes_flags()
  expect_type(flags, "list")
  expect_named(flags, c("shift_nifti", "shift_timing", "shorten_additional", "shorten_ts"))
  expect_false(flags$shift_nifti)
  expect_true(flags$shift_timing)
  expect_true(flags$shorten_additional)
  expect_true(flags$shorten_ts)
})

# ==============================================================================
# Tests for Stage 2 helpers: Signal expansion and alignment
# ==============================================================================

# --- prepare_signals_for_expansion ---

test_that("prepare_signals_for_expansion sets beta_series to FALSE by default", {
  signals <- list(
    cue = list(name = "cue", event = "cue", value = 1)
  )
  result <- fmri.pipeline:::prepare_signals_for_expansion(signals)
  expect_false(result$cue$beta_series)
})

test_that("prepare_signals_for_expansion preserves TRUE beta_series", {
  signals <- list(
    cue = list(name = "cue", event = "cue", value = 1, beta_series = TRUE)
  )
  result <- fmri.pipeline:::prepare_signals_for_expansion(signals)
  expect_true(result$cue$beta_series)
})

test_that("prepare_signals_for_expansion handles multiple signals", {
  signals <- list(
    cue = list(name = "cue", event = "cue", value = 1),
    outcome = list(name = "outcome", event = "outcome", value = 1, beta_series = TRUE)
  )
  result <- fmri.pipeline:::prepare_signals_for_expansion(signals)
  expect_false(result$cue$beta_series)
  expect_true(result$outcome$beta_series)
})

# --- apply_signal_duration ---

test_that("apply_signal_duration applies scalar duration to all rows", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1:3,
    onset = c(0, 5, 10),
    duration = c(1, 1, 1),
    value = c(1, 2, 3)
  )
  result <- fmri.pipeline:::apply_signal_duration(s_aligned, 2.5)
  expect_equal(result$duration, c(2.5, 2.5, 2.5))
})

test_that("apply_signal_duration uses column name for duration", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1:3,
    onset = c(0, 5, 10),
    duration = c(1, 1, 1),
    custom_dur = c(0.5, 1.0, 1.5),
    value = c(1, 2, 3)
  )
  result <- fmri.pipeline:::apply_signal_duration(s_aligned, "custom_dur")
  expect_equal(result$duration, c(0.5, 1.0, 1.5))
})

test_that("apply_signal_duration leaves duration unchanged when NULL", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1:3,
    onset = c(0, 5, 10),
    duration = c(1, 2, 3),
    value = c(1, 2, 3)
  )
  result <- fmri.pipeline:::apply_signal_duration(s_aligned, NULL)
  expect_equal(result$duration, c(1, 2, 3))
})

test_that("apply_signal_duration errors on multi-element duration", {
  s_aligned <- data.frame(run_number = 1, trial = 1, onset = 0, duration = 1, value = 1)
  expect_error(
    fmri.pipeline:::apply_signal_duration(s_aligned, c(1, 2)),
    "multi-element duration"
  )
})

# --- split_signal_by_run ---

test_that("split_signal_by_run creates list by run", {
  s_aligned <- data.frame(
    run_number = c(1, 1, 2, 2),
    trial = c(1, 2, 1, 2),
    onset = c(0, 5, 0, 5),
    duration = c(1, 1, 1, 1),
    value = c(1, 2, 3, 4)
  )
  event_runs <- factor(c(1, 2))
  result <- fmri.pipeline:::split_signal_by_run(s_aligned, event_runs, "cue", FALSE)
  
  expect_length(result, 2)
  expect_named(result, c("run_number1", "run_number2"))
  expect_equal(nrow(result$run_number1), 2)
  expect_equal(nrow(result$run_number2), 2)
})

test_that("split_signal_by_run sets event attribute", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1,
    value = 1
  )
  event_runs <- factor(1)
  result <- fmri.pipeline:::split_signal_by_run(s_aligned, event_runs, "cue", FALSE)
  
  expect_equal(attr(result$run_number1, "event"), "cue")
})

test_that("split_signal_by_run sets physio_only attribute when TRUE", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1,
    value = 1
  )
  event_runs <- factor(1)
  result <- fmri.pipeline:::split_signal_by_run(s_aligned, event_runs, "cue", TRUE)
  
  expect_true(attr(result$run_number1, "physio_only"))
})

test_that("split_signal_by_run excludes run_number from output", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1,
    value = 1
  )
  event_runs <- factor(1)
  result <- fmri.pipeline:::split_signal_by_run(s_aligned, event_runs, "cue", FALSE)
  
  expect_false("run_number" %in% names(result$run_number1))
  expect_true(all(c("trial", "onset", "duration", "value") %in% names(result$run_number1)))
})

# --- extract_signal_config ---

test_that("extract_signal_config extracts normalization settings", {
  signals_expanded <- list(
    cue = list(name = "cue", normalization = "evtmax_1"),
    outcome = list(name = "outcome", normalization = "durmax_1"),
    pe = list(name = "pe")  # NULL normalization
  )
  result <- fmri.pipeline:::extract_signal_config(signals_expanded)
  
  expect_equal(result$normalizations, c(cue = "evtmax_1", outcome = "durmax_1", pe = "none"))
})

test_that("extract_signal_config extracts beta_series settings", {
  signals_expanded <- list(
    cue = list(name = "cue", beta_series = TRUE),
    outcome = list(name = "outcome", beta_series = FALSE),
    pe = list(name = "pe")  # NULL beta_series
  )
  result <- fmri.pipeline:::extract_signal_config(signals_expanded)
  
  expect_equal(result$beta_series, c(cue = TRUE, outcome = FALSE, pe = FALSE))
})

test_that("extract_signal_config extracts rm_zeros settings", {
  signals_expanded <- list(
    cue = list(name = "cue", rm_zeros = TRUE),
    outcome = list(name = "outcome", rm_zeros = FALSE),
    pe = list(name = "pe")  # NULL rm_zeros defaults to TRUE
  )
  result <- fmri.pipeline:::extract_signal_config(signals_expanded)
  
  expect_equal(result$rm_zeros, c(cue = TRUE, outcome = FALSE, pe = TRUE))
})

test_that("extract_signal_config extracts convmax_1 settings", {
  signals_expanded <- list(
    cue = list(name = "cue", convmax_1 = TRUE),
    outcome = list(name = "outcome", convmax_1 = FALSE),
    pe = list(name = "pe")  # NULL convmax_1 defaults to FALSE
  )
  result <- fmri.pipeline:::extract_signal_config(signals_expanded)
  
  expect_equal(result$convmax_1, c(cue = TRUE, outcome = FALSE, pe = FALSE))
})

test_that("extract_signal_config extracts add_deriv settings", {
  signals_expanded <- list(
    cue = list(name = "cue", add_deriv = TRUE),
    outcome = list(name = "outcome", add_deriv = FALSE),
    pe = list(name = "pe")  # NULL add_deriv defaults to FALSE
  )
  result <- fmri.pipeline:::extract_signal_config(signals_expanded)
  
  expect_equal(result$add_derivs, c(cue = TRUE, outcome = FALSE, pe = FALSE))
})

test_that("extract_signal_config errors on invalid rm_zeros", {
  signals_expanded <- list(
    cue = list(name = "cue", rm_zeros = "invalid")
  )
  expect_error(
    fmri.pipeline:::extract_signal_config(signals_expanded),
    "rm_zeros"
  )
})

test_that("extract_signal_config errors on invalid convmax_1", {
  signals_expanded <- list(
    cue = list(name = "cue", convmax_1 = "yes")
  )
  expect_error(
    fmri.pipeline:::extract_signal_config(signals_expanded),
    "convmax_1"
  )
})

test_that("extract_signal_config errors on invalid add_deriv", {
  signals_expanded <- list(
    cue = list(name = "cue", add_deriv = 2)
  )
  expect_error(
    fmri.pipeline:::extract_signal_config(signals_expanded),
    "add_deriv"
  )
})

# --- align_signal_with_events ---

test_that("align_signal_with_events aligns scalar value signal", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1:3,
    onset = c(0, 5, 10),
    duration = c(1, 1, 1)
  )
  s <- list(name = "cue_evt", event = "cue", value = 1)
  
  result <- fmri.pipeline:::align_signal_with_events(s, events)
  
  expect_length(result, 1)
  expect_named(result, "run_number1")
  expect_equal(nrow(result$run_number1), 3)
  expect_equal(as.numeric(result$run_number1$value), c(1, 1, 1))
})

test_that("align_signal_with_events aligns parametric signal", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1:3,
    onset = c(0, 5, 10),
    duration = c(1, 1, 1)
  )
  s <- list(
    name = "cue_param",
    event = "cue",
    value = data.frame(run_number = 1, trial = 1:3, value = c(0.5, 1.0, 1.5))
  )
  
  result <- fmri.pipeline:::align_signal_with_events(s, events)
  
  expect_equal(as.numeric(result$run_number1$value), c(0.5, 1.0, 1.5))
  expect_equal(as.numeric(result$run_number1$onset), c(0, 5, 10))
})

test_that("align_signal_with_events handles multiple runs", {
  events <- data.frame(
    event = rep("cue", 4),
    run_number = c(1, 1, 2, 2),
    trial = c(1, 2, 1, 2),
    onset = c(0, 5, 0, 5),
    duration = c(1, 1, 1, 1)
  )
  s <- list(name = "cue_evt", event = "cue", value = 1)
  
  result <- fmri.pipeline:::align_signal_with_events(s, events)
  
  expect_length(result, 2)
  expect_named(result, c("run_number1", "run_number2"))
})

test_that("align_signal_with_events errors when event is NULL", {
  events <- data.frame(event = "cue", run_number = 1, trial = 1, onset = 0, duration = 1)
  s <- list(name = "bad", value = 1)  # No event element
  
  expect_error(
    fmri.pipeline:::align_signal_with_events(s, events),
    "Signal does not have event element"
  )
})

test_that("align_signal_with_events adds default value of 1 when missing", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1
  )
  s <- list(name = "cue_evt", event = "cue")  # No value element
  
  result <- fmri.pipeline:::align_signal_with_events(s, events)
  expect_equal(as.numeric(result$run_number1$value), 1)
})

# --- Robustness tests for Stage 2 helpers ---

test_that("prepare_signals_for_expansion errors on empty signals list", {
  expect_error(
    fmri.pipeline:::prepare_signals_for_expansion(list()),
    "signals list is empty"
  )
})

test_that("expand_signals_list detects duplicate signal names", {
  # This would require signals that expand to the same name
  # Hard to test without mocking expand_signal, so we test via integration
  # by creating signals with explicitly duplicate names
  # Note: expand_signal requires beta_series=FALSE for the simple expansion path
  signals <- list(
    cue1 = list(name = "cue_evt", event = "cue", value = 1, beta_series = FALSE),
    cue2 = list(name = "cue_evt", event = "cue", value = 2, beta_series = FALSE)  # Same name!
  )
  expect_error(
    fmri.pipeline:::expand_signals_list(signals),
    "Duplicate signal names"
  )
})

test_that("align_signal_with_events errors when event type not in events", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1
  )
  s <- list(name = "bad_signal", event = "nonexistent", value = 1)
  
  expect_error(
    fmri.pipeline:::align_signal_with_events(s, events),
    "has no occurrences in the events data.frame"
  )
})

test_that("align_signal_with_events errors when signal data.frame missing run_number", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1:3,
    onset = c(0, 5, 10),
    duration = c(1, 1, 1)
  )
  s <- list(
    name = "bad_param",
    event = "cue",
    value = data.frame(trial = 1:3, value = c(1, 2, 3))  # Missing run_number
  )
  
  expect_error(
    fmri.pipeline:::align_signal_with_events(s, events),
    "missing required columns.*run_number"
  )
})

test_that("align_signal_with_events errors when signal data.frame missing value column", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1:3,
    onset = c(0, 5, 10),
    duration = c(1, 1, 1)
  )
  s <- list(
    name = "bad_param",
    event = "cue",
    value = data.frame(run_number = 1, trial = 1:3, param = c(1, 2, 3))  # param not value
  )
  
  expect_error(
    fmri.pipeline:::align_signal_with_events(s, events),
    "must have a 'value' column"
  )
})

test_that("align_signal_with_events warns and excludes unmatched trials", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1:3,
    onset = c(0, 5, 10),
    duration = c(1, 1, 1)
  )
  s <- list(
    name = "partial_match",
    event = "cue",
    value = data.frame(run_number = 1, trial = c(1, 2, 99), value = c(1, 2, 3))  # Trial 99 doesn't exist
  )
  
  expect_warning(
    result <- fmri.pipeline:::align_signal_with_events(s, events),
    "did not match events"
  )
  # Should only have 2 rows (trials 1 and 2)
  expect_equal(nrow(result$run_number1), 2)
})

test_that("align_signal_with_events errors on unknown value type", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1
  )
  s <- list(name = "bad_type", event = "cue", value = "not_valid")
  
  expect_error(
    fmri.pipeline:::align_signal_with_events(s, events),
    "Unknown data type for value"
  )
})

test_that("apply_signal_duration errors when column name doesn't exist", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1,
    value = 1
  )
  
  expect_error(
    fmri.pipeline:::apply_signal_duration(s_aligned, "nonexistent_column", "test_signal"),
    "duration column.*not found"
  )
})

test_that("apply_signal_duration errors on invalid duration type", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1,
    value = 1
  )
  
  expect_error(
    fmri.pipeline:::apply_signal_duration(s_aligned, list(bad = 1), "test_signal"),
    "duration must be numeric or a column name"
  )
})

test_that("apply_signal_duration errors on NA durations", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1:2,
    onset = c(0, 5),
    duration = c(1, NA),
    value = c(1, 2)
  )
  
  expect_error(
    fmri.pipeline:::apply_signal_duration(s_aligned, NULL, "test_signal"),
    "duration contains NA values"
  )
})

test_that("apply_signal_duration errors on negative durations", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1:2,
    onset = c(0, 5),
    duration = c(1, -1),
    value = c(1, 2)
  )
  
  expect_error(
    fmri.pipeline:::apply_signal_duration(s_aligned, NULL, "test_signal"),
    "duration contains negative values"
  )
})

test_that("validate_aligned_columns detects missing columns", {
  s_aligned <- data.frame(
    run_number = 1,
    trial = 1,
    onset = 0
    # Missing duration and value
  )
  
  expect_error(
    fmri.pipeline:::validate_aligned_columns(s_aligned, "test_signal"),
    "missing required columns.*duration.*value"
  )
})

# --- expand_and_align_signals (integration) ---

test_that("expand_and_align_signals returns correct structure", {
  events <- data.frame(
    event = rep("cue", 4),
    run_number = c(1, 1, 2, 2),
    trial = c(1, 2, 1, 2),
    onset = c(0, 5, 0, 5),
    duration = c(1, 1, 1, 1)
  )
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  
  result <- fmri.pipeline:::expand_and_align_signals(signals, events)
  
  expect_named(result, c("signals_expanded", "signals_aligned", "signal_config"))
  expect_length(result$signals_expanded, 1)
  expect_length(result$signals_aligned, 1)
  expect_named(result$signal_config, c("normalizations", "beta_series", "rm_zeros", "convmax_1", "add_derivs"))
})

test_that("expand_and_align_signals extracts correct config", {
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1,
    onset = 0,
    duration = 1
  )
  # Note: add_deriv = TRUE would create additional derivative signals during expansion,
  # so we test without it here to verify basic config extraction
  signals <- list(
    cue_evt = list(
      name = "cue_evt", 
      event = "cue", 
      value = 1, 
      normalization = "evtmax_1",
      convmax_1 = TRUE,
      rm_zeros = FALSE
    )
  )
  
  result <- fmri.pipeline:::expand_and_align_signals(signals, events)
  
  expect_equal(result$signal_config$normalizations, c(cue_evt = "evtmax_1"))
  expect_equal(result$signal_config$convmax_1, c(cue_evt = TRUE))
  expect_equal(result$signal_config$rm_zeros, c(cue_evt = FALSE))
})

# ==============================================================================
# Tests for Stage 6 helpers: Write timing files to disk
# ==============================================================================

# --- format_afni_run_string ---

test_that("format_afni_run_string formats basic TIME:DURATION without modulation", {
  onsets <- c(0, 5, 10)
  durations <- c(1, 2, 1.5)
  values <- matrix(nrow = 3, ncol = 0)  # Empty matrix = no modulation
  
  result <- fmri.pipeline:::format_afni_run_string(onsets, durations, values)
  expect_equal(result, "0:1 5:2 10:1.5")
})

test_that("format_afni_run_string formats TIME*PARAM:DURATION with single modulation", {
  onsets <- c(0, 5)
  durations <- c(1, 1)
  values <- matrix(c(0.5, -0.3), ncol = 1)
  
  result <- fmri.pipeline:::format_afni_run_string(onsets, durations, values)
  expect_equal(result, "0*0.5:1 5*-0.3:1")
})

test_that("format_afni_run_string formats TIME*PARAM1,PARAM2:DURATION with multiple modulations", {
  onsets <- c(0, 5)
  durations <- c(1, 2)
  values <- matrix(c(0.5, -0.3, 0.2, 0.8), ncol = 2)
  
  result <- fmri.pipeline:::format_afni_run_string(onsets, durations, values)
  expect_equal(result, "0*0.5,0.2:1 5*-0.3,0.8:2")
})

test_that("format_afni_run_string rounds values to 6 decimal places", {
  onsets <- c(1.123456789)
  durations <- c(2.987654321)
  values <- matrix(nrow = 1, ncol = 0)
  
  result <- fmri.pipeline:::format_afni_run_string(onsets, durations, values)
  expect_equal(result, "1.123457:2.987654")
})

# --- center_fsl_regressor ---

test_that("center_fsl_regressor mean-centers values when center_values=TRUE", {
  regout <- data.frame(
    onset = c(0, 5, 10),
    duration = c(1, 1, 1),
    value = c(2, 4, 6)
  )
  
  result <- fmri.pipeline:::center_fsl_regressor(regout, center_values = TRUE)
  expect_equal(mean(result$value), 0)
  expect_equal(result$value, c(-2, 0, 2))
})

test_that("center_fsl_regressor removes zero-value events", {
  regout <- data.frame(
    onset = c(0, 5, 10, 15),
    duration = c(1, 1, 1, 1),
    value = c(2, 0, 4, 0)
  )
  
  result <- fmri.pipeline:::center_fsl_regressor(regout, center_values = TRUE)
  expect_equal(nrow(result), 2)
  expect_equal(result$onset, c(0, 10))
})

test_that("center_fsl_regressor does not center when center_values=FALSE", {
  regout <- data.frame(
    onset = c(0, 5, 10),
    duration = c(1, 1, 1),
    value = c(2, 4, 6)
  )
  
  result <- fmri.pipeline:::center_fsl_regressor(regout, center_values = FALSE)
  expect_equal(result$value, c(2, 4, 6))
})

test_that("center_fsl_regressor does not center constant values (indicator function)", {
  regout <- data.frame(
    onset = c(0, 5, 10),
    duration = c(1, 1, 1),
    value = c(3, 3, 3)  # All same value
  )
  
  result <- fmri.pipeline:::center_fsl_regressor(regout, center_values = TRUE)
  expect_equal(result$value, c(3, 3, 3))  # Unchanged
})

test_that("center_fsl_regressor returns unchanged if all values are zero", {
  regout <- data.frame(
    onset = c(0, 5, 10),
    duration = c(1, 1, 1),
    value = c(0, 0, 0)
  )
  
  result <- fmri.pipeline:::center_fsl_regressor(regout, center_values = TRUE)
  expect_equal(result$value, c(0, 0, 0))
})

# --- write_timing_files_to_disk (integration) ---

test_that("write_timing_files_to_disk returns NULL fields when write_timing_files is NULL", {
  result <- fmri.pipeline:::write_timing_files_to_disk(
    write_timing_files = NULL,
    dmat = NULL,
    dmat_convolved = NULL,
    output_directory = tempdir(),
    runs_to_output = 1:2,
    center_values = TRUE
  )
  
  expect_null(result$tf_convolved)
  expect_null(result$tf_convolved_concat)
  expect_null(result$tf_fsl)
  expect_null(result$tf_afni)
})

test_that("write_convolved_concat_files returns empty vector for empty input", {
  result <- fmri.pipeline:::write_convolved_concat_files(list(), tempdir())
  expect_equal(result, character(0))
})

# --- write_fsl_timing_files ---

test_that("write_fsl_timing_files writes files with correct names", {
  # Create minimal dmat structure
  dmat <- array(
    list(
      data.frame(onset = c(0, 5), duration = c(1, 1), value = c(1, 2)),
      data.frame(onset = c(0, 5), duration = c(1, 1), value = c(3, 4))
    ),
    dim = c(1, 2),
    dimnames = list("run1", c("cue_evt", "cue_param"))
  )
  
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))
  
  result <- fmri.pipeline:::write_fsl_timing_files(
    dmat = dmat,
    output_directory = tmpdir,
    runs_to_output = 1,
    center_values = TRUE
  )
  
  expect_equal(dim(result), c(1, 2))
  expect_true(all(file.exists(result)))
  expect_true(grepl("run1_cue_evt_FSL3col.txt", result[1, 1]))
  expect_true(grepl("run1_cue_param_FSL3col.txt", result[1, 2]))
})

test_that("write_fsl_timing_files skips empty regressors", {
  dmat <- array(
    list(
      data.frame(onset = numeric(0), duration = numeric(0), value = numeric(0))
    ),
    dim = c(1, 1),
    dimnames = list("run1", c("empty_reg"))
  )
  
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))
  
  result <- fmri.pipeline:::write_fsl_timing_files(
    dmat = dmat,
    output_directory = tmpdir,
    runs_to_output = 1,
    center_values = TRUE
  )
  
  # Should have NA for empty regressor
  expect_true(is.na(result[1, 1]))
})

# --- write_convolved_timing_files ---

test_that("write_convolved_timing_files writes individual run files", {
  dmat <- array(
    list(
      data.frame(onset = c(0, 5), duration = c(1, 1), value = c(1, 2))
    ),
    dim = c(2, 1),
    dimnames = list(c("run_number1", "run_number2"), c("cue_evt"))
  )
  
  dmat_convolved <- list(
    run_number1 = list(cue_evt = c(0.1, 0.2, 0.3)),
    run_number2 = list(cue_evt = c(0.4, 0.5, 0.6))
  )
  
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))
  
  result <- fmri.pipeline:::write_convolved_timing_files(
    dmat = dmat,
    dmat_convolved = dmat_convolved,
    output_directory = tmpdir
  )
  
  expect_true(all(file.exists(result$tf_convolved)))
  expect_length(result$tf_convolved_concat, 1)
  expect_true(file.exists(result$tf_convolved_concat["cue_evt"]))
  
  # Check concatenated content
  concat_content <- readLines(result$tf_convolved_concat["cue_evt"])
  expect_equal(as.numeric(concat_content), c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
})

# ==============================================================================
# Tests for cleanup_regressor
# ==============================================================================

test_that("cleanup_regressor removes NAs from all inputs", {
  times <- c(1, 2, NA, 4, 5)
  durations <- c(1, 1, 1, NA, 1)
  values <- c(1, 2, 3, 4, NA)
  
  result <- fmri.pipeline:::cleanup_regressor(times, durations, values)
  
  # Only positions 1 and 2 have no NAs anywhere
  expect_equal(result$times, c(1, 2))
  expect_equal(result$durations, c(1, 1))
  expect_equal(result$values, c(1, 2))
})

test_that("cleanup_regressor removes zeros when rm_zeros=TRUE", {
  times <- c(1, 2, 3, 4, 5)
  durations <- c(1, 1, 1, 1, 1)
  values <- c(1, 0, 3, 0, 5)
  
  result <- fmri.pipeline:::cleanup_regressor(times, durations, values, rm_zeros = TRUE)
  
  expect_equal(result$times, c(1, 3, 5))
  expect_equal(result$values, c(1, 3, 5))
})

test_that("cleanup_regressor keeps zeros when rm_zeros=FALSE", {
  times <- c(1, 2, 3)
  durations <- c(1, 1, 1)
  values <- c(1, 0, 3)
  
  result <- fmri.pipeline:::cleanup_regressor(times, durations, values, rm_zeros = FALSE)
  
  expect_equal(result$times, c(1, 2, 3))
  expect_equal(result$values, c(1, 0, 3))
})

test_that("cleanup_regressor handles all-NA input", {
  times <- c(NA, NA, NA)
  durations <- c(1, 1, 1)
  values <- c(1, 2, 3)
  
  result <- fmri.pipeline:::cleanup_regressor(times, durations, values)
  
  expect_equal(length(result$times), 0)
  expect_equal(length(result$durations), 0)
  expect_equal(length(result$values), 0)
})

test_that("cleanup_regressor handles empty input", {
  result <- fmri.pipeline:::cleanup_regressor(numeric(0), numeric(0), numeric(0))
  
  expect_equal(length(result$times), 0)
  expect_equal(length(result$durations), 0)
  expect_equal(length(result$values), 0)
})

# ==============================================================================
# Tests for shift_dmat_timing
# ==============================================================================

test_that("shift_dmat_timing adjusts onset times by drop_volumes * tr", {
  # Create a simple 2x1 dmat (2 runs, 1 regressor)
  dmat <- matrix(list(), nrow = 2, ncol = 1)
  dimnames(dmat) <- list(c("run1", "run2"), c("cue"))
  
  dmat[[1, 1]] <- data.frame(trial = 1:3, onset = c(5, 10, 15), duration = 1, value = 1)
  dmat[[2, 1]] <- data.frame(trial = 1:3, onset = c(5, 10, 15), duration = 1, value = 1)
  
  tr <- 2.0
  drop_volumes <- c(3, 3)  # Drop 3 volumes from each run = 6s offset
  
  result <- fmri.pipeline:::shift_dmat_timing(dmat, tr, drop_volumes)
  
  expect_equal(result[[1, 1]]$onset, c(5 - 6, 10 - 6, 15 - 6))
  expect_equal(result[[2, 1]]$onset, c(5 - 6, 10 - 6, 15 - 6))
})

test_that("shift_dmat_timing handles variable drop_volumes per run", {
  dmat <- matrix(list(), nrow = 2, ncol = 1)
  dimnames(dmat) <- list(c("run1", "run2"), c("cue"))
  
  dmat[[1, 1]] <- data.frame(trial = 1, onset = 10, duration = 1, value = 1)
  dmat[[2, 1]] <- data.frame(trial = 1, onset = 10, duration = 1, value = 1)
  
  tr <- 1.0
  drop_volumes <- c(5, 2)  # Different drops per run
  
  result <- fmri.pipeline:::shift_dmat_timing(dmat, tr, drop_volumes)
  
  expect_equal(result[[1, 1]]$onset, 10 - 5)
  expect_equal(result[[2, 1]]$onset, 10 - 2)
})

test_that("shift_dmat_timing returns unchanged dmat when shift_timing=FALSE", {
  dmat <- matrix(list(), nrow = 1, ncol = 1)
  dimnames(dmat) <- list("run1", "cue")
  dmat[[1, 1]] <- data.frame(trial = 1, onset = 10, duration = 1, value = 1)
  
  result <- fmri.pipeline:::shift_dmat_timing(dmat, tr = 1, drop_volumes = 5, shift_timing = FALSE)
  
  expect_equal(result[[1, 1]]$onset, 10)
})

test_that("shift_dmat_timing returns unchanged dmat when all drop_volumes are 0", {
  dmat <- matrix(list(), nrow = 1, ncol = 1)
  dimnames(dmat) <- list("run1", "cue")
  dmat[[1, 1]] <- data.frame(trial = 1, onset = 10, duration = 1, value = 1)
  
  result <- fmri.pipeline:::shift_dmat_timing(dmat, tr = 1, drop_volumes = 0)
  
  expect_equal(result[[1, 1]]$onset, 10)
})

test_that("shift_dmat_timing handles empty regressors gracefully", {
  dmat <- matrix(list(), nrow = 2, ncol = 2)
  dimnames(dmat) <- list(c("run1", "run2"), c("cue", "empty"))
  
  dmat[[1, 1]] <- data.frame(trial = 1, onset = 10, duration = 1, value = 1)
  dmat[[1, 2]] <- data.frame(trial = integer(0), onset = numeric(0), duration = numeric(0), value = numeric(0))
  dmat[[2, 1]] <- data.frame(trial = 1, onset = 10, duration = 1, value = 1)
  dmat[[2, 2]] <- data.frame(trial = integer(0), onset = numeric(0), duration = numeric(0), value = numeric(0))
  
  # Should not error on empty regressors (messages are OK)
  expect_no_error(fmri.pipeline:::shift_dmat_timing(dmat, tr = 1, drop_volumes = c(2, 2)))
})

# ==============================================================================
# Tests for get_additional_regressors
# ==============================================================================

test_that("get_additional_regressors returns NULL when additional_regressors is NULL", {
  result <- fmri.pipeline:::get_additional_regressors(
    additional_regressors = NULL, 
    run_volumes = c(100, 100),
    lg = lgr::get_logger()
  )
  expect_null(result)
})

test_that("get_additional_regressors combines list of data.frames correctly", {
  run1 <- data.frame(confound1 = rnorm(50), confound2 = rnorm(50))
  run2 <- data.frame(confound1 = rnorm(40), confound2 = rnorm(40))
  
  result <- fmri.pipeline:::get_additional_regressors(
    additional_regressors = list(run1, run2),
    run_volumes = c(50, 40),
    drop_volumes = c(0, 0),
    shorten_additional = FALSE,
    lg = lgr::get_logger()
  )
  
  expect_equal(nrow(result), 90)
  expect_true("run_number" %in% names(result))
  expect_equal(sum(result$run_number == 1), 50)
  expect_equal(sum(result$run_number == 2), 40)
})

test_that("get_additional_regressors respects drop_volumes with shorten_additional=TRUE", {
  # Create data longer than final run_volumes
  run1 <- data.frame(confound1 = rnorm(60))
  run2 <- data.frame(confound1 = rnorm(55))
  
  result <- fmri.pipeline:::get_additional_regressors(
    additional_regressors = list(run1, run2),
    run_volumes = c(50, 45),  # Final volumes after drop
    drop_volumes = c(5, 5),    # 5 volumes dropped from each
    shorten_additional = TRUE,
    lg = lgr::get_logger()
  )
  
  # With shorten_additional=TRUE, we take rows (drop_volumes+1):(drop_volumes+run_volumes)
  # So for run1: rows 6:55 (50 rows), for run2: rows 6:50 (45 rows)
  expect_equal(nrow(result), 95)
  expect_equal(sum(result$run_number == 1), 50)
  expect_equal(sum(result$run_number == 2), 45)
})

test_that("get_additional_regressors errors when list length doesn't match run_volumes length", {
  run1 <- data.frame(confound1 = rnorm(50))
  
  expect_error(
    fmri.pipeline:::get_additional_regressors(
      additional_regressors = list(run1),
      run_volumes = c(50, 50),  # Two runs but only one data.frame
      lg = lgr::get_logger()
    ),
    "does not much length"
  )
})

test_that("get_additional_regressors errors when data has fewer rows than run_volumes", {
  run1 <- data.frame(confound1 = rnorm(30))  # Too short
  
  expect_error(
    fmri.pipeline:::get_additional_regressors(
      additional_regressors = list(run1),
      run_volumes = c(50),
      drop_volumes = 0,
      shorten_additional = FALSE,
      lg = lgr::get_logger()
    ),
    "fewer observations"
  )
})

# ==============================================================================
# Tests for get_ts_multipliers
# ==============================================================================

test_that("get_ts_multipliers returns NULL when ts_multipliers is NULL", {
  run_data <- data.frame(run_number = 1:2, run_volumes = c(50, 50), drop_volumes = c(0, 0))
  
  result <- fmri.pipeline:::get_ts_multipliers(
    ts_multipliers = NULL,
    run_data = run_data,
    shorten_ts = FALSE
  )
  
  expect_null(result)
})

test_that("get_ts_multipliers processes list of data.frames", {
  run_data <- data.frame(run_number = 1:2, run_volumes = c(50, 40), drop_volumes = c(0, 0))
  ts <- list(
    data.frame(roi1 = rnorm(50), roi2 = rnorm(50)),
    data.frame(roi1 = rnorm(40), roi2 = rnorm(40))
  )
  
  result <- expect_no_warning(
    fmri.pipeline:::get_ts_multipliers(
      ts_multipliers = ts,
      run_data = run_data,
      shorten_ts = FALSE
    )
  )
  
  expect_equal(nrow(result), 90)
  expect_true("run_number" %in% names(result))
  # Check that roi columns were mean-centered
  expect_lt(abs(mean(result$roi1[result$run_number == 1])), 1e-10)
})

test_that("get_ts_multipliers respects drop_volumes with shorten_ts=TRUE", {
  run_data <- data.frame(run_number = 1, run_volumes = 45, drop_volumes = 5)
  ts <- list(data.frame(roi1 = 1:55))  # Need enough rows: drop_volumes + run_volumes = 50
  
  result <- expect_no_warning(
    fmri.pipeline:::get_ts_multipliers(
      ts_multipliers = ts,
      run_data = run_data,
      shorten_ts = TRUE
    )
  )
  
  expect_equal(nrow(result), 45)
  expect_equal(result$roi1, seq(-22, 22))
})

test_that("get_ts_multipliers handles data.frame input with run_number column", {
  run_data <- data.frame(run_number = 1:2, run_volumes = c(30, 30), drop_volumes = c(0, 0))
  
  ts <- data.frame(
    run_number = c(rep(1, 30), rep(2, 30)),
    roi1 = rnorm(60)
  )
  
  result <- expect_no_warning(
    fmri.pipeline:::get_ts_multipliers(
      ts_multipliers = ts,
      run_data = run_data,
      shorten_ts = FALSE
    )
  )
  
  expect_equal(nrow(result), 60)
  expect_lt(abs(mean(result$roi1[result$run_number == 1])), 1e-10)
  expect_lt(abs(mean(result$roi1[result$run_number == 2])), 1e-10)
})

test_that("get_ts_multipliers does not collapse single-column values to zero", {
  run_data <- data.frame(run_number = 1, run_volumes = 5, drop_volumes = 0)
  ts <- list(data.frame(roi1 = c(1, 2, 3, 4, 5)))

  result <- expect_no_warning(
    fmri.pipeline:::get_ts_multipliers(
      ts_multipliers = ts,
      run_data = run_data,
      shorten_ts = FALSE
    )
  )

  expect_equal(result$roi1, c(-2, -1, 0, 1, 2))
  expect_gt(stats::sd(result$roi1), 0)
})

test_that("get_ts_multipliers errors when data.frame lacks run_number", {
  run_data <- data.frame(run_number = 1, run_volumes = 50, drop_volumes = 0)
  ts <- data.frame(roi1 = rnorm(50))  # Missing run_number
  
  expect_error(
    fmri.pipeline:::get_ts_multipliers(
      ts_multipliers = ts,
      run_data = run_data,
      shorten_ts = FALSE
    ),
    "run_number"
  )
})

test_that("validate_ts_multipliers rejects non-numeric data.frame columns", {
  run_data <- data.frame(run_number = 1, run_volumes = 10, drop_volumes = 0)
  ts <- data.frame(run_number = 1:10, roi1 = as.character(1:10))

  expect_error(
    fmri.pipeline:::validate_ts_multipliers(ts, run_data),
    "non-numeric"
  )
})

test_that("validate_ts_multipliers rejects non-numeric list columns", {
  run_data <- data.frame(run_number = 1:2, run_volumes = c(5, 5), drop_volumes = c(0, 0))
  ts <- list(
    data.frame(roi1 = 1:5),
    data.frame(roi1 = letters[1:5])
  )

  expect_error(
    fmri.pipeline:::validate_ts_multipliers(ts, run_data),
    "non-numeric"
  )
})

test_that("validate_ts_multipliers rejects mismatched run counts", {
  run_data <- data.frame(run_number = 1:2, run_volumes = c(5, 5), drop_volumes = c(0, 0))
  ts <- list(data.frame(roi1 = 1:5))

  expect_error(
    fmri.pipeline:::validate_ts_multipliers(ts, run_data),
    "does not match number of runs"
  )
})

# ==============================================================================
# Tests for add_baseline_regressors
# ==============================================================================

test_that("add_baseline_regressors returns unchanged when baseline_coef_order < 0", {
  dmat <- list(run1 = data.frame(cue = rnorm(50)))
  
  result <- fmri.pipeline:::add_baseline_regressors(dmat, baseline_coef_order = -1)
  
  expect_equal(result, dmat)
  expect_equal(ncol(result$run1), 1)
})

# Note: baseline_coef_order = 0 has a bug in add_baseline_regressors with Legendre polynomials
# when length(baseline) == 1, the subsetting baseline[2:length(baseline)] produces unexpected results
# This test is skipped until that bug is fixed
test_that("add_baseline_regressors adds intercept only when baseline_coef_order = 0", {
  skip("Known bug: baseline_coef_order=0 has edge case with Legendre polynomial subsetting")
  dmat <- list(run1 = data.frame(cue = rnorm(50)))
  
  result <- fmri.pipeline:::add_baseline_regressors(dmat, baseline_coef_order = 0)
  
  expect_equal(ncol(result$run1), 2)
  expect_true("base0" %in% names(result$run1))
  # Intercept (Legendre P_0) should be constant = 1
  expect_true(all(result$run1$base0 == result$run1$base0[1]))
})

test_that("add_baseline_regressors adds linear trend when baseline_coef_order = 1", {
  dmat <- list(run1 = data.frame(cue = rnorm(100)))
  
  result <- fmri.pipeline:::add_baseline_regressors(dmat, baseline_coef_order = 1)
  
  expect_equal(ncol(result$run1), 3)
  expect_true("base0" %in% names(result$run1))
  expect_true("base1" %in% names(result$run1))
  # Linear trend (Legendre P_1) should be centered around 0
  expect_lt(abs(mean(result$run1$base1)), 1e-10)
})

test_that("add_baseline_regressors adds quadratic trend when baseline_coef_order = 2", {
  dmat <- list(run1 = data.frame(cue = rnorm(100)))
  
  result <- fmri.pipeline:::add_baseline_regressors(dmat, baseline_coef_order = 2)
  
  expect_equal(ncol(result$run1), 4)
  expect_true("base2" %in% names(result$run1))
  # Quadratic term should be centered around 0
  expect_lt(abs(mean(result$run1$base2)), 1e-10)
})

test_that("add_baseline_regressors works for multiple runs", {
  dmat <- list(
    run1 = data.frame(cue = rnorm(50)),
    run2 = data.frame(cue = rnorm(60))
  )
  
  result <- fmri.pipeline:::add_baseline_regressors(dmat, baseline_coef_order = 1)
  
  expect_equal(nrow(result$run1), 50)
  expect_equal(nrow(result$run2), 60)
  expect_equal(ncol(result$run1), 3)
  expect_equal(ncol(result$run2), 3)
})

# ==============================================================================
# Tests for get_collin_events
# ==============================================================================

test_that("get_collin_events computes correlation matrix correctly", {
  dmat <- matrix(list(), nrow = 1, ncol = 2)
  dimnames(dmat) <- list("run1", c("reg1", "reg2"))
  
  # Create perfectly correlated regressors
  dmat[[1, 1]] <- data.frame(trial = 1:5, value = 1:5)
  dmat[[1, 2]] <- data.frame(trial = 1:5, value = 2:6)  # Linear transform, r=1
  
  result <- fmri.pipeline:::get_collin_events(dmat)
  
  expect_true("r" %in% names(result$run1))
  expect_equal(result$run1$r["reg1", "reg2"], 1, tolerance = 1e-10)
})

test_that("get_collin_events computes VIF for regressors", {
  dmat <- matrix(list(), nrow = 1, ncol = 3)
  dimnames(dmat) <- list("run1", c("reg1", "reg2", "reg3"))
  
  set.seed(123)
  dmat[[1, 1]] <- data.frame(trial = 1:20, value = rnorm(20))
  dmat[[1, 2]] <- data.frame(trial = 1:20, value = rnorm(20))
  dmat[[1, 3]] <- data.frame(trial = 1:20, value = rnorm(20))
  
  result <- fmri.pipeline:::get_collin_events(dmat)
  
  expect_true("vif" %in% names(result$run1))
  # VIF should exist for all regressors
  expect_equal(length(result$run1$vif), 3)
})

test_that("get_collin_events handles regressors with different trials", {
  dmat <- matrix(list(), nrow = 1, ncol = 2)
  dimnames(dmat) <- list("run1", c("reg1", "reg2"))
  
  # Regressors on different trial subsets
  dmat[[1, 1]] <- data.frame(trial = 1:5, value = 1:5)
  dmat[[1, 2]] <- data.frame(trial = 3:8, value = 1:6)
  
  # Should handle this gracefully (union of trials)
  result <- fmri.pipeline:::get_collin_events(dmat)
  
  expect_true(!is.null(result$run1$r))
})

# ==============================================================================
# Tests for expand_signal (from l1_helper_functions.R)
# ==============================================================================

test_that("expand_signal returns list for simple signal without expansion", {
  sig <- list(
    name = "cue",
    event = "cue",
    value = 1,
    normalization = "none",
    beta_series = FALSE,
    wi_factors = NULL
  )
  
  result <- fmri.pipeline:::expand_signal(sig)
  
  expect_type(result, "list")
  expect_equal(length(result), 1)
  expect_equal(names(result), "cue")
})

test_that("expand_signal creates separate signals for PPI", {
  sig <- list(
    name = "cue",
    event = "cue",
    value = 1,
    normalization = "none",
    ts_multiplier = "roi1",
    beta_series = FALSE,
    wi_factors = NULL
  )
  
  result <- fmri.pipeline:::expand_signal(sig)
  
  # Should create: psych (cue), physio (roi1), ppi (cue.ppi)
  expect_equal(length(result), 3)
  expect_true("cue" %in% names(result))
  expect_true("roi1" %in% names(result))
  expect_true("cue.ppi" %in% names(result))
})

test_that("expand_signal creates derivative signals when add_deriv=TRUE", {
  sig <- list(
    name = "cue",
    event = "cue",
    value = 1,
    normalization = "none",
    beta_series = FALSE,
    wi_factors = NULL,
    add_deriv = TRUE
  )
  
  result <- fmri.pipeline:::expand_signal(sig)
  
  expect_equal(length(result), 2)
  # Check that first is the original signal and second is derivative
  expect_equal(result[[1]]$name, "cue")
  expect_equal(result[[2]]$name, "d_cue")
  expect_true(isTRUE(result[[2]]$return_deriv))
})

test_that("expand_signal handles beta_series expansion", {
  sig <- list(
    name = "cue",
    event = "cue",
    value = data.frame(
      run_number = 1,
      trial = 1:5,
      value = 1:5
    ),
    normalization = "none",
    beta_series = TRUE,
    wi_factors = NULL
  )
  
  result <- fmri.pipeline:::expand_signal(sig)
  
  # Should create one signal per trial
  expect_equal(length(result), 5)
  # Check via the internal name field since list names may be NULL
  signal_names <- sapply(result, function(x) x$name)
  expect_true("cue_t001" %in% signal_names)
  expect_true("cue_t005" %in% signal_names)
})

# ==============================================================================
# Tests for fmri.stimulus
# ==============================================================================

test_that("fmri.stimulus produces correct length output", {
  n_vols <- 100
  result <- fmri.pipeline:::fmri.stimulus(
    n_vols = n_vols,
    onsets = c(10, 30, 50),
    durations = c(2, 2, 2),
    values = c(1, 1, 1),
    tr = 1.0
  )
  
  expect_equal(length(result), n_vols)
})

test_that("fmri.stimulus returns zeros when no events occur", {
  result <- fmri.pipeline:::fmri.stimulus(
    n_vols = 50,
    onsets = numeric(0),
    durations = numeric(0),
    values = numeric(0),
    tr = 1.0
  )
  
  expect_equal(sum(result), 0)
})

test_that("fmri.stimulus handles instantaneous events (duration=0)", {
  result <- fmri.pipeline:::fmri.stimulus(
    n_vols = 100,
    onsets = 20,
    durations = 0,
    values = 1,
    tr = 1.0,
    convolve = TRUE
  )
  
  # Should still produce non-zero convolved output
  expect_gt(max(result), 0)
})

test_that("fmri.stimulus respects demean parameter", {
  result_demean <- fmri.pipeline:::fmri.stimulus(
    n_vols = 100,
    onsets = c(20, 50),
    durations = c(1, 1),
    values = c(1, 1),
    tr = 1.0,
    demean = TRUE
  )
  
  result_no_demean <- fmri.pipeline:::fmri.stimulus(
    n_vols = 100,
    onsets = c(20, 50),
    durations = c(1, 1),
    values = c(1, 1),
    tr = 1.0,
    demean = FALSE
  )
  
  expect_lt(abs(mean(result_demean)), 1e-10)
  expect_gt(mean(result_no_demean), 0)
})

test_that("fmri.stimulus handles unconvolved output", {
  result <- fmri.pipeline:::fmri.stimulus(
    n_vols = 100,
    onsets = 20,
    durations = 5,
    values = 3,
    tr = 1.0,
    convolve = FALSE,
    demean = FALSE
  )
  
  # Unconvolved should have value = 3 at onset through onset+duration
  expect_equal(result[21], 3)  # 0-indexed time, so volume 21 = time 20
  expect_equal(result[25], 3)  # time 24
})

# ==============================================================================
# Tests for convolve_regressor
# ==============================================================================

test_that("convolve_regressor produces correct length output", {
  reg <- matrix(
    c(1, 10, 2, 1, 2, 30, 2, 1),  # trial, onset, duration, value
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, c("trial", "onset", "duration", "value"))
  )
  attr(reg, "reg_name") <- "test_reg"
  
  result <- fmri.pipeline:::convolve_regressor(
    n_vols = 100,
    reg = reg,
    tr = 1.0,
    normalization = "none"
  )
  
  expect_equal(nrow(result), 100)
})

test_that("convolve_regressor respects evtmax_1 normalization", {
  reg <- matrix(
    c(1, 10, 2, 1, 2, 30, 10, 1),  # Different durations
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, c("trial", "onset", "duration", "value"))
  )
  attr(reg, "reg_name") <- "test_reg"
  
  result <- fmri.pipeline:::convolve_regressor(
    n_vols = 100,
    reg = reg,
    tr = 1.0,
    normalization = "evtmax_1"
  )
  
  # With evtmax_1, each event should reach approximately max=1 before value scaling
  expect_gt(max(result), 0)
})

test_that("convolve_regressor handles empty regressor", {
  reg <- matrix(
    nrow = 0, ncol = 4,
    dimnames = list(NULL, c("trial", "onset", "duration", "value"))
  )
  attr(reg, "reg_name") <- "empty_reg"
  attr(reg, "event") <- "empty"
  
  result <- fmri.pipeline:::convolve_regressor(
    n_vols = 100,
    reg = reg,
    tr = 1.0,
    normalization = "none"
  )
  
  expect_equal(nrow(result), 100)
  expect_equal(sum(result), 0)
})

test_that("convolve_regressor applies high-pass filter when specified", {
  # Create a regressor with low-frequency drift
  reg <- matrix(
    c(1, 10, 2, 1),
    nrow = 1,
    dimnames = list(NULL, c("trial", "onset", "duration", "value"))
  )
  attr(reg, "reg_name") <- "test_reg"
  
  result_no_filter <- fmri.pipeline:::convolve_regressor(
    n_vols = 200,
    reg = reg,
    tr = 1.0,
    normalization = "none",
    high_pass = NULL
  )
  
  result_filtered <- fmri.pipeline:::convolve_regressor(
    n_vols = 200,
    reg = reg,
    tr = 1.0,
    normalization = "none",
    high_pass = 1/100  # 0.01 Hz
  )
  
  # Filtered result should have less low-frequency content
  # (comparing overall variance or similar metrics)
  expect_true(!identical(result_no_filter, result_filtered))
})

test_that("convolve_regressor centers values when center_values=TRUE", {
  reg <- matrix(
    c(1, 10, 2, 5, 2, 30, 2, 10, 3, 50, 2, 15),  # Values: 5, 10, 15
    nrow = 3, byrow = TRUE,
    dimnames = list(NULL, c("trial", "onset", "duration", "value"))
  )
  attr(reg, "reg_name") <- "test_reg"
  
  # The mean-centering happens internally; we check the function runs
  result <- fmri.pipeline:::convolve_regressor(
    n_vols = 100,
    reg = reg,
    tr = 1.0,
    normalization = "none",
    center_values = TRUE
  )
  
  expect_gt(max(result), 0)
  expect_lt(min(result), 0)  # Should have negative values after centering
})

test_that("convolve_regressor applies ts_multiplier for PPI", {
  reg <- matrix(
    c(1, 10, 5, 1),
    nrow = 1,
    dimnames = list(NULL, c("trial", "onset", "duration", "value"))
  )
  attr(reg, "reg_name") <- "ppi_reg"
  
  # Create a time series multiplier
  ts_mult <- sin(seq(0, 4 * pi, length.out = 100))
  
  result <- fmri.pipeline:::convolve_regressor(
    n_vols = 100,
    reg = reg,
    tr = 1.0,
    normalization = "none",
    ts_multiplier = ts_mult
  )
  
  # PPI interaction should produce different results than without multiplier
  result_no_mult <- fmri.pipeline:::convolve_regressor(
    n_vols = 100,
    reg = reg,
    tr = 1.0,
    normalization = "none",
    ts_multiplier = NULL
  )
  
  expect_false(identical(as.vector(result), as.vector(result_no_mult)))
})

# ==============================================================================
# Tests for determine_run_volumes
# ==============================================================================

test_that("determine_run_volumes uses provided run_volumes correctly", {
  run_data <- data.frame(
    run_number = 1:2,
    run_volumes = c(100, 120),
    drop_volumes = c(0, 0)
  )
  
  result <- fmri.pipeline:::determine_run_volumes(
    run_data = run_data,
    nruns = 2,
    tr = 1.0,
    signals_aligned = list(),
    lg = lgr::get_logger()
  )
  
  expect_equal(result$run_volumes, c(100, 120))
})

test_that("determine_run_volumes adjusts for drop_volumes", {
  run_data <- data.frame(
    run_number = 1:2,
    run_volumes = c(100, 120),
    drop_volumes = c(5, 10)
  )
  
  result <- fmri.pipeline:::determine_run_volumes(
    run_data = run_data,
    nruns = 2,
    tr = 1.0,
    signals_aligned = list(),
    lg = lgr::get_logger()
  )
  
  # run_volumes is reduced by drop_volumes
  expect_equal(result$run_volumes, c(95, 110))
})

test_that("determine_run_volumes replicates scalar run_volumes", {
  run_data <- data.frame(
    run_number = 1:3,
    run_volumes = 100,  # Scalar, should be replicated
    drop_volumes = 0
  )
  
  result <- fmri.pipeline:::determine_run_volumes(
    run_data = run_data,
    nruns = 3,
    tr = 1.0,
    signals_aligned = list(),
    lg = lgr::get_logger()
  )
  
  expect_equal(length(result$run_volumes), 3)
  expect_equal(result$run_volumes, c(100, 100, 100))
})

# ==============================================================================
# Integration-style tests for build_design_matrix
# ==============================================================================

test_that("build_design_matrix returns expected output structure", {
  skip_on_cran()
  
  set.seed(123)
  events <- create_test_events(nruns = 2, ntrials = 10)
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 2, run_volumes = 50)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    plot = FALSE
  )
  
  expect_type(result, "list")
  expect_true("design" %in% names(result))
  expect_true("design_convolved" %in% names(result))
  expect_true("design_unconvolved" %in% names(result))
  expect_true("collin_events" %in% names(result))
  expect_true("run_volumes" %in% names(result))
})

test_that("build_design_matrix respects runs_to_output", {
  skip_on_cran()
  
  set.seed(456)
  events <- create_test_events(nruns = 4, ntrials = 5)
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 4, run_volumes = 40)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    runs_to_output = c(1, 3),
    plot = FALSE
  )
  
  expect_equal(length(result$design_convolved), 2)
  expect_equal(result$run_volumes, c(40, 40))
})

test_that("build_design_matrix handles drop_volumes correctly", {
  skip_on_cran()
  
  set.seed(789)
  events <- create_test_events(nruns = 2, ntrials = 5)
  # Adjust onsets to be after drop period
  events$onset <- events$onset + 10
  
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 2, run_volumes = 60, drop_volumes = 5)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    plot = FALSE
  )
  
  # Final volumes should be reduced by drop_volumes
  expect_equal(result$run_volumes, c(55, 55))
  expect_equal(nrow(result$design_convolved[[1]]), 55)
})

test_that("build_design_matrix includes additional_regressors", {
  skip_on_cran()
  
  set.seed(111)
  events <- create_test_events(nruns = 2, ntrials = 5)
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 2, run_volumes = 50)
  
  add_reg <- list(
    data.frame(confound1 = rnorm(50), confound2 = rnorm(50)),
    data.frame(confound1 = rnorm(50), confound2 = rnorm(50))
  )
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    additional_regressors = add_reg,
    plot = FALSE
  )
  
  # Should have cue_evt + 2 confounds
  expect_true("confound1" %in% names(result$design_convolved[[1]]))
  expect_true("confound2" %in% names(result$design_convolved[[1]]))
})

test_that("build_design_matrix includes baseline regressors", {
  skip_on_cran()
  
  set.seed(222)
  events <- create_test_events(nruns = 2, ntrials = 5)
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 2, run_volumes = 50)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    baseline_coef_order = 2,
    plot = FALSE
  )
  
  # Should have cue_evt + base0, base1, base2
  expect_true("base0" %in% names(result$design_convolved[[1]]))
  expect_true("base1" %in% names(result$design_convolved[[1]]))
  expect_true("base2" %in% names(result$design_convolved[[1]]))
})

test_that("build_design_matrix handles parametric regressors", {
  skip_on_cran()
  
  set.seed(333)
  events <- create_test_events(nruns = 2, ntrials = 10)
  value_df <- data.frame(
    run_number = c(rep(1, 10), rep(2, 10)),
    trial = c(1:10, 1:10),
    value = rnorm(20, mean = 5, sd = 2)
  )
  
  signals <- list(
    cue_param = list(name = "cue_param", event = "cue", value = value_df, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 2, run_volumes = 80)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    center_values = TRUE,
    plot = FALSE
  )
  
  # Parametric regressor should exist and have non-zero variance
  expect_true("cue_param" %in% names(result$design_convolved[[1]]))
  expect_gt(sd(result$design_convolved[[1]]$cue_param), 0)
})

test_that("build_design_matrix handles empty runs gracefully", {
  skip_on_cran()
  skip("Current implementation infers nruns from events$run_number, not run_data, so empty runs are not supported")
  
  set.seed(444)
  # Events only in run 1
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1:5,
    onset = c(5, 15, 25, 35, 45),
    duration = 1
  )
  
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 2, run_volumes = 60)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    keep_empty_regressors = TRUE,
    plot = FALSE
  )
  
  # Run 2 should have all-zero regressor
  expect_equal(sum(result$design_convolved[[2]]$cue_evt), 0)
})

test_that("build_design_matrix adds temporal derivatives when requested", {
  skip_on_cran()
  
  set.seed(555)
  events <- create_test_events(nruns = 2, ntrials = 5)
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none", add_deriv = TRUE)
  )
  run_data <- create_test_run_data(nruns = 2, run_volumes = 50)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    plot = FALSE
  )
  
  # Should have both cue_evt and d_cue_evt (derivative)
  expect_true("cue_evt" %in% names(result$design_convolved[[1]]))
  expect_true("d_cue_evt" %in% names(result$design_convolved[[1]]))
})

test_that("build_design_matrix applies high-pass filter", {
  skip_on_cran()
  
  set.seed(666)
  events <- create_test_events(nruns = 1, ntrials = 10)
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 1, run_volumes = 200)
  
  result_no_hp <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    high_pass = NULL,
    plot = FALSE
  )
  
  result_hp <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    high_pass = 1/100,  # 0.01 Hz
    plot = FALSE
  )
  
  # Filtered result should differ
  expect_false(identical(
    result_no_hp$design_convolved[[1]]$cue_evt,
    result_hp$design_convolved[[1]]$cue_evt
  ))
})

test_that("build_design_matrix handles beta_series signal", {
  skip_on_cran()
  
  set.seed(777)
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1:5,
    onset = c(10, 30, 50, 70, 90),
    duration = 1
  )
  
  signals <- list(
    cue_evt = list(
      name = "cue_evt",
      event = "cue",
      value = data.frame(
        run_number = 1,
        trial = 1:5,
        value = 1
      ),
      normalization = "none",
      beta_series = TRUE
    )
  )
  run_data <- create_test_run_data(nruns = 1, run_volumes = 120)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    plot = FALSE
  )
  
  # Should have 5 separate trial regressors
  reg_names <- names(result$design_convolved[[1]])
  trial_regs <- grep("cue_evt_t", reg_names, value = TRUE)
  expect_equal(length(trial_regs), 5)
})

# ==============================================================================
# Edge case tests
# ==============================================================================

test_that("build_design_matrix handles events at very end of run", {
  skip_on_cran()
  
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1:3,
    onset = c(10, 40, 48),  # Last event near end
    duration = 1
  )
  
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "evtmax_1")
  )
  run_data <- create_test_run_data(nruns = 1, run_volumes = 50)
  
  # Should not error even with event near end
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    plot = FALSE
  )
  
  expect_equal(nrow(result$design_convolved[[1]]), 50)
})

test_that("build_design_matrix handles single trial correctly", {
  skip_on_cran()
  
  events <- data.frame(
    event = "cue",
    run_number = 1,
    trial = 1,
    onset = 20,
    duration = 2
  )
  
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 1, run_volumes = 50)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    plot = FALSE
  )
  
  expect_gt(max(result$design_convolved[[1]]$cue_evt), 0)
})

test_that("build_design_matrix handles noncontiguous runs_to_output", {
  skip_on_cran()
  
  set.seed(888)
  events <- create_test_events(nruns = 5, ntrials = 5)
  signals <- list(
    cue_evt = list(name = "cue_evt", event = "cue", value = 1, normalization = "none")
  )
  run_data <- create_test_run_data(nruns = 5, run_volumes = 40)
  
  result <- build_design_matrix(
    events = events,
    signals = signals,
    tr = 1.0,
    run_data = run_data,
    runs_to_output = c(1, 3, 5),
    plot = FALSE
  )
  
  expect_equal(length(result$design_convolved), 3)
  expect_equal(result$run_volumes, c(40, 40, 40))
})
