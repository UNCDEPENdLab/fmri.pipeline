# Unit tests for collinearity_diagnostics.R functions
# Uses mock factories from helpers-mock_factories.R

# Note: create_mock_collinearity_summary() and create_empty_collinearity_summary()
# are defined in helpers-mock_factories.R and loaded automatically by testthat

# ==============================================================================
# Tests for l1_collinearity_summary object structure
# ==============================================================================

test_that("mock collinearity summary has correct structure", {
  mock <- create_mock_collinearity_summary()
  
  expect_s3_class(mock, "l1_collinearity_summary")
  expect_s3_class(mock, "list")
  
  # Check all expected components exist
  expect_true("vif_summary" %in% names(mock))
  expect_true("correlation_summary" %in% names(mock))
  expect_true("high_vif" %in% names(mock))
  expect_true("high_correlation" %in% names(mock))
  expect_true("summary_stats" %in% names(mock))
  expect_true("thresholds" %in% names(mock))
})

test_that("mock vif_summary has correct columns", {
  mock <- create_mock_collinearity_summary()
  
  expected_cols <- c("id", "session", "run_number", "l1_model", "regressor", "vif", "high_vif")
  expect_true(all(expected_cols %in% names(mock$vif_summary)))
})

test_that("mock correlation_summary has correct columns", {
  mock <- create_mock_collinearity_summary()
  
  expected_cols <- c("id", "session", "run_number", "l1_model", 
                     "regressor1", "regressor2", "correlation", 
                     "abs_correlation", "high_correlation")
  expect_true(all(expected_cols %in% names(mock$correlation_summary)))
})

# ==============================================================================
# Tests for print.l1_collinearity_summary
# ==============================================================================

test_that("print.l1_collinearity_summary prints without error", {
  mock <- create_mock_collinearity_summary()
  
  expect_output(print(mock), "L1 Collinearity Diagnostics Summary")
  expect_output(print(mock), "Variance Inflation Factors")
  expect_output(print(mock), "Regressor Correlations")
})

test_that("print.l1_collinearity_summary shows correct counts", {
  mock <- create_mock_collinearity_summary(n_subjects = 5, n_runs = 3)
  
  expect_output(print(mock), "Subjects analyzed: 5")
  # Total runs = n_subjects * n_runs = 15
  expect_output(print(mock), "Total runs: 15")
})

test_that("print.l1_collinearity_summary handles empty data", {
  mock <- create_empty_collinearity_summary()
  
  expect_output(print(mock), "L1 Collinearity Diagnostics Summary")
  expect_output(print(mock), "Subjects analyzed: 0")
})

test_that("print.l1_collinearity_summary returns invisible self", {
  mock <- create_mock_collinearity_summary()
  
  result <- capture.output(ret_val <- print(mock))
  expect_identical(ret_val, mock)
})

# ==============================================================================
# Tests for summarize_vif_by_regressor
# ==============================================================================

test_that("summarize_vif_by_regressor returns data.frame", {
  mock <- create_mock_collinearity_summary()
  
  result <- summarize_vif_by_regressor(mock)
  
  expect_s3_class(result, "data.frame")
})

test_that("summarize_vif_by_regressor has correct columns", {
  mock <- create_mock_collinearity_summary()
  
  result <- summarize_vif_by_regressor(mock)
  
  expected_cols <- c("regressor", "mean_vif", "sd_vif", "min_vif", "max_vif",
                     "n_runs", "n_flagged", "pct_flagged")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("summarize_vif_by_regressor counts runs correctly", {
  mock <- create_mock_collinearity_summary(n_subjects = 2, n_runs = 3)
  
  result <- summarize_vif_by_regressor(mock)
  
  # Each regressor should appear in all subject-run combos
  # n_runs per regressor = n_subjects * n_runs_per_subject = 2 * 3 = 6
  expect_true(all(result$n_runs == 6))
})

test_that("summarize_vif_by_regressor with by_model=TRUE includes model column", {
  mock <- create_mock_collinearity_summary()
  
  result <- summarize_vif_by_regressor(mock, by_model = TRUE)
  
  expect_true("l1_model" %in% names(result))
})

test_that("summarize_vif_by_regressor handles empty data", {
  mock <- create_empty_collinearity_summary()
  
  result <- summarize_vif_by_regressor(mock)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("summarize_vif_by_regressor rejects invalid input", {
  expect_error(summarize_vif_by_regressor(list()), "l1_collinearity_summary")
  expect_error(summarize_vif_by_regressor(NULL), "l1_collinearity_summary")
})

test_that("summarize_vif_by_regressor orders by mean_vif descending", {
  mock <- create_mock_collinearity_summary(include_high_vif = TRUE)
  
  result <- summarize_vif_by_regressor(mock)
  
  # First row should have highest mean_vif
  if (nrow(result) > 1) {
    expect_gte(result$mean_vif[1], result$mean_vif[2])
  }
})

# ==============================================================================
# Tests for summarize_correlations_by_pair
# ==============================================================================

test_that("summarize_correlations_by_pair returns data.frame", {
  mock <- create_mock_collinearity_summary()
  
  result <- summarize_correlations_by_pair(mock)
  
  expect_s3_class(result, "data.frame")
})

test_that("summarize_correlations_by_pair has correct columns", {
  mock <- create_mock_collinearity_summary()
  
  result <- summarize_correlations_by_pair(mock)
  
  expected_cols <- c("regressor1", "regressor2", "mean_cor", "sd_cor",
                     "min_cor", "max_cor", "mean_abs_cor", 
                     "n_runs", "n_flagged", "pct_flagged")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("summarize_correlations_by_pair with by_model=TRUE includes model column", {
  mock <- create_mock_collinearity_summary()
  
  result <- summarize_correlations_by_pair(mock, by_model = TRUE)
  
  expect_true("l1_model" %in% names(result))
})

test_that("summarize_correlations_by_pair handles empty data", {
  mock <- create_empty_collinearity_summary()
  
  result <- summarize_correlations_by_pair(mock)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("summarize_correlations_by_pair rejects invalid input", {
  expect_error(summarize_correlations_by_pair(list()), "l1_collinearity_summary")
  expect_error(summarize_correlations_by_pair(NULL), "l1_collinearity_summary")
})

test_that("summarize_correlations_by_pair orders by mean_abs_cor descending", {
  mock <- create_mock_collinearity_summary(include_high_cor = TRUE)
  
  result <- summarize_correlations_by_pair(mock)
  
  # First row should have highest mean_abs_cor
  if (nrow(result) > 1) {
    expect_gte(result$mean_abs_cor[1], result$mean_abs_cor[2])
  }
})

test_that("summarize_correlations_by_pair has correct number of pairs", {
  n_reg <- 4
  mock <- create_mock_collinearity_summary(n_regressors = n_reg)
  
  result <- summarize_correlations_by_pair(mock)
  
  # Number of unique pairs = n_reg choose 2 = n_reg * (n_reg - 1) / 2
  expected_pairs <- n_reg * (n_reg - 1) / 2
  expect_equal(nrow(result), expected_pairs)
})

# ==============================================================================
# Tests for export_collinearity_to_csv
# ==============================================================================

test_that("export_collinearity_to_csv writes files", {
  mock <- create_mock_collinearity_summary()
  temp_dir <- tempdir()
  
  files_written <- export_collinearity_to_csv(mock, output_dir = temp_dir, prefix = "test_")
  
  expect_gt(length(files_written), 0)
  expect_true(all(file.exists(files_written)))
  
  # Clean up
  file.remove(files_written)
})

test_that("export_collinearity_to_csv uses prefix correctly", {
  mock <- create_mock_collinearity_summary()
  temp_dir <- tempdir()
  
  files_written <- export_collinearity_to_csv(mock, output_dir = temp_dir, prefix = "myprefix_")
  
  expect_true(all(grepl("^myprefix_", basename(files_written))))
  
  # Clean up
  file.remove(files_written)
})

test_that("export_collinearity_to_csv writes expected files", {
  mock <- create_mock_collinearity_summary(include_high_vif = TRUE, include_high_cor = TRUE)
  temp_dir <- tempdir()
  
  files_written <- export_collinearity_to_csv(mock, output_dir = temp_dir, prefix = "export_test_")
  
  basenames <- basename(files_written)
  
  # Should include main summary files
  expect_true("export_test_vif_all.csv" %in% basenames)
  expect_true("export_test_correlations_all.csv" %in% basenames)
  
  # Should include flagged files (since we have high values)
  expect_true("export_test_high_vif.csv" %in% basenames)
  expect_true("export_test_high_correlations.csv" %in% basenames)
  
  # Should include regressor summaries
  expect_true("export_test_vif_by_regressor.csv" %in% basenames)
  expect_true("export_test_correlations_by_pair.csv" %in% basenames)
  
  # Clean up
  file.remove(files_written)
})

test_that("export_collinearity_to_csv handles empty data", {
  mock <- create_empty_collinearity_summary()
  temp_dir <- tempdir()
  
  files_written <- export_collinearity_to_csv(mock, output_dir = temp_dir, prefix = "empty_test_")
  
  # Should return empty vector for empty data

  expect_equal(length(files_written), 0)
})

test_that("export_collinearity_to_csv rejects invalid directory", {
  mock <- create_mock_collinearity_summary()
  
  expect_error(export_collinearity_to_csv(mock, output_dir = "/nonexistent/path"))
})

test_that("export_collinearity_to_csv rejects invalid input", {
  temp_dir <- tempdir()
  
  expect_error(export_collinearity_to_csv(list(), output_dir = temp_dir))
  expect_error(export_collinearity_to_csv(NULL, output_dir = temp_dir))
})

test_that("exported CSV files are readable and have correct structure", {
  mock <- create_mock_collinearity_summary(n_subjects = 2, n_runs = 2, n_regressors = 3)
  temp_dir <- tempdir()
  
  files_written <- export_collinearity_to_csv(mock, output_dir = temp_dir, prefix = "readable_test_")
  
  # Find and read the VIF all file
  vif_file <- files_written[grepl("vif_all.csv", files_written)]
  expect_length(vif_file, 1)
  
  vif_data <- read.csv(vif_file)
  expect_equal(nrow(vif_data), nrow(mock$vif_summary))
  expect_true("regressor" %in% names(vif_data))
  expect_true("vif" %in% names(vif_data))
  
  # Clean up
  file.remove(files_written)
})

# ==============================================================================
# Tests for plot_collinearity
# ==============================================================================

test_that("plot_collinearity returns ggplot for vif type", {
  skip_if_not_installed("ggplot2")
  
  mock <- create_mock_collinearity_summary()
  
  result <- plot_collinearity(mock, plot_type = "vif")
  
  expect_s3_class(result, "ggplot")
})

test_that("plot_collinearity returns ggplot for correlation type", {
  skip_if_not_installed("ggplot2")
  
  mock <- create_mock_collinearity_summary()
  
  result <- plot_collinearity(mock, plot_type = "correlation")
  
  expect_s3_class(result, "ggplot")
})

test_that("plot_collinearity returns list for both type", {
  skip_if_not_installed("ggplot2")
  
  mock <- create_mock_collinearity_summary()
  
  result <- plot_collinearity(mock, plot_type = "both")
  
  expect_type(result, "list")
  expect_true("vif" %in% names(result))
  expect_true("correlation" %in% names(result))
  expect_s3_class(result$vif, "ggplot")
  expect_s3_class(result$correlation, "ggplot")
})

test_that("plot_collinearity rejects invalid plot_type", {
  mock <- create_mock_collinearity_summary()
  
  expect_error(plot_collinearity(mock, plot_type = "invalid"))
  expect_error(plot_collinearity(mock, plot_type = "vifs"))
})

test_that("plot_collinearity rejects invalid input", {
  expect_error(plot_collinearity(list()))
  expect_error(plot_collinearity(NULL))
})

test_that("plot_collinearity handles empty vif data gracefully", {
  skip_if_not_installed("ggplot2")
  
  mock <- create_empty_collinearity_summary()
  
  # Should return empty list when no data
  result <- plot_collinearity(mock, plot_type = "both")
  
  expect_type(result, "list")
  expect_equal(length(result), 0)
})

# ==============================================================================
# Tests for flagging logic
# ==============================================================================

test_that("high VIF values are correctly flagged", {
  # Custom thresholds
  vif_thresh <- 3
  mock <- create_mock_collinearity_summary(
    n_subjects = 1, 
    n_runs = 1, 
    vif_threshold = vif_thresh,
    include_high_vif = TRUE
  )
  
  # We set first regressor to vif_threshold + 2 = 5 when include_high_vif = TRUE
  expect_true(any(mock$vif_summary$vif > vif_thresh))
  expect_true(any(mock$vif_summary$high_vif))
  expect_gt(nrow(mock$high_vif), 0)
})

test_that("high correlation values are correctly flagged", {
  cor_thresh <- 0.7
  mock <- create_mock_collinearity_summary(
    n_subjects = 1,
    n_runs = 1,
    cor_threshold = cor_thresh,
    include_high_cor = TRUE
  )
  
  # We set first pair to cor_threshold + 0.1 when include_high_cor = TRUE
  expect_true(any(mock$correlation_summary$abs_correlation > cor_thresh))
  expect_true(any(mock$correlation_summary$high_correlation))
  expect_gt(nrow(mock$high_correlation), 0)
})

test_that("no false positives when values are below threshold", {
  mock <- create_mock_collinearity_summary(
    n_subjects = 2,
    n_runs = 2,
    vif_threshold = 10,  # Very high threshold
    cor_threshold = 0.99,  # Very high threshold
    include_high_vif = FALSE,
    include_high_cor = FALSE
  )
  
  expect_equal(nrow(mock$high_vif), 0)
  expect_equal(nrow(mock$high_correlation), 0)
})

# ==============================================================================
# Tests for summary statistics
# ==============================================================================

test_that("summary_stats has correct counts", {
  n_subj <- 4
  n_run <- 3
  n_reg <- 5
  
  mock <- create_mock_collinearity_summary(
    n_subjects = n_subj,
    n_runs = n_run,
    n_regressors = n_reg
  )
  
  expect_equal(mock$summary_stats$n_subjects, n_subj)
  expect_equal(mock$summary_stats$n_runs, n_subj * n_run)
  expect_equal(mock$summary_stats$n_regressors, n_reg)
})

test_that("summary_stats VIF values are correct", {
  mock <- create_mock_collinearity_summary(n_subjects = 2, n_runs = 2)
  
  expected_mean <- mean(mock$vif_summary$vif, na.rm = TRUE)
  expected_max <- max(mock$vif_summary$vif, na.rm = TRUE)
  
  expect_equal(mock$summary_stats$vif_mean, expected_mean)
  expect_equal(mock$summary_stats$vif_max, expected_max)
})

test_that("summary_stats correlation values are correct", {
  mock <- create_mock_collinearity_summary(n_subjects = 2, n_runs = 2)
  
  expected_mean_abs <- mean(mock$correlation_summary$abs_correlation, na.rm = TRUE)
  expected_max_abs <- max(mock$correlation_summary$abs_correlation, na.rm = TRUE)
  
  expect_equal(mock$summary_stats$cor_mean_abs, expected_mean_abs)
  expect_equal(mock$summary_stats$cor_max_abs, expected_max_abs)
})

# ==============================================================================
# Tests for thresholds
# ==============================================================================

test_that("thresholds are stored correctly", {
  vif_t <- 7
  cor_t <- 0.75
  
  mock <- create_mock_collinearity_summary(
    vif_threshold = vif_t,
    cor_threshold = cor_t
  )
  
  expect_equal(mock$thresholds$vif, vif_t)
  expect_equal(mock$thresholds$correlation, cor_t)
})

# ==============================================================================
# Edge case tests
# ==============================================================================

test_that("single subject single run works correctly", {
  mock <- create_mock_collinearity_summary(n_subjects = 1, n_runs = 1)
  
  expect_s3_class(mock, "l1_collinearity_summary")
  expect_equal(mock$summary_stats$n_subjects, 1)
  expect_equal(mock$summary_stats$n_runs, 1)
  
  # Should be able to summarize
  vif_summ <- summarize_vif_by_regressor(mock)
  expect_s3_class(vif_summ, "data.frame")
  
  cor_summ <- summarize_correlations_by_pair(mock)
  expect_s3_class(cor_summ, "data.frame")
})

test_that("single regressor works correctly (no pairs)", {
  # With only 1 regressor, there are no pairs for correlation
  mock <- create_mock_collinearity_summary(
    n_subjects = 2,
    n_runs = 2,
    n_regressors = 1,
    include_high_vif = FALSE,
    include_high_cor = FALSE
  )
  
  # VIF should work fine
  expect_gt(nrow(mock$vif_summary), 0)
  
  # Correlation summary should be empty (no pairs with 1 regressor)
  expect_equal(nrow(mock$correlation_summary), 0)
})

test_that("two regressors produces exactly one pair", {
  mock <- create_mock_collinearity_summary(
    n_subjects = 2,
    n_runs = 2,
    n_regressors = 2,
    include_high_vif = FALSE,
    include_high_cor = FALSE
  )
  
  cor_summ <- summarize_correlations_by_pair(mock)
  
  # With 2 regressors, there should be exactly 1 unique pair
  expect_equal(nrow(cor_summ), 1)
})

# ==============================================================================
# Integration-style tests (using mock objects)
# ==============================================================================

test_that("full workflow executes without error", {
  skip_if_not_installed("ggplot2")
  
  # Create mock
  mock <- create_mock_collinearity_summary(
    n_subjects = 3,
    n_runs = 4,
    n_regressors = 5,
    include_high_vif = TRUE,
    include_high_cor = TRUE
  )
  
  # Print
  expect_output(print(mock), "L1 Collinearity Diagnostics Summary")
  
  # Summarize VIF
  vif_summ <- summarize_vif_by_regressor(mock)
  expect_s3_class(vif_summ, "data.frame")
  expect_gt(nrow(vif_summ), 0)
  
  vif_summ_model <- summarize_vif_by_regressor(mock, by_model = TRUE)
  expect_s3_class(vif_summ_model, "data.frame")
  expect_true("l1_model" %in% names(vif_summ_model))
  
  # Summarize correlations
  cor_summ <- summarize_correlations_by_pair(mock)
  expect_s3_class(cor_summ, "data.frame")
  expect_gt(nrow(cor_summ), 0)
  
  cor_summ_model <- summarize_correlations_by_pair(mock, by_model = TRUE)
  expect_s3_class(cor_summ_model, "data.frame")
  expect_true("l1_model" %in% names(cor_summ_model))
  
  # Export
  temp_dir <- tempdir()
  files <- export_collinearity_to_csv(mock, output_dir = temp_dir, prefix = "workflow_")
  expect_gt(length(files), 0)
  file.remove(files)
  
  # Plot
  plots <- plot_collinearity(mock, plot_type = "both")
  expect_type(plots, "list")
  expect_length(plots, 2)
})
