test_that("deconvolve_reglin recovers a continuous latent signal", {
  set.seed(1001)
  kernel <- c(0, 0.25, 0.5, 0.2, 0.05)
  n_time <- 80L
  n_latent <- n_time + length(kernel) - 1L
  latent <- as.numeric(stats::arima.sim(model = list(ar = 0.6), n = n_latent))
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, kernel) %*% latent

  estimated <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 1e-4,
    penalty = "diff1",
    normalize = FALSE,
    demean = FALSE
  )

  expect_equal(dim(estimated), c(n_time, 1L))
  expect_gt(stats::cor(estimated[, 1L], latent[length(kernel):n_latent]), 0.95)
})

test_that("deconvolve_reglin tunes shifted HRF candidates", {
  set.seed(1002)
  kernel <- c(0, 0.25, 0.5, 0.2, 0.05)
  true_lag <- 1L
  true_kernel <- fmri.pipeline:::shift_decon_kernel(true_lag, kernel)
  n_time <- 70L
  n_latent <- n_time + length(kernel) - 1L
  latent <- sin(seq(0, 6 * pi, length.out = n_latent)) + stats::rnorm(n_latent, sd = 0.05)
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, true_kernel) %*% latent

  result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 1e-4,
    penalty = "diff2",
    hrf_lags = c(-1L, 0L, 1L),
    tune_by = "global",
    normalize = FALSE,
    demean = FALSE,
    return_diagnostics = TRUE
  )

  expect_equal(result$kernel_index, 3L)
  expect_gt(stats::cor(result$activity[, 1L], latent[length(kernel):n_latent]), 0.9)
})

test_that("deconvolve_reglin C++ backend matches R backend", {
  set.seed(1003)
  kernel <- c(0, 0.25, 0.5, 0.2, 0.05)
  bold <- matrix(stats::rnorm(90), ncol = 3L)

  cpp_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 0.01,
    penalty = "diff2",
    hrf_lags = c(0L, 1L),
    normalize = TRUE,
    return_diagnostics = TRUE,
    backend = "cpp"
  )

  r_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 0.01,
    penalty = "diff2",
    hrf_lags = c(0L, 1L),
    normalize = TRUE,
    return_diagnostics = TRUE,
    backend = "r"
  )

  expect_equal(cpp_result$activity, r_result$activity, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(cpp_result$lambda, r_result$lambda, tolerance = 1e-12)
  expect_equal(cpp_result$kernel_index, r_result$kernel_index)
})

test_that("deconvolve_reglin prewhitened C++ backend matches R backend", {
  set.seed(1004)
  kernel <- c(0, 0.25, 0.5, 0.2, 0.05)
  bold <- matrix(stats::filter(stats::rnorm(120), filter = 0.6, method = "recursive"), ncol = 3L)

  cpp_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = c(0.001, 0.01, 0.1),
    penalty = "diff2",
    hrf_lags = c(0L, 1L),
    normalize = TRUE,
    prewhiten_gcv = TRUE,
    return_diagnostics = TRUE,
    backend = "cpp"
  )

  r_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = c(0.001, 0.01, 0.1),
    penalty = "diff2",
    hrf_lags = c(0L, 1L),
    normalize = TRUE,
    prewhiten_gcv = TRUE,
    return_diagnostics = TRUE,
    backend = "r"
  )

  expect_equal(cpp_result$activity, r_result$activity, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(cpp_result$ar1_rho, r_result$ar1_rho, tolerance = 1e-8)
  expect_true(cpp_result$prewhiten_gcv)
})

test_that("deconvolve_reglin Wiener backend recovers frequency-domain latent signal", {
  set.seed(1005)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 90L
  n_latent <- n_time + length(kernel) - 1L
  latent <- as.numeric(stats::arima.sim(model = list(ar = 0.5), n = n_latent))
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, kernel) %*% latent
  bold <- bold + stats::rnorm(n_time, sd = stats::sd(bold) / 20)

  result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda = 1e-4,
    normalize = FALSE,
    demean = FALSE,
    backend = "wiener",
    return_diagnostics = TRUE
  )

  expect_equal(dim(result$activity), c(n_time, 1L))
  expect_equal(result$penalty, "wiener")
  expect_gt(stats::cor(result$activity[, 1L], latent[length(kernel):n_latent]), 0.85)
})

test_that("deconvolve_reglin Wiener backend tunes lambda", {
  set.seed(1006)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 80L
  n_latent <- n_time + length(kernel) - 1L
  latent <- sin(seq(0, 4 * pi, length.out = n_latent))
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, kernel) %*% latent
  bold <- bold + stats::rnorm(n_time, sd = stats::sd(bold) / 4)

  result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = c(1e-5, 1e-3, 1e-1),
    normalize = FALSE,
    demean = FALSE,
    backend = "wiener",
    return_diagnostics = TRUE
  )

  expect_true(result$lambda %in% c(1e-5, 1e-3, 1e-1))
  expect_equal(nrow(result$tuning), 3L)
})

test_that("deconvolve_reglin PSD Wiener backend recovers latent signal", {
  set.seed(1007)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 90L
  n_latent <- n_time + length(kernel) - 1L
  latent <- as.numeric(stats::arima.sim(model = list(ar = 0.6), n = n_latent))
  bold <- fmri.pipeline:::decon_convolution_matrix(n_time, kernel) %*% latent
  bold <- bold + as.numeric(stats::arima.sim(model = list(ar = 0.5), n = n_time)) * stats::sd(bold) / 8

  result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = c(1e-4, 1e-2, 1),
    normalize = FALSE,
    demean = FALSE,
    backend = "wiener",
    wiener_shrinkage = "psd",
    return_diagnostics = TRUE
  )

  expect_equal(dim(result$activity), c(n_time, 1L))
  expect_equal(result$wiener_shrinkage, "psd")
  expect_equal(nrow(result$tuning), 3L)
  expect_gt(stats::cor(result$activity[, 1L], latent[length(kernel):n_latent]), 0.75)
})

test_that("deconvolve_reglin one-standard-error GCV rule selects smoother C++ candidate", {
  set.seed(1008)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 70L
  bold <- matrix(stats::rnorm(n_time * 12L), nrow = n_time)
  lambda_grid <- c(1e-4, 1e-2, 1)

  min_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = lambda_grid,
    normalize = TRUE,
    backend = "cpp",
    return_diagnostics = TRUE,
    gcv_rule = "min"
  )

  one_se_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = lambda_grid,
    normalize = TRUE,
    backend = "cpp",
    return_diagnostics = TRUE,
    gcv_rule = "1se"
  )

  expect_equal(dim(one_se_result$activity), c(n_time, ncol(bold)))
  expect_true("gcv_se" %in% names(one_se_result$tuning))
  expect_gte(one_se_result$lambda, min_result$lambda)
  expect_equal(one_se_result$gcv_rule, "1se")
})

test_that("deconvolve_reglin one-standard-error GCV rule works for Wiener backend", {
  set.seed(1009)
  kernel <- c(0.05, 0.25, 0.45, 0.2, 0.05)
  n_time <- 70L
  bold <- matrix(stats::rnorm(n_time * 10L), nrow = n_time)
  lambda_grid <- c(1e-4, 1e-2, 1)

  min_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = lambda_grid,
    normalize = TRUE,
    backend = "wiener",
    return_diagnostics = TRUE,
    gcv_rule = "min"
  )

  one_se_result <- deconvolve_reglin(
    BOLDobs = bold,
    kernel = kernel,
    lambda_grid = lambda_grid,
    normalize = TRUE,
    backend = "wiener",
    return_diagnostics = TRUE,
    gcv_rule = "1se"
  )

  expect_equal(dim(one_se_result$activity), c(n_time, ncol(bold)))
  expect_true("gcv_se" %in% names(one_se_result$tuning))
  expect_gte(one_se_result$lambda, min_result$lambda)
  expect_equal(one_se_result$gcv_rule, "1se")
})
