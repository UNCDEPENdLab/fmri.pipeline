test_that("jlmer validates call-site inputs before starting Julia", {
  testthat::skip_if_not_installed("JuliaCall")

  df <- data.frame(y = 1:4, x = 1:4, id = factor(c(1, 1, 2, 2)))

  expect_error(
    jlmer(~ x + (1 | id), df),
    "two-sided"
  )

  expect_error(
    jlmer(y ~ missing_x + (1 | id), df),
    "Variables not found in data"
  )

  expect_error(
    jlmer(y ~ x + (1 | id), df, bad.arg = 1),
    "valid Julia identifiers"
  )

  expect_error(
    jlmer(y ~ x + (1 | id), df, lmer_test = "yes"),
    "lmer_test"
  )

  df$y[1] <- NA
  expect_error(
    jlmer(y ~ x + (1 | id), df, na.action = stats::na.fail),
    "Missing values"
  )
})

test_that("mixed_by jlmer engine rejects custom jlmer functions", {
  df <- data.frame(
    id = factor(rep(1:2, each = 2)),
    split = rep(c("a", "b"), each = 2),
    x = rnorm(4),
    y = rnorm(4)
  )

  expect_error(
    mixed_by(
      df,
      outcomes = "y",
      rhs_model_formulae = list(m1 = ~ x + (1 | id)),
      split_on = "split",
      engine = "jlmer",
      engine_args = list(jlmer_fun = function(...) NULL),
      calculate = "parameter_estimates_reml",
      padjust_by = NULL
    ),
    "jlmer_fun is no longer supported"
  )
})
