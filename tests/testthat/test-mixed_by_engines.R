test_that("mixed_by rstanarm engine adds pd and posterior p-values", {
  testthat::skip_if_not_installed("rstanarm")
  testthat::skip_if_not_installed("bayestestR")

  set.seed(202)
  df <- data.frame(
    id = factor(rep(rep(1:4, each = 3), times = 2)),
    split = rep(c("a", "b"), each = 12),
    x = rnorm(24)
  )
  df$y <- 0.5 + 0.35 * df$x + rep(rnorm(8, sd = 0.2), each = 3) + rnorm(24, sd = 0.4)

  res <- suppressMessages(suppressWarnings(mixed_by(
    df,
    outcomes = "y",
    rhs_model_formulae = list(m1 = ~ x + (1 | id)),
    split_on = "split",
    engine = "rstanarm",
    engine_args = list(chains = 1, iter = 100, warmup = 50, seed = 202, refresh = 0),
    calculate = "parameter_estimates_reml"
  )))

  expect_true(all(c("pd", "p.value") %in% names(res$coef_df_reml)))
  expect_true("padj_BY_term" %in% names(res$coef_df_reml))
  expect_false(anyNA(res$coef_df_reml$pd))
  expect_false(anyNA(res$coef_df_reml$p.value))
  expect_true(all(res$coef_df_reml$pd >= 0.5 & res$coef_df_reml$pd <= 1))
  expect_true(all(res$coef_df_reml$p.value >= 0 & res$coef_df_reml$p.value <= 1))
  expect_equal(res$coef_df_reml$p.value, 2 * (1 - res$coef_df_reml$pd), tolerance = 1e-12)
})
