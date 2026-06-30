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

test_that("mixed_by rstanarm engine supports emmeans and emtrends pairs", {
  testthat::skip_if_not_installed("rstanarm")
  testthat::skip_if_not_installed("bayestestR")

  set.seed(203)
  df <- data.frame(
    id = factor(rep(rep(1:4, each = 3), times = 2)),
    split = rep(c("a", "b"), each = 12),
    g = factor(rep(c("low", "high", "low"), times = 8)),
    x = rnorm(24)
  )
  df$y <- 0.5 + 0.35 * df$x + 0.3 * (df$g == "high") +
    0.2 * df$x * (df$g == "high") +
    rep(rnorm(8, sd = 0.2), each = 3) + rnorm(24, sd = 0.4)

  res <- suppressMessages(suppressWarnings(mixed_by(
    df,
    outcomes = "y",
    rhs_model_formulae = list(m1 = ~ g * x + (1 | id)),
    split_on = "split",
    engine = "rstanarm",
    engine_args = list(chains = 1, iter = 100, warmup = 50, seed = 203, refresh = 0),
    calculate = character(0),
    emmeans_spec = list(
      ep = list(outcome = "y", model_name = "m1", specs = ~ g)
    ),
    emmeans_pairs = TRUE,
    emmeans_pairs_args = list(adjust = "tukey"),
    emtrends_spec = list(
      tp = list(outcome = "y", model_name = "m1", specs = ~ g, var = "x")
    ),
    emtrends_pairs = TRUE,
    emtrends_pairs_args = list(adjust = "tukey")
  )))

  expect_true(data.table::is.data.table(res$emmeans_list$ep))
  expect_true(data.table::is.data.table(res$emmeans_pairs_list$ep))
  expect_true(data.table::is.data.table(res$emtrends_list$tp))
  expect_true(data.table::is.data.table(res$emtrends_pairs_list$tp))
  expect_equal(nrow(res$emmeans_list$ep), 4L)
  expect_equal(nrow(res$emmeans_pairs_list$ep), 2L)
  expect_equal(nrow(res$emtrends_list$tp), 4L)
  expect_equal(nrow(res$emtrends_pairs_list$tp), 2L)
  expect_true(all(c("g", "estimate", "lower.HPD", "upper.HPD") %in% names(res$emmeans_list$ep)))
  expect_true(all(c("contrast", "estimate", "lower.HPD", "upper.HPD") %in% names(res$emmeans_pairs_list$ep)))
  expect_true(all(c("g", "x.trend", "lower.HPD", "upper.HPD") %in% names(res$emtrends_list$tp)))
  expect_true(all(c("contrast", "estimate", "lower.HPD", "upper.HPD") %in% names(res$emtrends_pairs_list$tp)))
})

test_that("mixed_by reapplies pair args when specs returns an emm_list", {
  set.seed(204)
  df <- data.frame(
    id = factor(rep(1:8, each = 4)),
    split = rep(c("a", "b"), each = 16),
    g = factor(rep(c("low", "high"), times = 16))
  )
  df$y <- 0.4 + 0.7 * (df$g == "high") +
    rep(rnorm(8, sd = 0.15), each = 4) + rnorm(32, sd = 0.2)

  res <- suppressMessages(mixed_by(
    df,
    outcomes = "y",
    rhs_model_formulae = list(m1 = ~ g + (1 | id)),
    split_on = "split",
    calculate = character(0),
    emmeans_spec = list(
      ep = list(outcome = "y", model_name = "m1", specs = pairwise ~ g)
    ),
    emmeans_pairs = TRUE,
    emmeans_pairs_args = list(reverse = TRUE)
  ))

  expect_equal(nrow(res$emmeans_list$ep), 4L)
  expect_equal(nrow(res$emmeans_pairs_list$ep), 2L)
  expect_true(all(res$emmeans_pairs_list$ep$contrast == "low - high"))
})
