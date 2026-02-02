test_that("intercept-only L2/L3 models auto-add diagonal contrast", {
  df <- data.frame(y = rnorm(6))
  lmfit <- lm(y ~ 1, data = df)

  mobj <- list(level = 2L, lmfit = lmfit)
  class(mobj) <- c("list", "hi_model_spec")

  res <- fmri.pipeline:::specify_contrasts(mobj, spec_list = NULL)
  expect_true(isTRUE(res$contrast_spec$diagonal))
  expect_true(is.matrix(res$contrasts))
  expect_equal(dim(res$contrasts), c(1L, 1L))
  expect_equal(res$contrasts[1, 1], 1)
})

test_that("intercept-only models with spec_list warn and keep diagonal only", {
  df <- data.frame(y = rnorm(6))
  lmfit <- lm(y ~ 1, data = df)

  mobj <- list(level = 3L, lmfit = lmfit)
  class(mobj) <- c("list", "hi_model_spec")

  spec_list <- list(
    diagonal = TRUE,
    cond_means = "group"
  )

  res <- fmri.pipeline:::specify_contrasts(mobj, spec_list = spec_list)
  expect_true(isTRUE(res$contrast_spec$diagonal))
  expect_equal(dim(res$contrasts), c(1L, 1L))
})
