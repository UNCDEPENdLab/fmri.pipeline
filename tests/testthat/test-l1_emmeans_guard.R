test_that("L1 spec drops emmeans fields when no model is available", {
  mobj <- list(
    level = 1L,
    regressors = c("reg1"),
    signals = c("sig1")
  )
  class(mobj) <- c("list", "l1_model_spec")

  spec_list <- list(
    diagonal = TRUE,
    cond_means = "cond",
    pairwise_diffs = "cond",
    cell_means = TRUE,
    overall_response = TRUE,
    weights = "cells"
  )

  out <- fmri.pipeline:::specify_contrasts(mobj, signals = list(sig1 = list()), spec_list = spec_list)

  expect_true(out$contrast_spec$diagonal)
  expect_length(out$contrast_spec$cond_means, 0)
  expect_length(out$contrast_spec$pairwise_diffs, 0)
  expect_false(out$contrast_spec$cell_means)
  expect_false(out$contrast_spec$overall_response)
  expect_true(is.matrix(out$contrasts))
  expect_equal(ncol(out$contrasts), 1L)
})
