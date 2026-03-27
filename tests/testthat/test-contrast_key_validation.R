test_that("contrast key aliases are mapped and unknown keys error", {
  lg <- lgr::get_logger("test-contrast-key-validation")
  lg$set_threshold("fatal")

  mobj_diag <- list(level = 1L, regressors = "stim")
  class(mobj_diag) <- c("list", "l1_model_spec")

  alias_spec_diag <- list(
    name = "basic",
    signals = list("stim"),
    contrasts = list(
      include_diagonal = TRUE
    )
  )

  res <- fmri.pipeline:::specify_contrasts(mobj_diag, signals = list(), spec_list = alias_spec_diag)
  expect_true(isTRUE(res$contrast_spec$diagonal))

  df <- data.frame(y = rnorm(6), cond = factor(rep(c("A", "B"), 3)))
  lmfit <- lm(y ~ cond, data = df)
  mobj_pair <- list(level = 2L, regressors = names(coef(lmfit)), lmfit = lmfit)
  class(mobj_pair) <- c("list", "hi_model_spec")

  alias_spec_pair <- list(
    name = "basic",
    contrasts = list(
      pairwise_diff = "cond"
    )
  )

  res_pair <- fmri.pipeline:::specify_contrasts(mobj_pair, signals = list(), spec_list = alias_spec_pair)
  expect_equal(res_pair$contrast_spec$pairwise_diffs, "cond")

  bad_spec <- list(
    name = "basic",
    signals = list("stim"),
    contrasts = list(
      diagonal = TRUE,
      pairwise_diffz = "condA"
    )
  )

  expect_error(
    fmri.pipeline:::specify_contrasts(mobj_diag, signals = list(), spec_list = bad_spec),
    "Unknown contrast field"
  )
})
