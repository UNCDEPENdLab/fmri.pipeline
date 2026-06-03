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

test_that("model-level backend fields do not collide with contrast validation", {
  df <- data.frame(y = rnorm(6), session_label = factor(rep(c("A", "B"), 3)))
  lmfit <- lm(y ~ session_label, data = df)
  mobj <- list(level = 3L, regressors = names(coef(lmfit)), lmfit = lmfit)
  class(mobj) <- c("list", "hi_model_spec")

  backend_spec <- list(
    name = "l3_pooled_subject_ev",
    level = 3L,
    l3_input_mode = "pooled_sessions_subject_ev",
    execution_backend = "fsl",
    producer_backend = "fsl",
    contrasts = list(
      diagonal = TRUE
    )
  )

  expect_no_error(
    res <- fmri.pipeline:::specify_contrasts(mobj, signals = list(), spec_list = backend_spec)
  )
  expect_true(isTRUE(res$contrast_spec$diagonal))

  bad_nested_spec <- backend_spec
  bad_nested_spec$execution_backend <- NULL
  bad_nested_spec$producer_backend <- NULL
  bad_nested_spec$contrasts$execution_backend <- "fsl"

  expect_error(
    fmri.pipeline:::specify_contrasts(mobj, signals = list(), spec_list = bad_nested_spec),
    "Unknown contrast field"
  )
})
