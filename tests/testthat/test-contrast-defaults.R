test_that("get_contrasts_from_spec warns and adds diagonal when no contrasts specified", {
  data <- data.frame(
    id = c("sub1", "sub2", "sub3"),
    session = c("ses1", "ses1", "ses1"),
    run_number = c(1, 2, 3),
    stringsAsFactors = FALSE
  )

  mobj <- list(level = 2L)
  class(mobj) <- c("list", "hi_model_spec")
  mobj <- fmri.pipeline:::mobj_fit_lm(
    mobj = mobj,
    model_formula = "~1",
    data = data,
    id_cols = c("id", "session", "run_number")
  )

  mobj$contrast_spec <- list(
    diagonal = FALSE,
    cond_means = character(0),
    pairwise_diffs = character(0),
    cell_means = FALSE,
    overall_response = FALSE,
    simple_slopes = list(),
    weights = "cells",
    cat_vars = character(0),
    regressors = mobj$regressors,
    delete = NULL
  )
  mobj$contrast_list <- list()

  res <- NULL
  expect_warning(
    res <- fmri.pipeline:::get_contrasts_from_spec(mobj),
    "No valid contrasts were specified; adding diagonal contrasts by default."
  )

  expect_true(inherits(res, "hi_model_spec"))
  expect_equal(nrow(res$contrasts), 1L)
  expect_equal(rownames(res$contrasts), "EV_Intercept")
  expect_equal(as.numeric(res$contrasts[1, 1]), 1)
})

test_that("fsl_generate_fsf_ev_syntax sanitizes intercept EV titles", {
  dmat <- matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  inputs <- c("feat1.feat", "feat2.feat")

  syntax <- fmri.pipeline::fsl_generate_fsf_ev_syntax(inputs, dmat)

  expect_true(any(grepl("set fmri\\(evtitle1\\) \"Intercept\"", syntax)))
  expect_false(any(grepl("set fmri\\(evtitle1\\) \"\\(Intercept\\)\"", syntax)))
})

test_that("malformed delete specs are ignored without dropping contrasts", {
  data <- data.frame(
    id = c("sub1", "sub2", "sub3"),
    session = c("ses1", "ses1", "ses1"),
    run_number = c(1, 2, 3),
    stringsAsFactors = FALSE
  )

  base_mobj <- list(level = 2L)
  class(base_mobj) <- c("list", "hi_model_spec")
  base_mobj <- fmri.pipeline:::mobj_fit_lm(
    mobj = base_mobj,
    model_formula = "~1",
    data = data,
    id_cols = c("id", "session", "run_number")
  )

  mobj_bad_delete <- base_mobj
  mobj_bad_delete$contrast_spec <- list(
    diagonal = FALSE,
    cond_means = character(0),
    pairwise_diffs = character(0),
    cell_means = FALSE,
    overall_response = FALSE,
    simple_slopes = list(),
    weights = "cells",
    cat_vars = character(0),
    regressors = base_mobj$regressors,
    delete = 1
  )
  mobj_bad_delete$contrast_list <- list()

  res_bad_delete <- NULL
  expect_warning(
    expect_warning(
      res_bad_delete <- fmri.pipeline:::get_contrasts_from_spec(mobj_bad_delete),
      "No valid contrasts were specified; adding diagonal contrasts by default."
    ),
    "contrast_spec\\$delete must be a character vector"
  )
  expect_equal(nrow(res_bad_delete$contrasts), 1L)
  expect_equal(rownames(res_bad_delete$contrasts), "EV_Intercept")

  mobj_empty_delete <- base_mobj
  mobj_empty_delete$contrast_spec <- list(
    diagonal = FALSE,
    cond_means = character(0),
    pairwise_diffs = character(0),
    cell_means = FALSE,
    overall_response = FALSE,
    simple_slopes = list(),
    weights = "cells",
    cat_vars = character(0),
    regressors = base_mobj$regressors,
    delete = character(0)
  )
  mobj_empty_delete$contrast_list <- list()

  res_empty_delete <- NULL
  expect_warning(
    res_empty_delete <- fmri.pipeline:::get_contrasts_from_spec(mobj_empty_delete),
    "No valid contrasts were specified; adding diagonal contrasts by default."
  )
  expect_equal(nrow(res_empty_delete$contrasts), 1L)
  expect_equal(rownames(res_empty_delete$contrasts), "EV_Intercept")
})
