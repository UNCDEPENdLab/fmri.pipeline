library(testthat)

# Helper to create a minimal hi_model_spec for respecify_l3_model testing
make_l3_mobj <- function(formula_str = "~ 1", l3_input_mode = "per_session",
                         subject_data = NULL) {
  if (is.null(subject_data)) {
    subject_data <- data.frame(
      id = rep(c("sub1", "sub2", "sub3"), each = 2),
      session = rep(c(1L, 2L), 3),
      age = c(25, 26, 30, 31, 40, 41),
      stringsAsFactors = FALSE
    )
  }

  # Build initial model object using mobj_fit_lm
  mobj <- list(level = 3L, l3_input_mode = l3_input_mode)
  class(mobj) <- c("list", "hi_model_spec")

  # add dummy DV
  subject_data$dummy <- seq_len(nrow(subject_data))
  id_cols <- c("id", "session")

  mobj <- mobj_fit_lm(
    mobj = mobj,
    model_formula = formula_str,
    data = subject_data,
    id_cols = id_cols
  )

  # Add diagonal contrasts by default
  mobj$contrast_spec <- list(diagonal = TRUE)
  mobj$contrast_list <- list()
  mobj <- get_contrasts_from_spec(mobj)

  return(mobj)
}


test_that("respecify_l3_model auto-injects id EVs for pooled_sessions_subject_ev", {
  mobj <- make_l3_mobj(
    formula_str = "~ age",
    l3_input_mode = "pooled_sessions_subject_ev"
  )

  # Verify id is NOT in model variables before respecification
  expect_false("id" %in% mobj$model_variables)

  # Create new_data representing available subjects/sessions
  new_data <- data.frame(
    id = rep(c("sub1", "sub2", "sub3"), each = 2),
    session = rep(c(1L, 2L), 3),
    stringsAsFactors = FALSE
  )

  result <- respecify_l3_model(mobj, new_data)

  # Verify that id columns are present in the model matrix
  mm_cols <- colnames(result$model_matrix)
  expect_true(any(grepl("^id", mm_cols)),
    info = "Model matrix should contain id indicator columns after auto-injection"
  )

  # Verify there is no global intercept (using ~ 0 + id + ...)
  expect_false("(Intercept)" %in% mm_cols,
    info = "Auto-injected subject EVs should suppress the global intercept"
  )

  # Verify there are per-subject columns (one per subject)
  id_cols <- mm_cols[grepl("^id", mm_cols)]
  expect_equal(length(id_cols), 3L,
    info = "Should have one indicator EV per subject"
  )

  # Verify the age covariate is still present
  expect_true("age" %in% mm_cols,
    info = "Original covariates should be preserved alongside injected subject EVs"
  )
})


test_that("respecify_l3_model does NOT inject id EVs for per_session mode", {
  mobj <- make_l3_mobj(
    formula_str = "~ 1",
    l3_input_mode = "per_session"
  )

  # Create new_data for a single session
  new_data <- data.frame(
    id = c("sub1", "sub2", "sub3"),
    session = rep(1L, 3),
    stringsAsFactors = FALSE
  )

  result <- respecify_l3_model(mobj, new_data)

  # Should NOT have id columns — just intercept
  mm_cols <- colnames(result$model_matrix)
  expect_false(any(grepl("^id", mm_cols)),
    info = "per_session mode should not inject id EVs"
  )
  expect_true("(Intercept)" %in% mm_cols)
})


test_that("respecify_l3_model skips injection if id already in formula", {
  # Build a model that already has id in the formula
  subject_data <- data.frame(
    id = factor(rep(c("sub1", "sub2", "sub3"), each = 2)),
    session = rep(c(1L, 2L), 3),
    age = c(25, 26, 30, 31, 40, 41),
    stringsAsFactors = FALSE
  )

  mobj <- make_l3_mobj(
    formula_str = "~ 0 + id + age",
    l3_input_mode = "pooled_sessions_subject_ev",
    subject_data = subject_data
  )

  new_data <- data.frame(
    id = factor(rep(c("sub1", "sub2", "sub3"), each = 2)),
    session = rep(c(1L, 2L), 3),
    stringsAsFactors = FALSE
  )

  result <- respecify_l3_model(mobj, new_data)

  # Should still have id columns but no double-injection
  mm_cols <- colnames(result$model_matrix)
  id_cols <- mm_cols[grepl("^id", mm_cols)]
  expect_equal(length(id_cols), 3L,
    info = "Should have exactly 3 subject indicator columns (no duplication)"
  )
})


test_that("respecify_l3_model handles intercept-only formula with pooled_sessions_subject_ev", {
  mobj <- make_l3_mobj(
    formula_str = "~ 1",
    l3_input_mode = "pooled_sessions_subject_ev"
  )

  new_data <- data.frame(
    id = rep(c("sub1", "sub2", "sub3"), each = 2),
    session = rep(c(1L, 2L), 3),
    stringsAsFactors = FALSE
  )

  result <- respecify_l3_model(mobj, new_data)

  # Should have per-subject EVs replacing the intercept
  mm_cols <- colnames(result$model_matrix)
  expect_false("(Intercept)" %in% mm_cols,
    info = "Global intercept should be replaced by per-subject intercepts"
  )
  id_cols <- mm_cols[grepl("^id", mm_cols)]
  expect_equal(length(id_cols), 3L)
})
