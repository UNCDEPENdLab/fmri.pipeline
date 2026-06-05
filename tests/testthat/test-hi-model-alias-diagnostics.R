test_that("create_new_hi_model stops before contrasts for aliased higher-level terms", {
  subject_data <- data.frame(
    id = rep(sprintf("sub%d", 1:4), each = 3),
    session = rep(1:3, times = 4),
    session_label = rep(c("analgesic", "placebo", "opioid"), times = 4),
    stringsAsFactors = FALSE
  )
  session_days <- c(analgesic = 0, placebo = 7, opioid = 14)
  subject_data$days_since_baseline <- as.numeric(session_days[subject_data$session_label])

  spec <- list(
    name = "aliased_l3",
    model_formula = "~ session_label + days_since_baseline",
    l3_input_mode = "3dlmer",
    random_effects = "(1 | Subj)",
    reference_level = c(session_label = "analgesic"),
    covariate_transform = NULL
  )

  err <- expect_error(
    fmri.pipeline:::create_new_hi_model(
      data = subject_data,
      level = 3L,
      spec_list = spec
    ),
    "aliased coefficients"
  )

  msg <- conditionMessage(err)
  expect_match(msg, "L3 model 'aliased_l3'")
  expect_match(msg, "days_since_baseline")
  expect_match(msg, "session_label")
  expect_match(msg, "rank deficient")
  expect_match(msg, "one value within each session_label level", fixed = TRUE)
  expect_match(msg, "Please remove or recode redundant predictors and respecify the model", fixed = TRUE)
})

test_that("format_hi_model_alias_message explains design-matrix dependencies", {
  data <- data.frame(
    id = rep(sprintf("sub%d", 1:4), each = 3),
    session = rep(1:3, times = 4),
    session_label = rep(c("analgesic", "placebo", "opioid"), times = 4),
    stringsAsFactors = FALSE
  )
  data$days_since_baseline <- as.numeric(c(0, 7, 14)[data$session])

  mobj <- list(level = 3L, name = "aliased_l3")
  class(mobj) <- c("list", "hi_model_spec")
  mobj$reference_level <- c(session_label = "analgesic")

  mobj <- fmri.pipeline:::mobj_fit_lm(
    mobj = mobj,
    model_formula = "~ session_label + days_since_baseline",
    data = data,
    id_cols = c("id", "session")
  )

  msg <- fmri.pipeline:::format_hi_model_alias_message(mobj)
  expect_match(msg, "Perfect redundancy detected")
  expect_match(msg, "days_since_baseline")
  expect_match(msg, "session_label")
})
