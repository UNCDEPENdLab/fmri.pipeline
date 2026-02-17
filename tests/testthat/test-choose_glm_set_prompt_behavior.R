test_that("choose_glm_set skips confirmation prompt when model names are explicit", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$multi_run <- TRUE
  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")

  choose_glm_set_local <- fmri.pipeline:::choose_glm_set
  override_env <- new.env(parent = environment(choose_glm_set_local))
  override_env$interactive <- function() TRUE
  override_env$menu <- function(...) stop("menu should not be called")
  environment(choose_glm_set_local) <- override_env

  result <- choose_glm_set_local(
    gpa,
    l1_model_names = "model1",
    l2_model_names = "l2_model1",
    l3_model_names = "l3_model1",
    lg = lg
  )

  expect_equal(result$l1_model_names, "model1")
  expect_equal(result$l2_model_names, "l2_model1")
  expect_equal(result$l3_model_names, "l3_model1")
})
