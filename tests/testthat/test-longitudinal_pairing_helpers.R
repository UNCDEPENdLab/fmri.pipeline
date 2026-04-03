library(testthat)

make_signature_test_gpa <- function() {
  gpa <- create_mock_gpa(
    n_subjects = 3,
    n_runs = 2,
    include_l2_models = TRUE,
    include_l3_models = TRUE
  )

  gpa$l2_models$models$l2_model2 <- gpa$l2_models$models$l2_model1
  gpa$l2_models$models$l2_model2$name <- "l2_model2"

  gpa$l3_models$models$l3_model2 <- gpa$l3_models$models$l3_model1
  gpa$l3_models$models$l3_model2$name <- "l3_model2"
  gpa$l3_models$models$l3_model3 <- gpa$l3_models$models$l3_model1
  gpa$l3_models$models$l3_model3$name <- "l3_model3"

  class(gpa$l2_models) <- c("list", "hi_model_set")
  class(gpa$l3_models) <- c("list", "hi_model_set")

  gpa
}

test_that("normalize_longitudinal_model_signatures sets defaults", {
  gpa <- make_signature_test_gpa()

  gpa2 <- normalize_longitudinal_model_signatures(gpa)
  l2_scopes <- vapply(gpa2$l2_models$models, function(mm) mm$l2_scope, character(1))
  l3_modes <- vapply(gpa2$l3_models$models, function(mm) mm$l3_input_mode, character(1))

  expect_true(all(l2_scopes == "id_session"))
  expect_true(all(l3_modes == "per_session"))
})

test_that("normalize_longitudinal_model_signatures rejects invalid L2 scope", {
  gpa <- make_signature_test_gpa()
  gpa$l2_models$models$l2_model1$l2_scope <- "bad_scope"

  expect_error(
    normalize_longitudinal_model_signatures(gpa),
    "Unknown l2_scope"
  )
})

test_that("normalize_longitudinal_model_signatures rejects invalid L3 input mode", {
  gpa <- make_signature_test_gpa()
  gpa$l3_models$models$l3_model1$l3_input_mode <- "bad_mode"

  expect_error(
    normalize_longitudinal_model_signatures(gpa),
    "Unknown l3_input_mode"
  )
})

test_that("normalize_longitudinal_model_signatures rejects legacy separate_sessions mode", {
  gpa <- make_signature_test_gpa()
  gpa$l3_models$models$l3_model1$l3_input_mode <- "separate_sessions"

  expect_error(
    normalize_longitudinal_model_signatures(gpa),
    "Unknown l3_input_mode"
  )
})

test_that("resolve_l2_l3_compatible_pairs enforces signature compatibility", {
  gpa <- make_signature_test_gpa()
  gpa$l2_models$models$l2_model1$l2_scope <- "id_session"
  gpa$l2_models$models$l2_model2$l2_scope <- "id"
  gpa$l3_models$models$l3_model1$l3_input_mode <- "per_session"
  gpa$l3_models$models$l3_model2$l3_input_mode <- "subject_rows"
  gpa$l3_models$models$l3_model3$l3_input_mode <- "pooled_sessions_subject_ev"

  pairs <- resolve_l2_l3_compatible_pairs(gpa)
  expect_equal(nrow(pairs), 3L)
  expect_true(any(pairs$l2_model == "l2_model1" & pairs$l3_model == "l3_model1"))
  expect_true(any(pairs$l2_model == "l2_model1" & pairs$l3_model == "l3_model3"))
  expect_true(any(pairs$l2_model == "l2_model2" & pairs$l3_model == "l3_model2"))
})

test_that("resolve_l2_l3_compatible_pairs respects model subset filters", {
  gpa <- make_signature_test_gpa()
  gpa$l2_models$models$l2_model1$l2_scope <- "id_session"
  gpa$l2_models$models$l2_model2$l2_scope <- "id"
  gpa$l3_models$models$l3_model1$l3_input_mode <- "per_session"
  gpa$l3_models$models$l3_model2$l3_input_mode <- "subject_rows"
  gpa$l3_models$models$l3_model3$l3_input_mode <- "pooled_sessions_subject_ev"

  pairs <- resolve_l2_l3_compatible_pairs(
    gpa,
    l2_model_names = "l2_model2",
    l3_model_names = c("l3_model1", "l3_model2")
  )
  expect_equal(nrow(pairs), 1L)
  expect_equal(unique(pairs$l2_model), "l2_model2")
  expect_equal(unique(pairs$l3_model), "l3_model2")
})

test_that("enumerate_l2_l3_signature_pairs reports incompatibility reasons", {
  gpa <- make_signature_test_gpa()
  gpa$l2_models$models$l2_model1$l2_scope <- "id_session"
  gpa$l2_models$models$l2_model2$l2_scope <- "id"
  gpa$l3_models$models$l3_model1$l3_input_mode <- "per_session"
  gpa$l3_models$models$l3_model2$l3_input_mode <- "subject_rows"
  gpa$l3_models$models$l3_model3$l3_input_mode <- "pooled_sessions_subject_ev"

  all_pairs <- enumerate_l2_l3_signature_pairs(gpa)
  expect_equal(nrow(all_pairs), 6L)
  expect_true(all(c("compatible", "reason", "pair_id") %in% names(all_pairs)))

  bad_row <- all_pairs[all_pairs$l2_model == "l2_model1" & all_pairs$l3_model == "l3_model2", , drop = FALSE]
  expect_false(bad_row$compatible[[1]])
  expect_match(bad_row$reason[[1]], "requires l2_scope='id'")
})

test_that("format_l2_l3_incompatibilities summarizes pair diagnostics", {
  gpa <- make_signature_test_gpa()
  gpa$l2_models$models$l2_model1$l2_scope <- "id"
  gpa$l2_models$models$l2_model2$l2_scope <- "id"
  gpa$l3_models$models$l3_model1$l3_input_mode <- "per_session"
  gpa$l3_models$models$l3_model2$l3_input_mode <- "subject_rows"
  gpa$l3_models$models$l3_model3$l3_input_mode <- "pooled_sessions_subject_ev"

  all_pairs <- enumerate_l2_l3_signature_pairs(gpa)
  txt <- format_l2_l3_incompatibilities(all_pairs, max_items = 2L)
  expect_match(txt, "requires l2_scope='id_session'")
  expect_match(txt, "\\+2 more incompatible pair")
})
