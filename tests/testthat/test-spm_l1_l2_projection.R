test_that("extract_spm_l2_projection aligns runs and drops intercept/constant columns", {
  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = TRUE, include_l2_models = FALSE)
  gpa$glm_settings <- list(spm = list())

  gpa$l2_models <- list(
    models = list(
      l2_main = list(
        name = "l2_main",
        l2_scope = "id_session",
        metadata = data.frame(
          id = c("sub1", "sub1"),
          session = c(1L, 1L),
          run_number = c(1L, 2L),
          stringsAsFactors = FALSE
        ),
        model_matrix = cbind(
          "(Intercept)" = c(1, 1),
          drug = c(1, 0),
          const_term = c(5, 5)
        )
      )
    )
  )
  class(gpa$l2_models) <- c("hi_model_set", "list")

  proj <- fmri.pipeline:::extract_spm_l2_projection(
    gpa = gpa,
    l2_model_name = "l2_main",
    id = "sub1",
    session = 1L,
    run_numbers = c(2L, 1L),
    lg = lgr::get_logger("test")
  )

  expect_true(is.data.frame(proj))
  expect_identical(names(proj), c("run_number", "drug"))
  expect_equal(proj$run_number, c(2L, 1L))
  expect_equal(proj$drug, c(0, 1))
})

test_that("extract_spm_l2_projection supports pooled L1 session mode across sessions", {
  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = TRUE, include_l2_models = FALSE)
  gpa$glm_settings <- list(spm = list())

  gpa$l2_models <- list(
    models = list(
      l2_main = list(
        name = "l2_main",
        l2_scope = "id",
        metadata = data.frame(
          id = rep("sub1", 4L),
          session = c(1L, 1L, 2L, 2L),
          run_number = c(1L, 2L, 1L, 2L),
          stringsAsFactors = FALSE
        ),
        model_matrix = cbind(
          "(Intercept)" = rep(1, 4L),
          drug = c(0, 1, 0, 1),
          session_effect = c(-1, -1, 1, 1)
        )
      )
    )
  )
  class(gpa$l2_models) <- c("hi_model_set", "list")

  proj <- fmri.pipeline:::extract_spm_l2_projection(
    gpa = gpa,
    l2_model_name = "l2_main",
    id = "sub1",
    session = 0L,
    run_numbers = c(4L, 2L, 1L, 3L),
    run_sessions = c(2L, 1L, 1L, 2L),
    source_run_numbers = c(2L, 2L, 1L, 1L),
    l1_session_mode = "pooled",
    lg = lgr::get_logger("test")
  )

  expect_true(is.data.frame(proj))
  expect_identical(names(proj), c("run_number", "drug", "session_effect"))
  expect_equal(proj$run_number, c(4L, 2L, 1L, 3L))
  expect_equal(proj$drug, c(1, 1, 0, 0))
  expect_equal(proj$session_effect, c(1, -1, -1, 1))
})

test_that("resolve_spm_l1_session_mode normalizes aliases", {
  expect_equal(
    fmri.pipeline:::resolve_spm_l1_session_mode(NULL, lg = lgr::get_logger("test")),
    "separate"
  )
  expect_equal(
    fmri.pipeline:::resolve_spm_l1_session_mode("id", lg = lgr::get_logger("test")),
    "pooled"
  )
  expect_equal(
    fmri.pipeline:::resolve_spm_l1_session_mode("session", lg = lgr::get_logger("test")),
    "separate"
  )
  expect_error(
    fmri.pipeline:::resolve_spm_l1_session_mode("bad_mode", lg = lgr::get_logger("test")),
    "Unknown spm\\$l1_session_mode"
  )
})

test_that("resolve_spm_l2_projection_model honors explicit and fallback choices", {
  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = TRUE, include_l2_models = FALSE)

  gpa$l2_models <- list(
    models = list(
      l2_a = list(name = "l2_a"),
      model1 = list(name = "model1")
    )
  )
  class(gpa$l2_models) <- c("hi_model_set", "list")

  gpa$glm_settings <- list(spm = list(l2_projection_model = NULL))
  expect_identical(
    fmri.pipeline:::resolve_spm_l2_projection_model(gpa, l1_model = "model1", lg = lgr::get_logger("test")),
    "model1"
  )

  gpa$glm_settings$spm$l2_projection_model <- "l2_a"
  expect_identical(
    fmri.pipeline:::resolve_spm_l2_projection_model(gpa, l1_model = "model1", lg = lgr::get_logger("test")),
    "l2_a"
  )

  gpa$glm_settings$spm$l2_projection_model <- "none"
  expect_null(
    fmri.pipeline:::resolve_spm_l2_projection_model(gpa, l1_model = "model1", lg = lgr::get_logger("test"))
  )
})

test_that("resolve_spm_l2_projection_interactions handles all/none/explicit", {
  terms <- c("face_emotion", "drug")

  expect_equal(
    fmri.pipeline:::resolve_spm_l2_projection_interactions(TRUE, terms, lg = lgr::get_logger("test")),
    terms
  )
  expect_equal(
    fmri.pipeline:::resolve_spm_l2_projection_interactions("all", terms, lg = lgr::get_logger("test")),
    terms
  )
  expect_equal(
    fmri.pipeline:::resolve_spm_l2_projection_interactions(c("drug"), terms, lg = lgr::get_logger("test")),
    "drug"
  )
  expect_equal(
    fmri.pipeline:::resolve_spm_l2_projection_interactions("none", terms, lg = lgr::get_logger("test")),
    character(0)
  )
  expect_error(
    fmri.pipeline:::resolve_spm_l2_projection_interactions(c("bad_term"), terms, lg = lgr::get_logger("test")),
    "Unknown spm\\$l2_projection_interactions"
  )
})

test_that("resolve_spm_l2_projection_interaction_contrast_modes handles aliases and defaults", {
  expect_equal(
    fmri.pipeline:::resolve_spm_l2_projection_interaction_contrast_modes(
      modes = NULL,
      legacy_unit_contrasts = TRUE,
      lg = lgr::get_logger("test")
    ),
    "pooled"
  )
  expect_equal(
    fmri.pipeline:::resolve_spm_l2_projection_interaction_contrast_modes(
      modes = "cell_means",
      legacy_unit_contrasts = TRUE,
      lg = lgr::get_logger("test")
    ),
    "session_specific"
  )
  expect_equal(
    fmri.pipeline:::resolve_spm_l2_projection_interaction_contrast_modes(
      modes = c("condition_means", "session_differences"),
      legacy_unit_contrasts = FALSE,
      lg = lgr::get_logger("test")
    ),
    c("pooled", "session_differences")
  )
  expect_equal(
    fmri.pipeline:::resolve_spm_l2_projection_interaction_contrast_modes(
      modes = "none",
      legacy_unit_contrasts = TRUE,
      lg = lgr::get_logger("test")
    ),
    character(0)
  )
  expect_error(
    fmri.pipeline:::resolve_spm_l2_projection_interaction_contrast_modes(
      modes = "not_a_mode",
      legacy_unit_contrasts = TRUE,
      lg = lgr::get_logger("test")
    ),
    "Unknown spm\\$l2_projection_interaction_contrast_modes"
  )
})

test_that("spm_l1_model augments contrast matrix with projected interaction terms", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_l1_int_con_")
  dir.create(tmp_dir, recursive = TRUE)
  nifti1 <- file.path(tmp_dir, "run1.nii")
  nifti2 <- file.path(tmp_dir, "run2.nii")
  img <- RNifti::asNifti(array(0, dim = c(2L, 2L, 2L, 10L)))
  RNifti::writeNifti(img, nifti1)
  RNifti::writeNifti(img, nifti2)

  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = TRUE, output_directory = tmp_dir)
  gpa$glm_software <- "spm"
  gpa$glm_settings <- list(spm = list(
    concatenate_runs = TRUE,
    spm_path = tmp_dir,
    l2_projection_model = "l2_main",
    l2_projection_interactions = "face_emotion",
    l2_projection_interaction_unit_contrasts = TRUE
  ))
  gpa$output_locations$spm_l1_directory <- file.path(tmp_dir, "{id}", "ses-{session}", "{l1_model}")
  gpa$parallel$compute_environment <- list(
    global = character(0), fsl = character(0), afni = character(0), spm = character(0), r = character(0)
  )
  gpa$l1_models$models <- list(
    model1 = list(
      name = "model1",
      signals = character(0),
      contrasts = matrix(
        1,
        nrow = 1,
        dimnames = list("base_face", "face_onset")
      )
    )
  )
  class(gpa$l1_models$models$model1) <- c("l1_model_spec", "list")

  gpa$l2_models <- list(
    models = list(
      l2_main = list(
        name = "l2_main",
        l2_scope = "id_session",
        metadata = data.frame(
          id = c("sub1", "sub1"),
          session = c(1L, 1L),
          run_number = c(1L, 2L),
          stringsAsFactors = FALSE
        ),
        model_matrix = cbind(
          "(Intercept)" = c(1, 1),
          face_emotion = c(1, -1)
        )
      )
    )
  )
  class(gpa$l2_models) <- c("hi_model_set", "list")

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  reg_run1 <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "face_onset")
  reg_run2 <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "face_onset")
  reg_concat <- mk_reg(c(1, 3, 11, 13), c(1, 1, 1, 1), c(1, 1, 1, 1), "face_onset")

  d_obj <- list(
    design = array(
      list(reg_run1, reg_run2),
      dim = c(2L, 1L),
      dimnames = list(c("run1", "run2"), "face_reg")
    ),
    design_concat = list(face_reg = reg_concat),
    run_volumes = c(10, 10),
    run_niftis = c(nifti1, nifti2),
    tr = 1.0,
    runs_to_output = c(1L, 2L)
  )
  class(d_obj) <- c("bdm", "list")

  out <- spm_l1_model(
    id = "sub1", session = 1L, model_name = "model1",
    d_obj = d_obj, gpa = gpa,
    run_nifti = c(nifti1, nifti2),
    run_numbers = c(1L, 2L)
  )

  expect_true(is.data.frame(out))
  spec_path <- file.path(out$spm_dir[1], "spm_contrast_spec.rds")
  expect_true(file.exists(spec_path))
  spec <- readRDS(spec_path)
  expect_true("face_onset_x_l2_face_emotion" %in% colnames(spec$contrast_matrix))
  expect_true("proj_int_face_onset_x_l2_face_emotion" %in% rownames(spec$contrast_matrix))
  expect_equal(
    spec$contrast_matrix["proj_int_face_onset_x_l2_face_emotion", "face_onset_x_l2_face_emotion"],
    1
  )
})

test_that("spm_l1_model carries pooled non-concat projection weights into contrast spec", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_l1_proj_main_")
  dir.create(tmp_dir, recursive = TRUE)
  nifti_paths <- file.path(tmp_dir, paste0("run", 1:4, ".nii"))
  img <- RNifti::asNifti(array(0, dim = c(2L, 2L, 2L, 8L)))
  for (pp in nifti_paths) RNifti::writeNifti(img, pp)

  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 4, include_l1_models = TRUE, output_directory = tmp_dir)
  gpa$glm_software <- "spm"
  gpa$glm_settings <- list(spm = list(
    concatenate_runs = FALSE,
    spm_path = tmp_dir,
    l2_projection_model = "l2_main"
  ))
  gpa$output_locations$spm_l1_directory <- file.path(tmp_dir, "{id}", "ses-{session}", "{l1_model}")
  gpa$parallel$compute_environment <- list(
    global = character(0), fsl = character(0), afni = character(0), spm = character(0), r = character(0)
  )
  gpa$l1_models$models <- list(
    model1 = list(
      name = "model1",
      signals = character(0),
      contrasts = matrix(
        1,
        nrow = 1,
        dimnames = list("face_main", "face_onset")
      )
    )
  )
  class(gpa$l1_models$models$model1) <- c("l1_model_spec", "list")

  gpa$l2_models <- list(
    models = list(
      l2_main = list(
        name = "l2_main",
        l2_scope = "id",
        metadata = data.frame(
          id = rep("sub1", 4L),
          session = c(1L, 1L, 2L, 2L),
          run_number = c(1L, 2L, 1L, 2L),
          stringsAsFactors = FALSE
        ),
        model_matrix = cbind(
          "(Intercept)" = rep(1, 4L),
          drug = c(1, 1, 0, 0)
        )
      )
    )
  )
  class(gpa$l2_models) <- c("hi_model_set", "list")

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  reg <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "face_onset")
  d_obj <- list(
    design = array(
      list(reg, reg, reg, reg),
      dim = c(4L, 1L),
      dimnames = list(c("run1", "run2", "run3", "run4"), "face_reg")
    ),
    run_volumes = rep(8L, 4L),
    run_niftis = nifti_paths,
    tr = 1.0,
    runs_to_output = c(1L, 2L, 3L, 4L)
  )
  class(d_obj) <- c("bdm", "list")

  out <- NULL
  expect_warning(
    out <- spm_l1_model(
      id = "sub1", session = 1L, model_name = "model1",
      d_obj = d_obj, gpa = gpa,
      run_nifti = nifti_paths,
      run_numbers = c(1L, 2L, 3L, 4L),
      run_sessions = c(1L, 1L, 2L, 2L),
      source_run_numbers = c(1L, 2L, 1L, 2L),
      l1_session_mode = "pooled"
    ),
    "Skipping projected L2 main-effect regressors in non-concatenated SPM design"
  )

  expect_true(is.data.frame(out))
  spec_path <- file.path(out$spm_dir[1], "spm_contrast_spec.rds")
  expect_true(file.exists(spec_path))
  spec <- readRDS(spec_path)
  expect_identical(spec$projection_main_effect_terms, "drug")
  expect_equal(spec$projection_main_effect_weights$drug, c(1, 1, 0, 0))
})
