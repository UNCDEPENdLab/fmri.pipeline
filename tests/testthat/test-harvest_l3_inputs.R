library(testthat)

test_that("harvest_l3_inputs resolves AFNI requests through the FSL producer and preserves subsets", {
  gpa <- list(
    subject_data = data.frame(
      id = c("sub1", "sub2"),
      session = c(1L, 1L),
      exclude_subject = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    l1_models = structure(
      list(models = list(l1_model1 = list(), l1_model2 = list())),
      class = c("l1_model_set", "list")
    ),
    l2_models = structure(
      list(models = list(
        l2_model1 = list(l1_model_names = c("l1_model1", "l1_model2"))
      )),
      class = c("hi_model_set", "list")
    ),
    l3_models = structure(
      list(models = list(l3_model1 = list())),
      class = c("hi_model_set", "list")
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  capture <- new.env(parent = emptyenv())

  local_mocked_bindings(
    resolve_l2_l3_compatible_pairs = function(gpa, l2_model_names = NULL, l3_model_names = NULL, lg = NULL) {
      capture$l2_model_names <- l2_model_names
      capture$l3_model_names <- l3_model_names
      data.frame(l2_model = "l2_model1", l3_model = "l3_model1", stringsAsFactors = FALSE)
    },
    get_fsl_l3_model_df = function(gpa, model_df, subj_df) {
      capture$model_df <- model_df
      capture$subj_df <- subj_df
      data.frame(
        id = subj_df$id[1],
        session = subj_df$session[1],
        l1_model = model_df$l1_model[1],
        l2_model = model_df$l2_model[1],
        l3_model = model_df$l3_model[1],
        l1_cope_name = "cope1",
        l2_cope_name = "copeA",
        stringsAsFactors = FALSE
      )
    },
    get_feat_l3_inputs = function(gpa, l3_cope_config, lg = NULL) {
      list(
        contrast1 = cbind(
          l3_cope_config,
          data.frame(cope_file = "cope.nii.gz", stringsAsFactors = FALSE)
        )
      )
    },
    .package = "fmri.pipeline"
  )

  harvested <- fmri.pipeline:::harvest_l3_inputs(
    gpa = gpa,
    l3_backend = list(name = "afni", l3_input_provider_backends = "fsl"),
    l1_model_names = "l1_model2",
    l2_model_names = "l2_model1",
    l3_model_names = "l3_model1"
  )

  expect_identical(capture$l2_model_names, "l2_model1")
  expect_identical(capture$l3_model_names, "l3_model1")
  expect_identical(capture$model_df$l1_model, "l1_model2")
  expect_identical(capture$model_df$l2_model, "l2_model1")
  expect_identical(capture$model_df$l3_model, "l3_model1")
  expect_identical(capture$subj_df$id, c("sub1", "sub2"))
  expect_identical(capture$subj_df$session, c(1L, 1L))

  expect_true(is.list(harvested))
  expect_true("l3_model1" %in% names(harvested))
  expect_true("contrast1" %in% names(harvested$l3_model1))
  expect_identical(harvested$l3_model1$contrast1$InputFile, "cope.nii.gz")
  expect_identical(harvested$l3_model1$contrast1$source_backend, "fsl")
})

test_that("harvest_l3_inputs does not duplicate AFNI inputs across L3 contrasts", {
  gpa <- list(
    subject_data = data.frame(
      id = c("sub1", "sub2"),
      session = c(1L, 1L),
      exclude_subject = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    l1_models = structure(
      list(models = list(l1_model1 = list())),
      class = c("l1_model_set", "list")
    ),
    l2_models = structure(
      list(models = list(l2_model1 = list(l1_model_names = "l1_model1"))),
      class = c("hi_model_set", "list")
    ),
    l3_models = structure(
      list(models = list(l3_model1 = list())),
      class = c("hi_model_set", "list")
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  capture <- new.env(parent = emptyenv())

  local_mocked_bindings(
    resolve_l2_l3_compatible_pairs = function(gpa, l2_model_names = NULL, l3_model_names = NULL, lg = NULL) {
      data.frame(l2_model = "l2_model1", l3_model = "l3_model1", stringsAsFactors = FALSE)
    },
    get_fsl_l3_model_df = function(gpa, model_df, subj_df) {
      base <- data.frame(
        id = subj_df$id,
        session = subj_df$session,
        l1_model = model_df$l1_model[1],
        l2_model = model_df$l2_model[1],
        l3_model = model_df$l3_model[1],
        l1_cope_name = "cope1",
        l2_cope_name = "copeA",
        stringsAsFactors = FALSE
      )
      do.call(
        rbind,
        lapply(1:3, function(ii) {
          transform(base, l3_cope_number = ii, l3_cope_name = paste0("l3_contrast", ii))
        })
      )
    },
    get_feat_l3_inputs = function(gpa, l3_cope_config, lg = NULL) {
      capture$l3_cope_config <- l3_cope_config
      list(
        contrast1 = cbind(
          l3_cope_config,
          data.frame(cope_file = paste0("cope", seq_len(nrow(l3_cope_config)), ".nii.gz"), stringsAsFactors = FALSE)
        )
      )
    },
    .package = "fmri.pipeline"
  )

  harvested <- fmri.pipeline:::harvest_l3_inputs(
    gpa = gpa,
    l3_backend = list(name = "afni", l3_input_provider_backends = "fsl"),
    l1_model_names = "l1_model1",
    l2_model_names = "l2_model1",
    l3_model_names = "l3_model1"
  )

  expect_equal(nrow(capture$l3_cope_config), 2L)
  expect_false("l3_cope_number" %in% names(capture$l3_cope_config))
  expect_false("l3_cope_name" %in% names(capture$l3_cope_config))
  expect_equal(nrow(harvested$l3_model1$contrast1), 2L)
  expect_equal(nrow(unique(harvested$l3_model1$contrast1[, c("id", "session")])), 2L)
})

test_that("harvest_l3_inputs returns NULL and warns for unsupported producer backend", {
  gpa <- list(
    subject_data = data.frame(
      id = "sub1", session = 1L, exclude_subject = FALSE,
      stringsAsFactors = FALSE
    ),
    l1_models = structure(
      list(models = list(l1_model1 = list())),
      class = c("l1_model_set", "list")
    ),
    l2_models = structure(
      list(models = list(l2_model1 = list(l1_model_names = "l1_model1"))),
      class = c("hi_model_set", "list")
    ),
    l3_models = structure(
      list(models = list(l3_model1 = list())),
      class = c("hi_model_set", "list")
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  # Provide a backend whose producer is "spm" — not handled by the switch
  expect_error(
    fmri.pipeline:::harvest_l3_inputs(
      gpa = gpa,
      l3_backend = list(name = "custom", l3_input_provider_backends = "spm"),
      l1_model_names = "l1_model1",
      l2_model_names = "l2_model1",
      l3_model_names = "l3_model1"
    ),
    "does not support producer backend"
  )
})
