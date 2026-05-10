test_that("extract_fsl_betas resolves per-cope L2 cope-file outputs", {
  root <- tempfile("extract_l2_cope_files_")
  dir.create(root, recursive = TRUE)

  mask_file <- file.path(root, "mask.nii.gz")
  RNifti::writeNifti(array(1, dim = c(2, 2, 2)), mask_file)

  make_l2_row <- function(cope_number, cope_name, value) {
    feat_dir <- file.path(root, paste0("L1m-model1_l1c-", cope_number, ".gfeat"))
    stats_dir <- file.path(feat_dir, "cope1.feat", "stats")
    dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
    RNifti::writeNifti(array(value, dim = c(2, 2, 2)), file.path(stats_dir, "cope1.nii.gz"))

    row <- data.frame(
      id = "sub1",
      session = 1L,
      l1_model = "model1",
      l1_cope_number = cope_number,
      l1_cope_name = cope_name,
      l2_model = "l2_model1",
      l2_scope = "id_session",
      l2_input_mode = "cope_files",
      l2_passthrough = FALSE,
      n_l2_copes = 1L,
      n_input_files = 2L,
      passthrough_cope_file = NA_character_,
      feat_dir = feat_dir,
      feat_dir_exists = TRUE,
      feat_complete = TRUE,
      stringsAsFactors = FALSE
    )
    row$cope_list <- list(data.frame(
      id = "sub1",
      session = 1L,
      l2_cope_number = 1L,
      l2_cope_name = "intercept",
      stringsAsFactors = FALSE
    ))
    row
  }

  l2_setup <- dplyr::bind_rows(
    make_l2_row(1L, "copeA", 2),
    make_l2_row(2L, "copeB", 4)
  )

  gpa <- list(
    lgr_threshold = "warn",
    subject_data = data.frame(id = "sub1", session = 1L, exclude_subject = FALSE),
    l2_model_setup = structure(list(fsl = l2_setup), class = c("l2_setup", "list"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  res <- fmri.pipeline:::extract_fsl_betas(
    gpa,
    extract = data.frame(l1_model = "model1", l2_model = "l2_model1"),
    level = 2L,
    what = "cope",
    mask_file = mask_file,
    ncores = 1L,
    lg = lgr::get_logger("test/extract_l2_cope_files")
  )

  expect_equal(nrow(res), 2L)
  expect_true(all(c("l1_cope_number", "l1_cope_name", "l2_cope_number", "l2_cope_name") %in% names(res)))
  expect_equal(sort(res$l1_cope_name), c("copeA", "copeB"))
  expect_equal(res$value[match(c("copeA", "copeB"), res$l1_cope_name)], c(2, 4))
  expect_true(all(grepl(file.path("cope1.feat", "stats", "cope1.nii.gz"), res$img, fixed = TRUE)))

  explicit_res <- fmri.pipeline:::extract_fsl_betas(
    gpa,
    extract = data.frame(
      l1_model = "model1",
      l1_cope_number = 2L,
      l1_cope_name = "copeB",
      l2_model = "l2_model1"
    ),
    level = 2L,
    what = "cope",
    mask_file = mask_file,
    ncores = 1L,
    lg = lgr::get_logger("test/extract_l2_cope_files_explicit")
  )
  expect_equal(nrow(explicit_res), 1L)
  expect_equal(explicit_res$l1_cope_name[1], "copeB")
  expect_equal(explicit_res$value[1], 4)
})

test_that("extract_fsl_betas resolves one-run L2 pass-through cope files", {
  root <- tempfile("extract_l2_passthrough_")
  dir.create(root, recursive = TRUE)

  mask_file <- file.path(root, "mask.nii.gz")
  RNifti::writeNifti(array(1, dim = c(2, 2, 2)), mask_file)

  passthrough_cope <- file.path(root, "run1.feat", "stats", "cope1.nii.gz")
  dir.create(dirname(passthrough_cope), recursive = TRUE, showWarnings = FALSE)
  RNifti::writeNifti(array(7, dim = c(2, 2, 2)), passthrough_cope)

  l2_setup <- data.frame(
    id = "sub1",
    session = 1L,
    l1_model = "model1",
    l1_cope_number = 1L,
    l1_cope_name = "copeA",
    l2_model = "l2_model1",
    l2_scope = "id_session",
    l2_input_mode = "l1_cope_file_passthrough",
    l2_passthrough = TRUE,
    n_l2_copes = 1L,
    n_input_files = 1L,
    passthrough_cope_file = passthrough_cope,
    feat_dir = file.path(root, "run1.feat"),
    feat_dir_exists = TRUE,
    feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  l2_setup$cope_list <- list(data.frame(
    id = "sub1",
    session = 1L,
    l2_cope_number = 1L,
    l2_cope_name = "intercept",
    stringsAsFactors = FALSE
  ))

  gpa <- list(
    lgr_threshold = "warn",
    subject_data = data.frame(id = "sub1", session = 1L, exclude_subject = FALSE),
    l2_model_setup = structure(list(fsl = l2_setup), class = c("l2_setup", "list"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  res <- fmri.pipeline:::extract_fsl_betas(
    gpa,
    extract = data.frame(l1_model = "model1", l2_model = "l2_model1"),
    level = 2L,
    what = "cope",
    mask_file = mask_file,
    ncores = 1L,
    lg = lgr::get_logger("test/extract_l2_passthrough")
  )

  expect_equal(nrow(res), 1L)
  expect_equal(res$img[1], passthrough_cope)
  expect_equal(res$value[1], 7)
  expect_equal(res$l2_input_mode[1], "l1_cope_file_passthrough")
})
