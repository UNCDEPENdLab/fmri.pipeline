library(testthat)

test_that("afni_3dlmer_setup builds 3dLMEr tables from refit model covariates", {
  tmp_dir <- tempfile("afni_3dlmer_setup_")
  dir.create(tmp_dir, recursive = TRUE)

  mask_file <- file.path(tmp_dir, "mask.nii.gz")
  file.create(mask_file)

  harvested <- list(
    l3_model1 = list(
      contrast1 = data.frame(
        id = c("sub1", "sub1", "sub2", "sub2"),
        session = c(1L, 2L, 1L, 2L),
        l1_model = "l1_model1",
        l2_model = "l2_model1",
        l3_model = "l3_model1",
        l1_cope_name = "cope1",
        l2_cope_name = "copeA",
        InputFile = file.path(tmp_dir, paste0("cope", 1:4, ".nii.gz")),
        stringsAsFactors = FALSE
      )
    )
  )

  l3_obj <- list(
    l3_input_mode = "3dlmer",
    model_variables = c("age", "group"),
    lmer_formula = "age + group + (1 | Subj)",
    lmer_mask = NULL,
    lmer_njobs = NULL
  )
  class(l3_obj) <- c("hi_model_spec", "list")

  gpa <- list(
    output_directory = tmp_dir,
    output_locations = list(
      afni_3dlmer_directory = file.path(
        "{gpa$output_directory}", "afni_3dlmer",
        "L1m-{l1_model}", "l1c-{l1_cope_name}",
        "L2m-{l2_model}_l2c-{l2_cope_name}",
        "L3m-{l3_model}"
      )
    ),
    parallel = list(afni = list(l3_lmer_njobs = 2L)),
    mask_file = mask_file,
    l3_models = structure(list(models = list(l3_model1 = l3_obj)), class = c("hi_model_set", "list"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  lg <- lgr::get_logger("test/afni_3dlmer_setup")
  harvest_call <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    harvest_l3_inputs = function(gpa, l3_backend, l1_model_names = NULL, l2_model_names = NULL,
                                 l3_model_names = NULL, lg = NULL) {
      harvest_call$backend <- fmri.pipeline:::backend_name(l3_backend)
      harvest_call$l1_model_names <- l1_model_names
      harvest_call$l2_model_names <- l2_model_names
      harvest_call$l3_model_names <- l3_model_names
      harvested
    },
    respecify_l3_model = function(mobj, new_data) {
      list(
        metadata = new_data[, c("id", "session"), drop = FALSE],
        model_data = data.frame(
          age = c(20, 21, 30, 31),
          group = factor(c("control", "control", "patient", "patient")),
          dummy = seq_len(nrow(new_data))
        ),
        model_variables = c("age", "group")
      )
    },
    emmeans_to_3dlmer_glt = function(mobj, data) list(),
    .package = "fmri.pipeline"
  )

  res <- fmri.pipeline:::afni_3dlmer_setup(
    gpa = gpa,
    backend = NULL,
    lg = lg,
    l1_model_names = "l1_model1",
    l2_model_names = "l2_model1",
    l3_model_names = "l3_model1",
    l2_l3_pairs = NULL,
    subj_df = data.frame(id = c("sub1", "sub2"), session = c(1L, 1L)),
    requires_l2 = TRUE
  )

  expect_true(is.data.frame(res$data))
  expect_equal(nrow(res$data), 1L)
  expect_identical(harvest_call$backend, "afni")
  expect_identical(harvest_call$l1_model_names, "l1_model1")
  expect_identical(harvest_call$l2_model_names, "l2_model1")
  expect_identical(harvest_call$l3_model_names, "l3_model1")
  expect_true(file.exists(res$data$afni_script[1L]))
  expect_identical(res$data$mask_file[1L], mask_file)

  data_table_file <- file.path(dirname(res$data$afni_script[1L]), "dataTable.txt")
  expect_true(file.exists(data_table_file))

  written <- utils::read.delim(data_table_file, stringsAsFactors = FALSE, check.names = FALSE)
  expect_identical(names(written), c("Subj", "session", "age", "group", "InputFile"))
  expect_identical(written$Subj, harvested$l3_model1$contrast1$id)
  expect_identical(written$InputFile, harvested$l3_model1$contrast1$InputFile)
  expect_identical(written$age, c(20L, 21L, 30L, 31L))
  expect_identical(written$group, c("control", "control", "patient", "patient"))

  script_lines <- readLines(res$data$afni_script[1L])
  expect_true(any(grepl("script_dir=", script_lines, fixed = TRUE)))
  expect_true(any(grepl("cd \"$script_dir\"", script_lines, fixed = TRUE)))
  expect_true(any(grepl("-qVars 'age'", script_lines, fixed = TRUE)))
  expect_true(any(grepl("-model 'age+group+(1|Subj)'", script_lines, fixed = TRUE)))
  expect_true(any(grepl("-dataTable @dataTable.txt", script_lines, fixed = TRUE)))
})

test_that("afni_3dlmer_setup builds an intersection mask from FSL L2 inputs when no explicit mask is set", {
  tmp_dir <- tempfile("afni_3dlmer_setup_mask_")
  dir.create(tmp_dir, recursive = TRUE)

  feat_dirs <- file.path(tmp_dir, c("sub1_ses1.gfeat", "sub1_ses2.gfeat", "sub2_ses1.gfeat", "sub2_ses2.gfeat"))
  mask_arrays <- list(
    array(1L, dim = c(2, 2, 2)),
    array(c(0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), dim = c(2, 2, 2)),
    array(c(1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L), dim = c(2, 2, 2)),
    array(1L, dim = c(2, 2, 2))
  )
  for (ii in seq_along(feat_dirs)) {
    dir.create(feat_dirs[ii], recursive = TRUE)
    RNifti::writeNifti(mask_arrays[[ii]], file.path(feat_dirs[ii], "mask.nii.gz"))
  }

  harvested <- list(
    l3_model1 = list(
      contrast1 = data.frame(
        id = c("sub1", "sub1", "sub2", "sub2"),
        session = c(1L, 2L, 1L, 2L),
        l1_model = "l1_model1",
        l2_model = "l2_model1",
        l3_model = "l3_model1",
        l1_cope_name = "cope1",
        l2_cope_name = "copeA",
        feat_dir = feat_dirs,
        InputFile = file.path(tmp_dir, paste0("cope", 1:4, ".nii.gz")),
        source_backend = "fsl",
        stringsAsFactors = FALSE
      )
    )
  )

  l3_obj <- list(
    l3_input_mode = "3dlmer",
    model_variables = c("age", "group"),
    lmer_formula = "age + group + (1 | Subj)",
    lmer_mask = NULL,
    lmer_njobs = NULL
  )
  class(l3_obj) <- c("hi_model_spec", "list")

  gpa <- list(
    output_directory = tmp_dir,
    output_locations = list(
      afni_3dlmer_directory = file.path(
        "{gpa$output_directory}", "afni_3dlmer",
        "L1m-{l1_model}", "l1c-{l1_cope_name}",
        "L2m-{l2_model}_l2c-{l2_cope_name}",
        "L3m-{l3_model}"
      )
    ),
    parallel = list(afni = list(l3_lmer_njobs = 2L)),
    mask_file = NULL,
    l3_models = structure(list(models = list(l3_model1 = l3_obj)), class = c("hi_model_set", "list"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  testthat::local_mocked_bindings(
    harvest_l3_inputs = function(gpa, l3_backend, l1_model_names = NULL, l2_model_names = NULL,
                                 l3_model_names = NULL, lg = NULL) {
      harvested
    },
    respecify_l3_model = function(mobj, new_data) {
      list(
        metadata = new_data[, c("id", "session"), drop = FALSE],
        model_data = data.frame(
          age = c(20, 21, 30, 31),
          group = factor(c("control", "control", "patient", "patient")),
          dummy = seq_len(nrow(new_data))
        ),
        model_variables = c("age", "group")
      )
    },
    emmeans_to_3dlmer_glt = function(mobj, data) list(),
    .package = "fmri.pipeline"
  )

  res <- fmri.pipeline:::afni_3dlmer_setup(
    gpa = gpa,
    backend = NULL,
    lg = lgr::get_logger("test/afni_3dlmer_setup_mask"),
    l1_model_names = "l1_model1",
    l2_model_names = "l2_model1",
    l3_model_names = "l3_model1",
    l2_l3_pairs = NULL,
    subj_df = data.frame(id = c("sub1", "sub2"), session = c(1L, 1L)),
    requires_l2 = TRUE
  )

  expect_true(file.exists(res$data$mask_file[1L]))
  expect_match(res$data$mask_file[1L], "intersection_mask\\.nii\\.gz$")

  mask_img <- RNifti::readNifti(res$data$mask_file[1L])
  expect_identical(dim(mask_img), c(2L, 2L, 2L))
  expect_identical(sum(mask_img != 0), 6L)

  script_lines <- readLines(res$data$afni_script[1L])
  expect_true(any(grepl(paste0("-mask ", res$data$mask_file[1L]), script_lines, fixed = TRUE)))
})
