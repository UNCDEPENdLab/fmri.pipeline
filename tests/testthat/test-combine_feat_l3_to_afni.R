make_fake_l3_cope_feat <- function(feat_dir, contrast_names) {
  cope_feat <- file.path(feat_dir, "cope1.feat")
  stats_dir <- file.path(cope_feat, "stats")
  dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

  for (ii in seq_along(contrast_names)) {
    file.create(file.path(stats_dir, sprintf("cope%d.nii.gz", ii)))
    file.create(file.path(stats_dir, sprintf("varcope%d.nii.gz", ii)))
    file.create(file.path(stats_dir, sprintf("zstat%d.nii.gz", ii)))
    file.create(file.path(stats_dir, sprintf("tstat%d.nii.gz", ii)))
  }

  writeLines(
    paste(sprintf("/ContrastName%d %s", seq_along(contrast_names), contrast_names), collapse = "\n"),
    file.path(cope_feat, "design.con")
  )
}

test_that("combine_feat_l3_to_afni groups per-cope L2 contrasts into one useful AFNI file", {
  tmp_dir <- tempfile("feat_l3_afni_")
  dir.create(tmp_dir, recursive = TRUE)

  l2_copes <- c("session_mean", "session_linear")
  feat_dirs <- file.path(tmp_dir, "feat_l3", paste0("l2c-", l2_copes), "FEAT_l1c-face.gfeat")
  for (feat_dir in feat_dirs) {
    make_fake_l3_cope_feat(feat_dir, c("group_mean", "age_slope"))
  }

  gpa <- create_mock_gpa(output_directory = tmp_dir)
  gpa$multi_run <- TRUE
  gpa$l3_model_setup <- list(
    fsl = data.frame(
      l1_model = "facehouse",
      l1_cope_name = "face",
      l2_model = "session_model",
      l2_cope_name = l2_copes,
      l3_model = "group_model",
      feat_dir = feat_dirs,
      stringsAsFactors = FALSE
    )
  )
  gpa$output_locations$feat_l3_combined_filename <- file.path(
    "{gpa$output_directory}", "feat_l3_combined", "L3m-{l3_model}",
    "L1m-{l1_model}", "L2m-{l2_model}", "l1c-{l1_cope_name}", "stats"
  )
  gpa$output_locations$feat_l3_combined_briknames <- "l2c-{l2_cope_name}_l3c-{l3_cope_name}"

  afni_cmds <- character(0)
  expected_out <- file.path(
    tmp_dir, "feat_l3_combined", "L3m-group_model", "L1m-facehouse",
    "L2m-session_model", "l1c-face", "stats+tlrc"
  )
  testthat::local_mocked_bindings(
    run_afni_command = function(args, ...) {
      afni_cmds <<- c(afni_cmds, args)
      if (startsWith(args, "3dTcat")) file.create(paste0(expected_out, ".HEAD"))
      0L
    },
    .package = "fmri.pipeline"
  )

  res <- fmri.pipeline::combine_feat_l3_to_afni(gpa)
  tcat_cmds <- grep("^3dTcat", afni_cmds, value = TRUE)
  refit_cmds <- grep("^3drefit", afni_cmds, value = TRUE)

  expect_equal(nrow(res), 4L)
  expect_length(unique(res$afni_out), 1L)
  expect_false(any(grepl("l2c-", unique(res$afni_out), fixed = TRUE)))
  expect_length(tcat_cmds, 1L)
  expect_equal(length(gregexpr("\\.nii\\.gz", tcat_cmds, fixed = FALSE)[[1]]), 8L)
  expect_true(any(grepl("l2c-session_mean_l3c-group_mean_cope", refit_cmds, fixed = TRUE)))
  expect_true(any(grepl("l2c-session_linear_l3c-age_slope_z", refit_cmds, fixed = TRUE)))
})

test_that("combine_feat_l3_to_afni relabels the actual AFNI view created by 3dTcat", {
  tmp_dir <- tempfile("feat_l3_afni_view_")
  dir.create(tmp_dir, recursive = TRUE)

  feat_dir <- file.path(tmp_dir, "feat_l3", "l2c-overall", "FEAT_l1c-face.gfeat")
  make_fake_l3_cope_feat(feat_dir, "group_mean")

  gpa <- create_mock_gpa(output_directory = tmp_dir)
  gpa$multi_run <- TRUE
  gpa$l3_model_setup <- list(
    fsl = data.frame(
      l1_model = "facehouse",
      l1_cope_name = "face",
      l2_model = "session_model",
      l2_cope_name = "overall",
      l3_model = "group_model",
      feat_dir = feat_dir,
      stringsAsFactors = FALSE
    )
  )
  gpa$output_locations$feat_l3_combined_filename <- file.path(
    "{gpa$output_directory}", "feat_l3_combined", "L3m-{l3_model}",
    "L1m-{l1_model}", "L2m-{l2_model}", "l1c-{l1_cope_name}", "stats"
  )
  gpa$output_locations$feat_l3_combined_briknames <- "l2c-{l2_cope_name}_l3c-{l3_cope_name}"

  requested_out <- file.path(
    tmp_dir, "feat_l3_combined", "L3m-group_model", "L1m-facehouse",
    "L2m-session_model", "l1c-face", "stats+tlrc"
  )
  actual_out <- sub("\\+tlrc$", "+orig", requested_out)

  afni_cmds <- character(0)
  testthat::local_mocked_bindings(
    run_afni_command = function(args, ...) {
      afni_cmds <<- c(afni_cmds, args)
      if (startsWith(args, "3dTcat")) file.create(paste0(actual_out, ".HEAD"))
      0L
    },
    .package = "fmri.pipeline"
  )

  fmri.pipeline::combine_feat_l3_to_afni(gpa)
  refit_cmds <- grep("^3drefit", afni_cmds, value = TRUE)

  expect_length(refit_cmds, 1L)
  expect_true(grepl(actual_out, refit_cmds, fixed = TRUE))
  expect_false(grepl(requested_out, refit_cmds, fixed = TRUE))
})
