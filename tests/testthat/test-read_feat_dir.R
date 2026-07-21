make_partial_feat_dir <- function(feat_dir) {
  stats_dir <- file.path(feat_dir, "stats")
  dir.create(stats_dir, recursive = TRUE)
  file.create(file.path(stats_dir, "cope1.nii.gz"))
  file.create(file.path(stats_dir, "zstat2.nii.gz"))
  writeLines(
    c(
      "/ContrastName1\tface-house interaction",
      "/ContrastName2\tcontrol",
      "/NumContrasts\t2"
    ),
    file.path(feat_dir, "design.con")
  )
  invisible(feat_dir)
}

test_that("read_feat_dir retains its return structure for partial outputs", {
  feat_dir <- make_partial_feat_dir(file.path(tempfile("read_feat_dir_"), "model.feat"))
  feat_dir <- normalizePath(feat_dir)
  stats_dir <- file.path(feat_dir, "stats")

  info <- read_feat_dir(feat_dir, statistics = "cope")

  expect_named(
    info,
    c(
      "mask_file", "smoothness_file", "dof_file", "pe_files",
      "cope_files", "varcope_files", "z_files", "zthresh_files",
      "zptfce_files", "rendered_zthresh_files", "t_files",
      "contrast_names", "aux_files", "design_files", "parsed_txt",
      "cope_df"
    )
  )
  expect_named(
    info$cope_df,
    c(
      "cope_number", "contrast_name", "cope", "varcope", "z",
      "zthresh", "rendered_zthresh", "t", "zptfce"
    )
  )
  expect_equal(info$contrast_names, c("face-house_interaction", "control"))
  expect_equal(info$cope_files, c(file.path(stats_dir, "cope1.nii.gz"), NA_character_))
  expect_true(all(is.na(info$varcope_files)))
  expect_true(all(is.na(info$z_files)))
  expect_true(all(is.na(info$t_files)))
})

test_that("read_feat_dir can include expected paths for selected statistics", {
  feat_dir <- make_partial_feat_dir(file.path(tempfile("read_feat_dir_"), "model.feat"))
  feat_dir <- normalizePath(feat_dir)
  stats_dir <- file.path(feat_dir, "stats")

  info <- read_feat_dir(
    feat_dir,
    statistics = c("cope", "zstat"),
    include_missing = TRUE
  )

  expect_equal(
    info$cope_files,
    file.path(stats_dir, paste0("cope", 1:2, ".nii.gz"))
  )
  expect_equal(
    info$z_files,
    file.path(stats_dir, paste0("zstat", 1:2, ".nii.gz"))
  )
  expect_true(all(is.na(info$varcope_files)))
  expect_true(all(is.na(info$t_files)))
  expect_equal(info$cope_df$cope, info$cope_files)
  expect_equal(info$cope_df$z, info$z_files)
})

test_that("read_gfeat_dir passes statistic selection through", {
  gfeat_dir <- file.path(tempfile("read_gfeat_dir_"), "model.gfeat")
  cope_dir <- make_partial_feat_dir(file.path(gfeat_dir, "cope1.feat"))
  cope_dir <- normalizePath(cope_dir)

  info <- read_gfeat_dir(gfeat_dir, statistics = "zstat")

  expect_s3_class(info, "gfeat_info")
  expect_named(info$cope_dirs, "cope1.feat")
  expect_equal(
    info$cope_dirs[[1L]]$z_files,
    c(NA_character_, file.path(cope_dir, "stats", "zstat2.nii.gz"))
  )
  expect_true(all(is.na(info$cope_dirs[[1L]]$cope_files)))
})
