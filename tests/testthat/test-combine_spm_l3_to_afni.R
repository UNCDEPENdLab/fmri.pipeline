test_that("combine_spm_l3_to_afni forces +tlrc output view", {
  tmp_dir <- tempfile("spm_l3_afni_")
  dir.create(tmp_dir, recursive = TRUE)

  spm_dir <- file.path(tmp_dir, "spm_l3", "L1m-m1", "L3m-l3")
  dir.create(spm_dir, recursive = TRUE)

  con_file <- file.path(spm_dir, "con_0001.nii")
  stat_file <- file.path(spm_dir, "spmT_0001.nii")
  file.create(con_file)
  file.create(stat_file)

  gpa <- create_mock_gpa(output_directory = tmp_dir)
  gpa$l3_model_setup <- list(
    spm = data.frame(
      l1_model = "m1",
      l1_cope_name = "c1",
      l3_model = "l3",
      spm_dir = spm_dir,
      stringsAsFactors = FALSE
    )
  )
  gpa$output_locations$spm_l3_combined_filename <- "{gpa$output_directory}/spm_l3_combined/L1m-{l1_model}/L3m-{l3_model}/l1c-{l1_cope_name}_stats"
  gpa$output_locations$spm_l3_combined_briknames <- "l3c-{l3_cope_name}"

  afni_cmds <- character(0)
  expected_out <- file.path(tmp_dir, "spm_l3_combined", "L1m-m1", "L3m-l3", "l1c-c1_stats+tlrc")

  testthat::local_mocked_bindings(
    spm_l3_collect_contrasts = function(spm_dir, lg = NULL) {
      data.frame(
        l3_cope_number = 1L,
        l3_cope_name = "contrast_1",
        stat = "T",
        df1 = 10,
        df2 = NA_real_,
        con_file = con_file,
        stat_file = stat_file,
        stringsAsFactors = FALSE
      )
    },
    run_afni_command = function(args, ...) {
      afni_cmds <<- c(afni_cmds, args)
      0L
    },
    .package = "fmri.pipeline"
  )

  res <- fmri.pipeline:::combine_spm_l3_to_afni(gpa)
  expect_equal(nrow(res), 1L)
  expect_true(any(grepl(paste0("3dTcat -overwrite -prefix ", expected_out), afni_cmds, fixed = TRUE)))
  expect_true(any(grepl("3drefit -fbuc", afni_cmds, fixed = TRUE)))
  expect_false(any(grepl("\\+orig\\b", afni_cmds)))
})

test_that("combine_spm_l3_to_afni errors when AFNI command fails", {
  tmp_dir <- tempfile("spm_l3_afni_fail_")
  dir.create(tmp_dir, recursive = TRUE)

  spm_dir <- file.path(tmp_dir, "spm_l3", "L1m-m1", "L3m-l3")
  dir.create(spm_dir, recursive = TRUE)

  con_file <- file.path(spm_dir, "con_0001.nii")
  stat_file <- file.path(spm_dir, "spmT_0001.nii")
  file.create(con_file)
  file.create(stat_file)

  gpa <- create_mock_gpa(output_directory = tmp_dir)
  gpa$l3_model_setup <- list(
    spm = data.frame(
      l1_model = "m1",
      l1_cope_name = "c1",
      l3_model = "l3",
      spm_dir = spm_dir,
      stringsAsFactors = FALSE
    )
  )
  gpa$output_locations$spm_l3_combined_filename <- "{gpa$output_directory}/spm_l3_combined/L1m-{l1_model}/L3m-{l3_model}/l1c-{l1_cope_name}_stats"
  gpa$output_locations$spm_l3_combined_briknames <- "l3c-{l3_cope_name}"

  testthat::local_mocked_bindings(
    spm_l3_collect_contrasts = function(spm_dir, lg = NULL) {
      data.frame(
        l3_cope_number = 1L,
        l3_cope_name = "contrast_1",
        stat = "T",
        df1 = 10,
        df2 = NA_real_,
        con_file = con_file,
        stat_file = stat_file,
        stringsAsFactors = FALSE
      )
    },
    run_afni_command = function(args, ...) {
      if (startsWith(args, "3dTcat")) return(1L)
      0L
    },
    .package = "fmri.pipeline"
  )

  expect_error(fmri.pipeline:::combine_spm_l3_to_afni(gpa), "AFNI command failed")
})
