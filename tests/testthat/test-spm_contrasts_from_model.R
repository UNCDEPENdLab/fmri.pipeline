test_that("generate_spm_contrasts_from_model maps L1 contrasts to SPM columns", {
  skip_if_not_installed("R.matlab")

  tmp_dir <- tempfile("spm_contrast_model_")
  dir.create(tmp_dir, recursive = TRUE)

  # Fake SPM design columns (two runs, two conditions)
  mnames <- c(
    "Sn(1) feedback.A*bf(1)",
    "Sn(1) feedback.B*bf(1)",
    "Sn(2) feedback.A*bf(1)",
    "Sn(2) feedback.B*bf(1)",
    "Sn(1) constant",
    "Sn(2) constant"
  )
  cpos <- 1:4
  bpos <- 5:6
  npos <- integer(0)

  matfile <- file.path(tmp_dir, "design_columns.mat")
  R.matlab::writeMat(matfile, mnames = mnames, cpos = cpos, bpos = bpos, npos = npos)

  # Minimal l1_model_spec with contrasts
  cmat <- matrix(c(1, -1), nrow = 1, dimnames = list("A_gt_B", c("feedback.A", "feedback.B")))
  mobj <- list(contrasts = cmat)
  class(mobj) <- c("list", "l1_model_spec")

  # Generate spec + commands
  cmds <- fmri.pipeline:::generate_spm_contrasts_from_model(
    output_dir = tmp_dir,
    mobj = mobj,
    spm_path = tmp_dir,
    execute = FALSE,
    average_across_runs = TRUE
  )

  spec_path <- file.path(tmp_dir, "spm_contrast_spec.rds")
  expect_true(file.exists(spec_path))
  expect_true(file.exists(file.path(tmp_dir, "extract_design_columns.m")))

  # Run contrast script directly with the prewritten .mat file
  setup_script <- system.file("Rscript", "setup_spm_contrasts_from_model.R", package = "fmri.pipeline")
  expect_true(file.exists(setup_script))
  cmd <- paste(
    "Rscript --no-save --no-restore",
    shQuote(setup_script),
    "-mat_file", shQuote(matfile),
    "-contrast_rds", shQuote(spec_path),
    "-average_across_runs TRUE",
    "-spm_path", shQuote(tmp_dir)
  )
  system(cmd)

  mfile <- file.path(tmp_dir, "estimate_glm_contrasts.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  convec_line <- lines[grepl("tcon.convec", lines)][1]
  expect_true(nchar(convec_line) > 0)

  # Extract numeric weights from convec line
  weights <- as.numeric(strsplit(gsub(".*\\[|\\].*", "", convec_line), ",")[[1]])
  expect_equal(length(weights), length(mnames))
  # Average across runs => 0.5 for feedback.A in each run, -0.5 for feedback.B in each run
  expect_equal(weights[1:4], c(0.5, -0.5, 0.5, -0.5), tolerance = 1e-8)
  expect_equal(weights[5:6], c(0, 0))
})

test_that("generate_spm_contrasts_from_model errors on missing contrasts", {
  mobj <- list(contrasts = NULL)
  class(mobj) <- c("list", "l1_model_spec")
  expect_error(
    fmri.pipeline:::generate_spm_contrasts_from_model(output_dir = tempdir(), mobj = mobj),
    "missing or empty"
  )
})

test_that("generate_spm_contrasts_from_model matches compact SPM pmod labels", {
  skip_if_not_installed("R.matlab")

  tmp_dir <- tempfile("spm_contrast_pmod_")
  dir.create(tmp_dir, recursive = TRUE)

  # Compact pmod form from SPM (no spaces around x)
  mnames <- c(
    "Sn(1) choice*bf(1)",
    "Sn(1) choicexchoice_trial-pmod^1*bf(1)",
    "Sn(1) feedback*bf(1)",
    "Sn(1) constant"
  )
  cpos <- 1:3
  bpos <- 4
  npos <- integer(0)

  matfile <- file.path(tmp_dir, "design_columns.mat")
  R.matlab::writeMat(matfile, mnames = mnames, cpos = cpos, bpos = bpos, npos = npos)

  cmat <- matrix(
    c(0, 0, 1),
    nrow = 1,
    dimnames = list("EV_choice_trial", c("choice", "feedback", "choice_trial"))
  )
  mobj <- list(contrasts = cmat)
  class(mobj) <- c("list", "l1_model_spec")

  fmri.pipeline:::generate_spm_contrasts_from_model(
    output_dir = tmp_dir,
    mobj = mobj,
    spm_path = tmp_dir,
    execute = FALSE,
    average_across_runs = TRUE
  )

  spec_path <- file.path(tmp_dir, "spm_contrast_spec.rds")
  setup_script <- system.file("Rscript", "setup_spm_contrasts_from_model.R", package = "fmri.pipeline")
  cmd <- paste(
    "Rscript --no-save --no-restore",
    shQuote(setup_script),
    "-mat_file", shQuote(matfile),
    "-contrast_rds", shQuote(spec_path),
    "-average_across_runs TRUE",
    "-spm_path", shQuote(tmp_dir)
  )
  system(cmd)

  mfile <- file.path(tmp_dir, "estimate_glm_contrasts.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  expect_true(any(grepl("EV_choice_trial", lines, fixed = TRUE)))

  convec_line <- lines[grepl("tcon.convec", lines)][1]
  weights <- as.numeric(strsplit(gsub(".*\\[|\\].*", "", convec_line), ",")[[1]])
  expect_equal(weights, c(0, 1, 0, 0), tolerance = 1e-8)
})
