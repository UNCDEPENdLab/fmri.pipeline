test_that("scheduler log filenames include scheduler, job id, job name, and extra token", {
  log_file <- fmri.pipeline:::scheduler_log_file(
    "/tmp/batch",
    scheduler = "slurm",
    job_name = "3dlmer_l3 session/EV face",
    job_id = "%j",
    extra = "jobabc"
  )

  expect_equal(log_file, "/tmp/batch/slurm-%j-3dlmer_l3_session_EV_face-jobabc.out")

  directives <- fmri.pipeline:::scheduler_output_directives(
    "slurm",
    "/tmp/batch",
    job_name = "featsep_l1_1",
    extra = "jobabc"
  )
  expect_true(any(grepl("^#SBATCH --output=/tmp/batch/slurm-%j-featsep_l1_1-jobabc\\.out$", directives)))
  expect_true(any(grepl("^#SBATCH --error=/tmp/batch/slurm-%j-featsep_l1_1-jobabc\\.out$", directives)))
})

test_that("job manifest shell helper writes a tsv row in the artifact directory", {
  tmp <- withr::local_tempdir()
  artifact_dir <- file.path(tmp, "model.feat")
  script <- file.path(tmp, "manifest_test.sh")

  lines <- c(
    "#!/bin/bash",
    "job_id=12345",
    "job_name=featsep_l1_1",
    "scheduler_name=slurm",
    paste0("scheduler_log=", shQuote(file.path(tmp, "slurm-12345-featsep_l1_1.out"))),
    paste0("batch_directory=", shQuote(tmp)),
    paste0("batch_file=", shQuote(file.path(tmp, "featsep_l1_1.sbatch"))),
    fmri.pipeline:::job_manifest_shell_function(),
    paste(
      "write_job_manifest",
      shQuote(artifact_dir),
      shQuote("fsl_l1_feat"),
      shQuote("COMPLETED"),
      shQuote(file.path(tmp, "model.fsf"))
    )
  )

  writeLines(lines, script)
  Sys.chmod(script, "0755")
  expect_equal(system2("bash", script), 0L)

  manifest <- file.path(artifact_dir, "job_manifest.tsv")
  expect_true(file.exists(manifest))
  df <- read.delim(manifest, stringsAsFactors = FALSE)
  expect_equal(nrow(df), 1L)
  expect_equal(as.character(df$job_id), "12345")
  expect_equal(df$job_name, "featsep_l1_1")
  expect_equal(df$artifact_type, "fsl_l1_feat")
  expect_equal(df$status, "COMPLETED")
  expect_equal(df$scheduler_log, file.path(tmp, "slurm-12345-featsep_l1_1.out"))
})

test_that("R_batch_job writes named scheduler output directives", {
  batch_dir <- withr::local_tempdir()

  job <- fmri.pipeline::R_batch_job$new(
    batch_directory = batch_dir,
    job_name = "setup run l3 afni",
    r_code = "x <- 1",
    scheduler = "slurm",
    print_session_info = FALSE
  )
  job$generate()

  batch_file <- list.files(batch_dir, pattern = "^submit_batch.*\\.sh$", full.names = TRUE)
  expect_length(batch_file, 1L)
  lines <- readLines(batch_file, warn = FALSE)

  expect_true(any(grepl("^#SBATCH -J setup_run_l3_afni$", lines)))
  expect_true(any(grepl("slurm-%j-setup_run_l3_afni-submit_batch_setup_run_l3_afni\\.out", lines)))
})
