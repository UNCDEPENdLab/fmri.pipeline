test_that("run_3dlmer_sepjobs wrapper reuses one terminal job status", {
  run_3dlmer_path <- source_tree_file("R", "run_3dlmer_sepjobs.R")
  lines <- readLines(run_3dlmer_path, warn = FALSE)

  expect_true(any(grepl("job_status=\\\"COMPLETED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("job_status=\\\"FAILED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("write_job_manifest", lines) & grepl("\\$\\{job_status\\}", lines)))
  expect_true(any(grepl("Rscript", lines) & grepl("--status", lines) & grepl("\\$\\{job_status\\}", lines)))
  expect_equal(sum(grepl("write_job_manifest", lines) & grepl("afni_3dlmer", lines)), 1L)
})

test_that("run_3dlmer_sepjobs writes slurm logs to the active pipeline batch directory", {
  tmp <- withr::local_tempdir()
  launch_dir <- file.path(tmp, "launch")
  batch_dir <- file.path(tmp, "scheduler_scripts", "batch_active")
  afni_dir <- file.path(tmp, "afni_3dlmer", "L3m-l3_model", "l1c-copeA_l2c-copeB")
  dir.create(launch_dir, recursive = TRUE)
  dir.create(batch_dir, recursive = TRUE)
  dir.create(afni_dir, recursive = TRUE)

  afni_script <- file.path(afni_dir, "run_3dlmer.sh")
  writeLines(c("#!/bin/bash", "exit 0"), afni_script)

  gpa <- create_mock_gpa(output_directory = tmp)
  gpa$scheduler <- "slurm"
  gpa$batch_run <- list(batch_directory = batch_dir)
  gpa$output_locations$sqlite_db <- file.path(tmp, "jobs.sqlite")
  gpa$parallel$sched_args <- NULL
  gpa$parallel$afni <- list(
    l3_lmer_time = "1:00:00",
    l3_lmer_memgb = "2",
    l3_lmer_njobs = 1L
  )
  gpa$l3_model_setup <- list(
    afni = data.frame(
      afni_script = afni_script,
      l3_model = "l3_model",
      l1_cope_name = "copeA",
      l2_cope_name = "copeB",
      afni_complete = FALSE,
      stringsAsFactors = FALSE
    )
  )

  captured_script <- NULL
  withr::local_envvar(c(SLURM_SUBMIT_DIR = launch_dir))
  job_ids <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_3dlmer_sepjobs(gpa),
    cluster_job_submit = function(script, ...) {
      captured_script <<- script
      "12345"
    },
    get_compute_environment = function(...) character(0),
    .package = "fmri.pipeline"
  )

  expect_equal(job_ids, "12345")
  expect_true(file.exists(captured_script))

  lines <- readLines(captured_script, warn = FALSE)
  output_line <- lines[grepl("^#SBATCH --output=", lines)]
  error_line <- lines[grepl("^#SBATCH --error=", lines)]
  runtime_line <- lines[grepl("^scheduler_log=", lines)]

  expect_length(output_line, 1L)
  expect_true(startsWith(output_line, paste0("#SBATCH --output=", batch_dir)))
  expect_true(startsWith(error_line, paste0("#SBATCH --error=", batch_dir)))
  expect_true(grepl(batch_dir, runtime_line, fixed = TRUE))
  expect_false(grepl(launch_dir, output_line, fixed = TRUE))
})
