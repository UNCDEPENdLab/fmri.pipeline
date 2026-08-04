make_fwe_execution_test_gpa <- function(root, create_smoothness = TRUE) {
  l3_feat <- file.path(
    root, "feat_l3", "L3m-l3_group", "L1m-facehouse",
    "L2m-l2_session", "l2c-overall", "FEAT_l1c-EV_face.gfeat"
  )
  image_feat_dir <- file.path(l3_feat, "cope1.feat")
  stats_dir <- file.path(image_feat_dir, "stats")
  dir.create(stats_dir, recursive = TRUE)
  file.create(file.path(stats_dir, c("zstat1.nii.gz", "zstat2.nii.gz")))
  file.create(file.path(image_feat_dir, "mask.nii.gz"))
  if (isTRUE(create_smoothness)) {
    writeLines(
      c("DLH 0.1", "VOLUME 1000", "RESELS 100"),
      file.path(stats_dir, "smoothness")
    )
  }

  l3_setup <- data.frame(
    l1_model = "facehouse",
    l1_cope_name = "EV_face",
    l2_model = "l2_session",
    l2_cope_name = "overall",
    l3_model = "l3_group",
    session = 0L,
    l3_input_mode = "all_sessions",
    feat_fsf = sub("\\.gfeat$", ".fsf", l3_feat),
    feat_dir = l3_feat,
    feat_complete = TRUE,
    feat_failed = FALSE,
    feat_dir_exists = TRUE,
    feat_fsf_exists = FALSE,
    stringsAsFactors = FALSE
  )
  l3_metadata <- data.frame(
    id = c("01", "01"),
    session = c(1L, 1L),
    l1_model = c("facehouse", "facehouse"),
    l1_cope_number = c(1L, 1L),
    l1_cope_name = c("EV_face", "EV_face"),
    l2_model = c("l2_session", "l2_session"),
    l2_cope_number = c(1L, 1L),
    l2_cope_name = c("overall", "overall"),
    l3_model = c("l3_group", "l3_group"),
    l3_cope_number = c(1L, 2L),
    l3_cope_name = c("group", "age"),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      output_directory = root,
      scheduler = "slurm",
      lgr_threshold = "warn",
      multi_run = TRUE,
      l3_model_setup = structure(
        list(metadata = l3_metadata, fsl = l3_setup),
        class = c("l3_setup", "list")
      )
    ),
    class = c("glm_pipeline_arguments", "list")
  )
}

make_fake_ptfce_worker <- function(root, fail = FALSE) {
  worker <- file.path(root, if (fail) "fail_ptfce.R" else "fake_ptfce.R")
  if (isTRUE(fail)) {
    lines <- c("#!/usr/bin/env Rscript", "quit(save = 'no', status = 3L)")
  } else {
    lines <- c(
      "#!/usr/bin/env Rscript",
      "args <- commandArgs(trailingOnly = TRUE)",
      "one <- function(flag) args[match(flag, args) + 1L]",
      "zstat <- one('--zstat')",
      "out_dir <- one('--out_dir')",
      "dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)",
      "filename <- basename(zstat)",
      "ext <- if (grepl('\\\\.nii\\\\.gz$', filename)) '.nii.gz' else '.nii'",
      "stem <- substr(filename, 1L, nchar(filename) - nchar(ext))",
      "base <- file.path(out_dir, paste0(stem, '_ptfce'))",
      "file.create(paste0(base, ext))",
      "file.create(paste0(base, '_zthresh.csv'))",
      "if ('--write_thresh_imgs' %in% args) {",
      "  start <- match('--fwep', args) + 1L",
      "  stop_at <- which(seq_along(args) > start & grepl('^--', args))[1L]",
      "  end <- if (is.na(stop_at)) length(args) else stop_at - 1L",
      "  alpha <- args[start:end]",
      "  for (value in alpha) file.create(paste0(base, '_fwep_', value, ext))",
      "}"
    )
  }
  writeLines(lines, worker)
  Sys.chmod(worker, mode = "0755")
  worker
}

test_that("pTFCE execution dry runs expose commands without side effects", {
  root <- tempfile("fwe_execute_dry_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root)
  spec <- fwe_spec(
    "all_maps", targets = list(l3_model = "l3_group"),
    fwe_alpha = c(0.05, 0.01)
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  run <- run_fwe_plan(
    plan, scheduler = "local", dry_run = TRUE, worker_script = worker
  )

  expect_s3_class(run, "fwe_run")
  expect_true(run$dry_run)
  expect_true(all(run$execution$execution_status == "planned"))
  expect_true(all(grepl("--out_dir", run$execution$command, fixed = TRUE)))
  expect_true(all(grepl("--fwep", run$execution$command, fixed = TRUE)))
  expect_false(dir.exists(plan$output_directory))
})

test_that("local pTFCE execution refreshes artifacts and task completion", {
  root <- tempfile("fwe_execute_local_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root)
  spec <- fwe_spec(
    "local_maps", targets = list(l3_model = "l3_group"),
    fwe_alpha = c(0.05, 0.01)
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  # A nested Rscript must not inherit the R CMD check test harness.
  r_tests <- file.path(root, "child_r_tests.R")
  writeLines("stop('R_TESTS leaked into the pTFCE worker')", r_tests)
  withr::local_envvar(c(R_TESTS = r_tests))
  run <- run_fwe_plan(plan, scheduler = "local", worker_script = worker)

  expect_true(all(run$execution$execution_status == "completed"))
  expect_true(all(run$execution$exit_status == 0L))
  expect_true(all(file.exists(run$execution$script_file)))
  expect_true(all(file.exists(run$execution$log_file)))
  expect_true(all(run$plan$tasks$status == "complete"))
  expect_true(all(run$plan$artifacts$exists))
  expect_s3_class(run$result, "fwe_result")
  expect_true(run$result$complete)
  expect_true(file.exists(run$result$manifest_file))
  expect_true(file.exists(run$result$result_file))

  restored <- read_fwe_result(plan$output_directory)
  manifest <- read_fwe_artifact_manifest(plan$output_directory)
  expect_s3_class(restored, "fwe_result")
  expect_true(restored$complete)
  expect_equal(nrow(manifest), 8L)
  expect_identical(anyDuplicated(manifest$artifact_id), 0L)
  expect_true(all(manifest$artifact_status == "available"))

  masks <- fwe_result_artifacts(
    restored,
    artifact_role = "thresholded_statistic",
    fwe_alpha = 0.05
  )
  expect_equal(nrow(masks), 2L)
  expect_true(all(masks$can_use_as_mask))
  expect_true(all(file.exists(masks$artifact_file)))

  skipped <- run_fwe_plan(run$plan, scheduler = "local")
  expect_true(all(skipped$execution$execution_status == "skipped_complete"))
  expect_length(skipped$job_ids, 0L)

  forced <- run_fwe_plan(
    run$plan, scheduler = "local", force = TRUE, dry_run = TRUE,
    worker_script = worker
  )
  expect_true(all(forced$execution$execution_status == "planned"))
})

test_that("blocked FWE tasks can be rejected or explicitly skipped", {
  root <- tempfile("fwe_execute_blocked_")
  gpa <- make_fwe_execution_test_gpa(root, create_smoothness = FALSE)
  spec <- fwe_spec("blocked", targets = list(l3_model = "l3_group"))
  plan <- plan_fwe_correction(
    gpa, spec, source = "setup", require_inputs = FALSE
  )

  expect_error(run_fwe_plan(plan, dry_run = TRUE), "blocked tasks")
  run <- run_fwe_plan(plan, dry_run = TRUE, allow_blocked = TRUE)
  expect_true(all(run$execution$execution_status == "skipped_blocked"))
  expect_true(all(is.na(run$execution$command)))
})

test_that("local execution failures return an inspectable run report", {
  root <- tempfile("fwe_execute_failure_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root, fail = TRUE)
  spec <- fwe_spec("failure", targets = list(l3_cope_name = "age"))
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  run <- run_fwe_plan(
    plan, scheduler = "local", worker_script = worker,
    stop_on_error = FALSE
  )
  expect_identical(run$execution$execution_status, "failed")
  expect_identical(run$execution$exit_status, 3L)

  condition <- tryCatch(
    run_fwe_plan(plan, scheduler = "local", worker_script = worker),
    fwe_execution_error = identity
  )
  expect_s3_class(condition, "fwe_execution_error")
  expect_s3_class(condition$run, "fwe_run")
})

test_that("refresh_fwe_plan detects artifacts created outside the runner", {
  root <- tempfile("fwe_execute_refresh_")
  gpa <- make_fwe_execution_test_gpa(root)
  spec <- fwe_spec(
    "refresh", targets = list(l3_cope_name = "group"),
    write_threshold_images = FALSE
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")
  dir.create(plan$tasks$output_directory, recursive = TRUE)
  file.create(plan$artifacts$artifact_file)

  refreshed <- refresh_fwe_plan(plan)

  expect_identical(refreshed$tasks$status, "complete")
  expect_true(all(refreshed$artifacts$exists))
})

test_that("scheduler execution records submitted job identifiers", {
  root <- tempfile("fwe_execute_scheduler_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root)
  spec <- fwe_spec("scheduler", targets = list(l3_cope_name = "group"))
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  testthat::local_mocked_bindings(
    cluster_submit_shell_jobs = function(job_list, log_file, ...) {
      utils::write.csv(
        data.frame(
          job_number = seq_along(job_list),
          job_id = rep("12345", length(job_list)),
          cmd = job_list,
          stringsAsFactors = FALSE
        ),
        log_file,
        row.names = FALSE
      )
      "12345"
    },
    .package = "fmri.pipeline"
  )

  run <- run_fwe_plan(plan, scheduler = "slurm", worker_script = worker)

  expect_identical(run$execution$execution_status, "submitted")
  expect_identical(run$execution$job_id, "12345")
  expect_identical(run$job_ids, "12345")
  expect_s3_class(run$result, "fwe_result")
  expect_false(run$result$complete)
  expect_true(file.exists(run$result$manifest_file))
  expect_true(all(run$result$manifest$execution_status == "submitted"))
  expect_true(all(run$result$manifest$job_id == "12345"))
})

test_that("local execution safely applies named environment variables", {
  expect_identical(
    fmri.pipeline:::normalize_fwe_run_environment(c(TEST_VALUE = "a b")),
    "TEST_VALUE='a b'"
  )
  expect_error(
    fmri.pipeline:::normalize_fwe_run_environment(c("BAD-NAME" = "x")),
    "invalid environment variable name"
  )
})

test_that("result collection updates a partial persistent manifest", {
  root <- tempfile("fwe_collect_partial_")
  gpa <- make_fwe_execution_test_gpa(root)
  spec <- fwe_spec(
    "partial", targets = list(l3_cope_name = "age"),
    write_threshold_images = FALSE
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  partial <- collect_fwe_results(plan)
  expect_s3_class(partial, "fwe_result")
  expect_false(partial$complete)
  expect_true(all(partial$manifest$artifact_status == "missing_required"))
  expect_error(
    collect_fwe_results(partial, require_complete = TRUE),
    "Required FWE artifacts are missing"
  )

  dir.create(plan$tasks$output_directory, recursive = TRUE)
  file.create(plan$artifacts$artifact_file)
  complete <- collect_fwe_results(partial, require_complete = TRUE)

  expect_true(complete$complete)
  expect_true(all(complete$manifest$exists))
  expect_true(all(complete$manifest$artifact_status == "available"))
  expect_equal(
    read_fwe_artifact_manifest(plan$output_directory),
    complete$manifest,
    ignore_attr = TRUE
  )
})
