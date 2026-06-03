test_that("run_3dlmer_sepjobs wrapper reuses one terminal job status", {
  run_3dlmer_path <- normalizePath(file.path(pkg_root, "R", "run_3dlmer_sepjobs.R"), mustWork = TRUE)
  lines <- readLines(run_3dlmer_path, warn = FALSE)

  expect_true(any(grepl("job_status=\\\"COMPLETED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("job_status=\\\"FAILED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("write_job_manifest", lines) & grepl("\\$\\{job_status\\}", lines)))
  expect_true(any(grepl("Rscript", lines) & grepl("--status", lines) & grepl("\\$\\{job_status\\}", lines)))
  expect_equal(sum(grepl("write_job_manifest", lines) & grepl("afni_3dlmer", lines)), 1L)
})
