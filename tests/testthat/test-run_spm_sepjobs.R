test_that("run_spm_sepjobs final tracking reuses one terminal job status", {
  spm_path <- normalizePath(file.path(pkg_root, "R", "spm_utils.R"), mustWork = TRUE)
  lines <- readLines(spm_path, warn = FALSE)

  expect_true(any(grepl("if [ $job_failed -ne 0 ]; then", lines, fixed = TRUE)))
  expect_true(any(grepl("job_status=\\\"FAILED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("job_status=\\\"COMPLETED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("--status", lines) & grepl("\\$\\{job_status\\}", lines)))
  expect_true(any(grepl("if [ \\\"${job_status}\\\" = \\\"FAILED\\\" ]; then", lines, fixed = TRUE)))
})
