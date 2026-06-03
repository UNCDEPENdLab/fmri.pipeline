test_that("run_feat_sepjobs source wires L3 execution to outlier retry wrapper", {
  run_feat_path <- normalizePath(file.path(pkg_root, "R", "run_feat_sepjobs.R"), mustWork = TRUE)
  lines <- readLines(run_feat_path, warn = FALSE)

  expect_true(any(grepl("feat_retry_binary <- system.file\\(\"bash/run_feat_with_outlier_retry\\.sh\"", lines)))
  expect_true(any(grepl("auto_retry_l3_excessive_outliers", lines)))
  expect_true(any(grepl("level == 3L && isTRUE\\(auto_retry_l3_excessive_outliers\\)", lines)))
  expect_true(any(grepl("feat_call_with_jobs", lines) & grepl("feat_retry_binary", lines)))
  expect_true(any(grepl("return \\$exit_code", lines)))
})

test_that("run_feat_sepjobs waits on each background FEAT child", {
  run_feat_path <- normalizePath(file.path(pkg_root, "R", "run_feat_sepjobs.R"), mustWork = TRUE)
  lines <- readLines(run_feat_path, warn = FALSE)

  expect_true(any(grepl("job_failed=0", lines, fixed = TRUE)))
  expect_true(any(grepl("pids=()", lines, fixed = TRUE)))
  expect_true(any(grepl("pids+=(\\\"$!\\\")", lines, fixed = TRUE)))
  expect_true(any(grepl("for pid in \\\"${pids[@]}\\\"; do", lines, fixed = TRUE)))
  expect_true(any(grepl("if ! wait \\\"$pid\\\"; then", lines, fixed = TRUE)))
  expect_true(any(grepl("job_failed=1", lines, fixed = TRUE)))
  expect_true(any(grepl("job_status=\\\"FAILED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("job_status=\\\"COMPLETED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("--status", lines) & grepl("\\$\\{job_status\\}", lines)))
})
