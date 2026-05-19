test_that("FLAME fallback TSV records are promoted to lgr estimation logs", {
  tmp <- withr::local_tempdir()
  feat_dir <- file.path(tmp, "cope1.feat")
  dir.create(feat_dir, recursive = TRUE)

  fallback_tsv <- file.path(feat_dir, "flame_runner_fallbacks.tsv")
  writeLines(c(
    "timestamp_utc\tpid\tflame12_status\tld_dir\tfailed_ld_dir\tfallback_command",
    paste(
      "2026-05-19T14:04:02Z",
      "75",
      "134",
      "stats0007",
      "stats0007.flame12_failed.75",
      "flameo --cope=tmpcope0007 --ld=stats0007 --runmode=flame1",
      sep = "\t"
    )
  ), fallback_tsv)

  log_txt <- file.path(tmp, "logs", "l3_estimation.txt")
  script <- normalizePath(
    file.path(pkg_root, "inst/bin/log_flame_runner_fallbacks.R"),
    mustWork = TRUE
  )

  status <- system2(
    "Rscript",
    c(
      script,
      "--feat_dir", feat_dir,
      "--level", "3",
      "--log_txt", log_txt,
      "--log_json", "NULL",
      "--threshold", "warn"
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  expect_false(inherits(status, "status"))
  expect_true(file.exists(log_txt))
  log_lines <- readLines(log_txt, warn = FALSE)
  expect_true(any(grepl("FSL FLAME12 failed", log_lines, fixed = TRUE)))
  expect_true(any(grepl("rerunning it with FLAME1", log_lines, fixed = TRUE)))
  expect_true(any(grepl("stats0007", log_lines, fixed = TRUE)))
})

test_that("run_feat_sepjobs source wires FLAME fallback logging into FEAT jobs", {
  run_feat_path <- normalizePath(file.path(pkg_root, "R", "run_feat_sepjobs.R"), mustWork = TRUE)
  lines <- readLines(run_feat_path, warn = FALSE)

  expect_true(any(grepl("log_flame_runner_fallbacks\\.R", lines)))
  expect_true(any(grepl("--feat_dir", lines) & grepl("\\$\\{odir\\}", lines)))
  expect_true(any(grepl("l%d_estimation.txt", lines, fixed = TRUE)))
  expect_true(any(grepl("l%d_estimation.json", lines, fixed = TRUE)))
})
