test_that("flame_runner uses parallel and clamps to allocation", {
  tmp <- withr::local_tempdir()
  stub_dir <- file.path(tmp, "stubs_parallel")
  dir.create(stub_dir, recursive = TRUE, showWarnings = FALSE)

  make_stub <- function(name, lines) {
    path <- file.path(stub_dir, name)
    writeLines(lines, path)
    Sys.chmod(path, "0755")
    path
  }

  make_stub("sh", c(
    "#!/bin/bash",
    "exec /bin/sh \"$@\""
  ))
  make_stub("grep", c(
    "#!/bin/bash",
    "exec /bin/grep \"$@\""
  ))

  parallel_log <- file.path(tmp, "parallel.log")
  make_stub("parallel", c(
    "#!/bin/bash",
    "jobs=1",
    "while [[ $# -gt 0 ]]; do",
    "  case \"$1\" in",
    "    --jobs|-j) jobs=\"$2\"; shift 2;;",
    "    *) shift;;",
    "  esac",
    "done",
    "echo \"$jobs\" > \"$PARALLEL_LOG\"",
    "while IFS= read -r line; do",
    "  if [ -n \"${line//[[:space:]]/}\" ]; then",
    "    sh -c \"$line\"",
    "  fi",
    "done"
  ))

  flame_list <- file.path(tmp, "flame.list")
  out_file <- file.path(tmp, "out.txt")
  writeLines(c(
    paste("echo A >>", out_file),
    paste("echo B >>", out_file),
    paste("echo C >>", out_file)
  ), flame_list)

  withr::local_envvar(c(
    PATH = stub_dir,
    PARALLEL_LOG = parallel_log,
    SLURM_CPUS_PER_TASK = NA_character_,
    SLURM_NTASKS = "2",
    SLURM_CPUS_ON_NODE = NA_character_,
    PBS_NP = NA_character_,
    NCPUS = NA_character_,
    PBS_NUM_PPN = NA_character_
  ))

  flame_runner <- normalizePath(file.path(pkg_root, "inst/bin/flame_runner"), mustWork = TRUE)
  system2(flame_runner, c(flame_list, "999"), stdout = TRUE, stderr = TRUE)

  expect_equal(readLines(parallel_log), "2")
  expect_true(file.exists(out_file))
  expect_equal(length(readLines(out_file)), 3)
})

test_that("flame_runner uses xargs when parallel is unavailable", {
  tmp <- withr::local_tempdir()
  stub_dir <- file.path(tmp, "stubs_xargs")
  dir.create(stub_dir, recursive = TRUE, showWarnings = FALSE)

  make_stub <- function(name, lines) {
    path <- file.path(stub_dir, name)
    writeLines(lines, path)
    Sys.chmod(path, "0755")
    path
  }

  make_stub("sh", c(
    "#!/bin/bash",
    "exec /bin/sh \"$@\""
  ))
  make_stub("grep", c(
    "#!/bin/bash",
    "exec /bin/grep \"$@\""
  ))

  xargs_log <- file.path(tmp, "xargs.log")
  make_stub("xargs", c(
    "#!/bin/bash",
    "if [[ \"$1\" == \"--help\" ]]; then",
    "  echo \"  -d delim\"",
    "  exit 0",
    "fi",
    "p=1",
    "while [[ $# -gt 0 ]]; do",
    "  case \"$1\" in",
    "    -P) p=\"$2\"; shift 2;;",
    "    -d) shift 2;;",
    "    -n) shift 2;;",
    "    -I) shift 2;;",
    "    *) shift;;",
    "  esac",
    "done",
    "echo \"$p\" > \"$XARGS_LOG\"",
    "while IFS= read -r line; do",
    "  if [ -n \"${line//[[:space:]]/}\" ]; then",
    "    sh -c \"$line\"",
    "  fi",
    "done"
  ))

  flame_list <- file.path(tmp, "flame.list")
  out_file <- file.path(tmp, "out.txt")
  writeLines(c(
    paste("echo A >>", out_file),
    paste("echo B >>", out_file)
  ), flame_list)

  withr::local_envvar(c(
    PATH = stub_dir,
    XARGS_LOG = xargs_log
  ))

  flame_runner <- normalizePath(file.path(pkg_root, "inst/bin/flame_runner"), mustWork = TRUE)
  system2(flame_runner, c(flame_list, "4"), stdout = TRUE, stderr = TRUE)

  expect_equal(readLines(xargs_log), "2")
  expect_true(file.exists(out_file))
  expect_equal(length(readLines(out_file)), 2)
})

test_that("flame_runner falls back when xargs lacks -d support", {
  tmp <- withr::local_tempdir()
  stub_dir <- file.path(tmp, "stubs_fallback")
  dir.create(stub_dir, recursive = TRUE, showWarnings = FALSE)

  make_stub <- function(name, lines) {
    path <- file.path(stub_dir, name)
    writeLines(lines, path)
    Sys.chmod(path, "0755")
    path
  }

  make_stub("sh", c(
    "#!/bin/bash",
    "exec /bin/sh \"$@\""
  ))
  make_stub("grep", c(
    "#!/bin/bash",
    "exec /bin/grep \"$@\""
  ))
  make_stub("xargs", c(
    "#!/bin/bash",
    "if [[ \"$1\" == \"--help\" ]]; then",
    "  echo \"xargs help\"",
    "  exit 0",
    "fi",
    "exit 1"
  ))

  flame_list <- file.path(tmp, "flame.list")
  out_file <- file.path(tmp, "out.txt")
  writeLines(c(
    paste("echo A >>", out_file),
    paste("echo B >>", out_file)
  ), flame_list)

  withr::local_envvar(c(
    PATH = stub_dir
  ))

  flame_runner <- normalizePath(file.path(pkg_root, "inst/bin/flame_runner"), mustWork = TRUE)
  system2(flame_runner, c(flame_list, "10"), stdout = TRUE, stderr = TRUE)

  expect_true(file.exists(out_file))
  expect_equal(length(readLines(out_file)), 2)
})
