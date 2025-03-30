#!/usr/bin/env Rscript
# 
# Command-line version of 'update_tracked_job_status'
# Used for jobs submitted to 'cluster_job_submit' that are command-line only (e.g., not nested in R_batch_job)
# retreive in code using:   upd_job_status_path <- system.file("bin/upd_job_status.R", package = "fmri.pipeline")
# example status update: paste("Rscript", upd_job_status_path, "--job_id" , JOBID, "--sqlite_db", SQLITE_DB, "--status", STATUS)

print_help <- function() {
  cat(paste("This script makes a call to the `update_tracked_job_status` command and is to be",
            "used internally in command-line only scripts submitted to `cluster_job_submit`.",
            "Options:",
            "  --job_id <job_id>: The job id of the job whose status is to be updated.",
            "  --sqlite_db <sqlite_db>: The path to the tracking SQLite database",
            "  --status <status>: The new status of the job specified by `--job_id`.",
            "                     One of QUEUED, STARTED, COMPLETED, FAILED, or FAILED_BY_EXT",
            "  --cascade: If specified, then new status will cascade to child jobs.",
            "             Only works for status 'FAILED'/'FAILED_BY_EXT'",
            "  --help: Print the help menu",
            "\n\n",
            sep = "\n"
  ))
}

#read in command line arguments
args <- commandArgs(trailingOnly = FALSE)

# check if call to "--args" is explicit
argpos <- grep("--args", args, fixed=TRUE)
if (length(argpos) > 0L) {
  args <- args[(argpos + 1):length(args)]
} else {
  args <- c()
}

# load required pacakges
if (!require("checkmate")) {
  install.packages("checkmate")
  library(checkmate)
}

# set default args
job_id <- NULL
sqlite_db <- NULL
status <- NULL
cascade <- FALSE

# gather arguments
argpos <- 1
while (argpos <= length(args)) {
  #print(args[argpos])
  if (args[argpos] == "--job_id") {
    job_id <- args[argpos + 1] # job_id of job to update
    argpos <- argpos + 2
  } else if (args[argpos] == "--sqlite_db") {
    sqlite_db <- args[argpos + 1] # path to the sqlite db
    argpos <- argpos + 2
  } else if (args[argpos] == "--status") {
    status <- args[argpos + 1] # new status
    argpos <- argpos + 2
  } else if (args[argpos] == "--help") {
    print_help()
    quit(save = "no", 0, FALSE)
  } else if (args[argpos] == "--cascade") {
    cascade <- TRUE
    argpos <- argpos + 1
  } else {
    stop("Not sure what to do with argument: ", args[argpos])
  }
}

# convert string versions of NULL to regular NULL
if (checkmate::test_set_equal(job_id, "NULL")) job_id <- NULL
if (checkmate::test_set_equal(sqlite_db, "NULL")) sqlite_db <- NULL
if (checkmate::test_set_equal(status, "NULL")) status <- NULL

fmri.pipeline::update_tracked_job_status(job_id = job_id, 
                                         sqlite_db = sqlite_db, 
                                         status = status, 
                                         cascade = cascade)