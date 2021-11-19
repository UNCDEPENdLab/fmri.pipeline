#' This function pauses execution of an R script while a scheduled qsub job is not yet complete.
#'
#' It is intended to give you control over job dependencies within R when the formal PBS
#' depend approach is insufficient, especially in the case of a script that spawns child jobs that
#' need to be scheduled or complete before the parent script should continue.
#'
#' @param job_ids One or more job ids of existing PBS or slurm jobs, or process ids of a local process for
#'   \code{scheduler="sh"}.
#' @param repolling_interval How often to recheck the job status, in seconds. Default: 30
#' @param max_wait How long to wait on the job before giving up, in seconds. Default: 24 hours (86,400 seconds)
#' @param scheduler What scheduler is used for job execution.
#'   Options: c("torque", "qsub", "slurm", "sbatch", "sh", "local")
#' @param quiet If \code{TRUE}, \code{wait_for_job} will not print out any status updates on jobs. If \code{FALSE},
#'   the function prints out status updates for each tracked job so that the user knows what's holding up progress.
#'
#' @return Nothing. Just returns when the blocking job completes.
#'
#' @details Note that for the \code{scheduler} argument, "torque" and "qsub" are the same;
#'   "slurm" and "sbatch" are the same, and "sh" and "local" are the same.
#' @examples
#' \dontrun{
#' # example on qsub/torque cluster
#' wait_for_job("7968857.torque01.util.production.int.aci.ics.psu.edu", scheduler = "torque")
#'
#' # example of waiting for two jobs on slurm cluster
#' wait_for_job(c("24147864", "24147876"), scheduler = "slurm")
#'
#' # example of waiting for two jobs on local machine
#' wait_for_job(c("9843", "9844"), scheduler = "local")
#' }
#'
#' @author Michael Hallquist
#' @importFrom dplyr full_join
#' @export
wait_for_job <- function(job_ids, repolling_interval = 60, max_wait = 60 * 60 * 24,
                         scheduler = "local", quiet = TRUE, stop_on_timeout = TRUE) {
  checkmate::assert_number(repolling_interval, lower = 0.1, upper = 2e5)
  checkmate::assert_number(max_wait, lower = 1, upper = 1814400) # 21 days
  scheduler <- tolower(scheduler) # ignore case
  checkmate::assert_subset(scheduler, c("torque", "qsub", "slurm", "sbatch", "sh", "local"))

  job_complete <- FALSE
  wait_start <- Sys.time()

  get_job_status <- function() { # use variables in parent environment
    if (scheduler %in% c("slurm", "sbatch")) {
      status <- slurm_job_status(job_ids)
      state <- sapply(status$State, function(x) {
        switch(x,
          "BOOT_FAIL" = "failed",
          "CANCELLED" = "cancelled",
          "COMPLETED" = "complete",
          "DEADLINE" = "failed",
          "FAILED" = "failed",
          "NODE_FAIL" = "failed",
          "OUT_OF_MEMORY" = "failed",
          "PENDING" = "queued",
          "PREEMPTED" = "failed",
          "RUNNING" = "running",
          "REQUEUED" = "queued",
          "REVOKED" = "failed",
          "SUSPENDED" = "suspended",
          "TIMEOUT" = "failed",
          "MISSING" = "missing" # scheduler has not registered the job
        )
      })
    } else if (scheduler %in% c("sh", "local")) {
      status <- local_job_status(job_ids)
      state <- sapply(status$STAT, function(x) {
        switch(x,
          "C" = "complete",
          "I" = "running", # idle/sleeping
          "R" = "running",
          "S" = "running", # sleeping
          "T" = "suspended",
          "U" = "running",
          "Z" = "failed", # zombie
          stop("Unable to understand job state: ", x)
        )
      })
    } else if (scheduler %in% c("torque", "qsub")) {
      # QSUB
      status <- system(paste("qstat -f", jj, "| grep -i 'job_state'"), intern = TRUE)
      job_state <- sub(".*job_state = ([A-z]).*", "\\1", status, perl = TRUE)
      state <- switch(status,
        "C" = "complete",
        "R" = "running",
        "Q" = "queued",
        "H" = "suspended",
        "other"
      )
    } else {
      stop("unknown scheduler: ", scheduler)
    }
  }

  while (job_complete == FALSE) {
    status <- get_job_status()

    # update wait time
    wait_total <- difftime(Sys.time(), wait_start, units = "sec")

    # Debugging
    # cat("Wait so far: ", wait_total, "\n")

    if (any(status == "running")) {
      if (isFALSE(quiet)) {
        cat("Job(s) still running: ", paste(job_ids[status == "running"], collapse = ","), "\n")
      }
    }

    if (any(status == "queued")) {
      if (isFALSE(quiet)) {
        cat("Job(s) still queued: ", paste(job_ids[status == "queued"], collapse = ","), "\n")
      }
    }

    if (any(status == "suspended")) {
      if (isFALSE(quiet)) {
        cat("Job(s) suspended: ", paste(job_ids[status == "suspended"], collapse = ","), "\n")
      }
    }

    if (any(status == "missing")) {
      if (isFALSE(quiet)) {
        cat("Job(s) missing from scheduler response: ", paste(job_ids[status == "missing"], collapse = ","), "\n")
      }
    }

    if (wait_total > max_wait) {
      if (isTRUE(stop_on_timeout)) {
        stop("Maximum wait time: ", max_wait, " exceeded. Stopping execution of parent script because something is wrong.")
      } else {
        return(FALSE)
      }
    } else if (all(status == "complete")) {
      job_complete <- TRUE # drop out of this loop
      if (isFALSE(quiet)) {
        cat("All jobs have finished.\n")
      }
    } else {
      Sys.sleep(repolling_interval) # wait and repoll jobs
    }
  }

  return(invisible(TRUE))
}

# calls sacct with a job list
slurm_job_status <- function(job_ids = NULL, user = NULL,
                             sacct_format = "jobid,submit,timelimit,start,end,state") {
  if (!is.null(job_ids)) {
    jstring <- paste0("-j", paste(job_ids, collapse = ","))
  } else {
    jstring <- ""
  }

  if (!is.null(user)) {
    ustring <- paste("-u", paste(user, collapse = ","))
  } else {
    ustring <- ""
  }

  # -P specifies a parsable output separated by pipes
  # -X avoids printing subsidiary jobs within each job id
  #cmd <- paste("sacct", jstring, ustring, "-X -P -o", sacct_format)
  cmd <- paste(jstring, ustring, "-X -P -o", sacct_format)
  #cat(cmd, "\n")
  res <- system2("sacct", args=cmd, stdout=TRUE)
  
  df_base <- data.frame(JobID = job_ids)
  df_empty <- df_base %>%
    mutate(
      Submit = NA_character_,
      Timelimit = NA_character_,
      Start = NA_character_,
      End = NA_character_,
      State = "MISSING"
    )

  # handle non-zero exit status -- return empty data
  if (!is.null(attr(res, "status"))) {
    return(df_empty)
  }
  
  if (length(res) == 1L) {
    # data.table::fread will break down (see Github issue )
    out <- readr::read_delim(I(res), delim = "|", show_col_types = FALSE)
  } else {
    out <- data.table::fread(text = res, data.table=FALSE)
  }

  if (!checkmate::test_subset(c("JobID", "State"), names(out))) {
    warning("Missing columns in sacct output")
    return(df_empty)
  }

  out$JobID <- as.character(out$JobID)
  df <- df_base %>%
    dplyr::left_join(out, by = "JobID") %>%
    mutate(State = if_else(is.na(State), "MISSING", State))

  return(df)
}

local_job_status <- function(job_ids = NULL, user = NULL,
                             ps_format = "user,pid,state,time,etime,%cpu,%mem,comm,xstat") {
  job_ids <- type.convert(job_ids, as.is = T) # convert to integers
  checkmate::assert_integerish(job_ids)

  if (!is.null(job_ids)) {
    jstring <- paste("-p", paste(job_ids, collapse = ","))
  } else {
    jstring <- ""
  }

  if (!is.null(user)) {
    ustring <- paste("-u", paste(user, collapse = ","))
  } else {
    ustring <- ""
  }

  # cat(paste("ps", jstring, ustring, "-o", ps_format), sep = "\n")
  res <- suppressWarnings(system2("ps", args = paste(jstring, ustring, "-o", ps_format), stdout = TRUE)) # intern=TRUE)

  # need to trap res of length 1 (just header row) to avoid data.table bug.
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0) {
    hrow <- strsplit(res, "\\s+")[[1]]
    dt <- data.frame(matrix(NA, nrow = length(job_ids), ncol = length(hrow)))
    names(dt) <- hrow
    dt$PID <- as.integer(job_ids)
  } else {
    stopifnot(length(res) > 1)
    # fread and any other parsing can break down with consecutive spaces in body of output.
    # This happens with lstart and start, avoid these for now.
    #header <- gregexpr("\\b", res[1], perl = T)
    #l2 <- gregexpr("\\b", res[2], perl=T)
    dt <- data.table::fread(text = res)
  }

  # fix difference in column naming between FreeBSD and *nux (make all like FreeBSD)
  data.table::setnames(dt, c("S", "COMMAND"), c("STAT", "COMM"), skip_absent = TRUE)

  # build df that fills in missing jobs (completed/killed)
  all_dt <- data.frame(PID = as.integer(job_ids)) %>%
    dplyr::full_join(dt, by = "PID") %>%
    mutate(STAT = substr(STAT, 1, 1)) # only care about first character of state

  all_dt$STAT[is.na(all_dt$STAT)] <- "C" # complete

  return(all_dt)
}