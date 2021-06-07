#' This function submits a single script to a high-performance cluster using a scheduler (Slurm or TORQUE).
#' It accepts a vector of arguments to be passed to the scheduler and
#' a vector of environment variables that should be passed to the compute node at job execution.
#'
#' The function returns the jobid of the scheduled job.
#'
#' @param script A script that should be executed by the scheduler. This can contain scheduler directives, but in the
#'      case of conflicts, the directives passed with \code{sched_args} will take precedence.
#' @param scheduler Which scheduler to use for job submission. Options are 'qsub', 'torque', 'sbatch', 'slurm', or 'sh'.
#'      The terms 'qsub' and 'torque' are aliases (where 'torque' submits via the qsub command). Likewise for 'sbatch'
#'      and 'slurm'. The scheduler 'sh' does not submit to any scheduler at all, but instead executes the command 
#'      immediately via sh.
#' @param sched_args A character vector of arguments to be included in the scheduling command. On TORQUE, these
#'      will typically begin with '-l' such as '-l walltime=10:00:00'.
#' @param env_variables A named character vector containing environment variables and their values to be passed
#'      to the \code{script} at execution time. This is handled by the '-v' directive on TORQUE clusters and
#'      by '--export' on Slurm clusters. The names of this vector are the environment variable names and
#'      the values of the vector are the environment variable values to be passed in.
#'      If you want to propagate the current value of an environment variable to the compute node at runtime,
#'      use NA as the value of the element in \code{env_variables}. See examples.
#' @param export_all Whether to export all environment variables to the compute node at runtime. Default: FALSE
#' @param echo Whether to echo the job submission command to the terminal at the time it is scheduled. Default: TRUE.
#' @param fail_on_error Whether to stop execution of the script (TRUE), or issue a warning (FALSE) if the job 
#'      submission fails. Defaults to FALSE (i.e., issue a warning).
#'
#' @return A character string containing the jobid of the scheduled job.
#'
#' @examples
#' \dontrun{
#'   #simple PBS submission
#'   cluster_job_submit('myscript.bash', scheduler="torque", sched_args=c('-l walltime=10:00:00', '-l nodes=1:ppn=20'),
#'      env_variables=c(RUN_INDEX=2, MODEL_NAME='FSE21'))
#'
#'   #To forward environment variables without explicitly providing values. Note that these must
#'   #  be in R's system environment (cf. Sys.getenv) at execution time to forward as expected.
#'   cluster_job_submit('myscript.sbatch', scheduler="slurm",
#'      sched_args=c('-p general', '-N 1', '-n 12', '--mem=10g', '-t 02-00:00:00'),
#'      env_variables=c(RUN_INDEX=2, R_HOME=NA, JAVA_HOME=NA))
#' }
#'
#' @author Michael Hallquist
#' @importFrom tools file_path_sans_ext
#' @importFrom checkmate assert_character assert_subset
#' @export
cluster_job_submit <- function(script, scheduler="slurm", sched_args=NULL,
                           env_variables=NULL, export_all=FALSE, echo=TRUE, fail_on_error=FALSE) {

  checkmate::assert_character(scheduler, max.len=1)
  checkmate::assert_subset(scheduler, c("qsub", "torque", "sbatch", "slurm", "sh"))
  checkmate::assert_logical(echo, max.len = 1L)
  checkmate::assert_logical(fail_on_error, max.len=1L)

  if (scheduler == "torque") {
    scheduler <- "qsub" # simpler internal tracking
    if (isTRUE(export_all)) {
      sched_args <- c("-V", sched_args)
    } # directive to export all environment variables to script
  } else if (scheduler == "slurm") {
    scheduler <- "sbatch"
    if (isTRUE(export_all)) {
      env_variables <- c(ALL = NA, env_variables)
    } # directive to export all environment variables to script
  } else if (scheduler == "sh") {
    if (!is.null(sched_args)) {
      message("Omitting scheduler arguments for sh execution")
    }
    sched_args <- NULL # not relevant
  }

  # Scheduler arguments are just pasted together with spaces.
  # Thus, arguments like '--mem=5g' and '-n 12' are not handled differently
  if (!is.null(sched_args)) { sched_args <- paste(sched_args, collapse=" ") }

  #subfunction to handle variable=value and variable combinations
  paste_args <- function(str_vec) {
    nms <- names(str_vec)
    sapply(seq_along(str_vec), function(x) {
      if (is.na(str_vec[x])) {
        return(nms[x]) #just the name of the env variable (forwards from environment)
      } else {
        #force argument to be quoted to avoid problems with spaces
        val <- ifelse(grepl("^[\"'].*[\"']$", str_vec[x], perl=TRUE), str_vec[x], paste0("\"", str_vec[x], "\""))
        return(paste0(nms[x], "=", val))
      }
    })
  }

  #pass through environment variables
  if (!is.null(env_variables)) {
    env_variables <- paste_args(env_variables) #convert to appropriate name-value pairs
    if (scheduler == "qsub") {
      env_variables <- paste("-v", paste(env_variables, collapse=","))
    } else if (scheduler == "sbatch") {
      env_variables <- paste0("--export=", paste(env_variables, collapse=","))
    } else if (scheduler == "sh") {
      env_variables <- paste(sapply(env_variables, function(x) {
        ifelse(grepl("=", x, fixed = TRUE), x, paste0(x, "=\"$", x, "\""))
      }), collapse = " ")
    }

    sched_args <- paste(sched_args, env_variables)
  }

  stopifnot(file.exists(script))
  # use unique temp files to avoid parallel collisions in job tracking
  sub_stdout <- paste0(tempfile(), "_", tools::file_path_sans_ext(basename(script)), "_stdout") 
  sub_stderr <- paste0(tempfile(), "_", tools::file_path_sans_ext(basename(script)), "_stderr")

  if (scheduler == "sh") {
    # for direct execution, need to pass environment variables by prepending
    if (isTRUE(echo)) cat(paste(env_variables, scheduler, script), "\n")
    # submit the job script and return the jobid
    jobres <- system(paste(env_variables, scheduler, script, ">", sub_stdout, "2>", sub_stderr))

  } else {
    if (isTRUE(echo)) cat(paste(scheduler, script, sched_args), "\n")
    # submit the job script and return the jobid
    jobres <- system2(scheduler, args=paste(script, sched_args), stdout=sub_stdout, stderr=sub_stderr)
  }

  jobid <- if (file.exists(sub_stdout)) scan(file = sub_stdout, what = "char", sep = "\n", quiet = TRUE) else ""
  joberr <- if (file.exists(sub_stderr)) {
    paste(scan(file = sub_stderr, what = "char", sep = "\n", quiet = TRUE), collapse = ". ")
  } else {
    ""
  }

  if (jobres != 0) {
    jobid <- NULL
    if (isTRUE(fail_on_error)) {
      stop("Job submission failed: ", script, ", error: ", joberr, ", errcode: ", jobres)
    } else {
      warning("Job submission failed: ", script, ", error: ", joberr, ", errcode: ", jobres)
    }
  } else {
    jobid <- sub("Submitted batch job ", "", jobid, fixed = TRUE) # replace irrelevant details if needed
  }

  return(jobid)
}
