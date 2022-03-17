#' function to submit a set of jobs on a cluster to estimate many Feat level 1 models
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#' @param model_names an optional level 1 model name used to subset the runs to submit to feat.
#'   If not provided, all level 1 models in \code{gpa} will be submitted for feat estimation.
#' @param rerun a logical indicating whether to re-run an existing directory. Default: FALSE
#' @param wait_for a parent job that should complete before these jobs commence
#' @return a vector of job ids for all files that were submitted to the cluster
#' @importFrom glue glue
#' @export
run_feat_sepjobs <- function(gpa, level=1L, model_names=NULL, rerun=FALSE, wait_for="") {
  # this version of the FSL LVL1 feat estimation creates multiple qsub scripts in a temporary directory
  # where each script has a number of feat calls that are forked, then the script waits for completion
  ## This circumvents the ICS limit on multiple nodes in R using a SNOW cluster.
  ## The primary calculation is how many files there are to run relative to processors and files per processor (chunking)

  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower=1L, upper=3L, len=1L)
  checkmate::assert_character(model_names, null.ok=TRUE)
  checkmate::assert_logical(rerun)

  lg <- lgr::get_logger(paste0("glm_pipeline/l", level, "_estimation"))
  lg$set_threshold(gpa$lgr_threshold)

  # TODO: support named subsetting in model_names, like list(l1_model="pe_only", l2_model=c("l2_1", "l2_2"))
  # level-specific copy-paste is a bit clunky here, but at least it's clear...
  if (level == 1) {
    if (!checkmate::test_class(gpa$l1_model_setup, "l1_setup")) {
      lg$error("In run_feat_sepjobs, did not find an l1_setup object in gpa$l1_model_setup.")
      lg$error("Make sure to run setup_l1_models before running run_feat_sepjobs")
      return(NULL)
    }

    if (!is.null(model_names[1L])) {
      feat_queue <- gpa$l1_model_setup$fsl %>% dplyr::filter(l1_model %in% !!model_names)
    } else {
      feat_queue <- gpa$l1_model_setup$fsl
    }

    feat_time   <- gpa$parallel$fsl$l1_feat_time
    feat_memgb  <- gpa$parallel$fsl$l1_feat_memgb
    feat_cpus   <- gpa$parallel$fsl$l1_feat_cpus_per_job # number of cpus per scheduled job
    runsperproc <- gpa$parallel$fsl$l1_feat_runs_per_cpu # number of feat calls per processor
  } else if (level == 2) {
    if (!checkmate::test_class(gpa$l2_model_setup, "l2_setup")) {
      lg$error("In run_feat_sepjobs, did not find an l2_setup object in gpa$l2_model_setup.")
      lg$error("Make sure to run setup_l2_models before running run_feat_sepjobs")
      return(NULL)
    }

    if (!is.null(model_names[1L])) {
      feat_queue <- gpa$l2_model_setup$fsl %>% dplyr::filter(l2_model %in% !!model_names)
    } else {
      feat_queue <- gpa$l2_model_setup$fsl
    }

    feat_time <- gpa$parallel$fsl$l2_feat_time
    feat_memgb <- gpa$parallel$fsl$l2_feat_memgb
    feat_cpus <- 8 # TODO: populate parallel settings for l2
    runsperproc <- 2 # number of feat calls per processor
  } else if (level == 3) {
    if (!checkmate::test_class(gpa$l3_model_setup, "l3_setup")) {
      lg$error("In run_feat_sepjobs, did not find an gpa$l3_model_setup object.")
      lg$error("Make sure to run setup_l3_models before running run_feat_sepjobs")
      return(NULL)
    }

    if (!is.null(model_names[1L])) {
      feat_queue <- gpa$l3_model_setup$fsl %>% dplyr::filter(l3_model %in% !!model_names)
    } else {
      feat_queue <- gpa$l3_model_setup$fsl
    }

    feat_time <- gpa$parallel$fsl$l3_feat_time
    # memory is for total job, not per cpu at l3
    feat_memgb <- ceiling(as.numeric(gpa$parallel$fsl$l3_feat_memgb) / as.numeric(gpa$parallel$fsl$l3_feat_cpusperjob))
    feat_cpus <- gpa$parallel$fsl$l3_feat_cpusperjob
    runsperproc <- 1 # number of feat calls per processor
  }

  fail_action <- gpa$glm_settings$fsl[[glue("failed_l{level}_folder_action")]]
  incomplete_action <- gpa$glm_settings$fsl[[glue("incomplete_l{level}_folder_action")]]
  if (is.null(fail_action)) fail_action <- "delete" #should be populated in finalize step, but just in case
  if (is.null(incomplete_action)) incomplete_action <- "delete"

  if (!"feat_fsf" %in% names(feat_queue)) {
    stop(paste0(
      "Fatal error with queue of feat jobs to be run. There is no feat_fsf column in gpa$l",
      level, "_model_setup$fsl. Make sure that the setup step has been run successfully first."
    ))
  }

  feat_job_df <- feat_queue %>%
    dplyr::select(feat_fsf, feat_dir, feat_complete, feat_failed, to_run)

  # location of scheduler scripts
  feat_output_directory <- file.path(gpa$output_locations$scheduler_scripts, paste0("feat_l", level))

  # TODO: probably use a jobs | wc -l approach to throttling jobs within a submission
  # and we need to make this more flexible

  handle_delete_archive_keep <- function(df, action="delete", type="failed") {
    if (nrow(df) > 0L) {
      if (action == "delete") {
        for (ff in df$feat_dir) {
          if (dir.exists(ff)) {
            cmd <- glue("rm -rf \"{ff}\"")
            lg$info("Removing %s l%d directory: %s", type, level, cmd)
            system(cmd)
          }
        }
      } else if (action == "archive") {
        for (ff in df$feat_dir) {
          if (dir.exists(ff)) {
            new_dir <- sub("(.*)(\\.g?feat)$", paste0("\\1-", format(Sys.time(), "%F-%H%M%S"), "\\2"), ff, perl = TRUE)
            cmd <- glue("mv \"{ff}\" \"{new_dir}\"")
            lg$info("Archiving %s l%d directory: %s", type, level, cmd)
            system(cmd)
          }
        }
      } else { # "keep"
        lg$info("Keeping %s l%d directory: %s", type, level, df$feat_dir)
      }
    }
  }

  feat_failed <- feat_job_df %>% dplyr::filter(feat_failed == TRUE)
  feat_incomplete <- feat_job_df %>% dplyr::filter(feat_complete == FALSE)

  # if a directory is both failed and incomplete, it will be deleted under the failed set
  handle_delete_archive_keep(feat_failed, action=fail_action, type="failed")
  handle_delete_archive_keep(feat_incomplete, action=incomplete_action, type="incomplete")

  if (isTRUE(rerun)) {
    lg$info("rerun = TRUE in run_feat_sepjobs. All fsfs will be marked for job execution.")
    feat_job_df$to_run <- TRUE
  }

  # TODO: keep this as a data.frame and return an amended lXX_model_setup to the caller that includes the job id and batch script
  to_run <- feat_job_df %>%
    dplyr::filter(to_run == TRUE) %>%
    dplyr::pull(feat_fsf)

  if (length(to_run) == 0L) {
    lg$warn("No Feat level %d .fsf files to execute.", level)
    return(NULL)
  }

  lg$info("About to run the following fsf files in parallel:")
  lg$info("  File: %s", to_run)

  if (gpa$scheduler == "slurm") {
    file_suffix <- ".sbatch"
    preamble <- c(
      "#!/bin/bash",
      "#SBATCH -N 1", #always single node for now
      paste0("#SBATCH -n ", feat_cpus),
      paste0("#SBATCH --time=", feat_time),
      paste0("#SBATCH --mem-per-cpu=", feat_memgb, "G"),
      ifelse(wait_for != "", paste0("#SBATCH --dependency=afterok:", paste(wait_for, collapse=":")), ""), # allow job dependency on upstream setup
      sched_args_to_header(gpa), # analysis-level SBATCH directives
      "",
      "",
      gpa$parallel$fsl$compute_environment,
      "",
      "cd $SLURM_SUBMIT_DIR"
    )

  } else if (gpa$scheduler == "torque") {
    file_suffix <- ".pbs"
    preamble <- c(
      "#!/bin/bash",
      paste0("#PBS -l nodes=1:ppn=", feat_cpus),
      paste0("#PBS -l pmem=", feat_memgb, "gb"),
      ifelse(wait_for != "", paste0("#PBS -W depend=afterok:", paste(wait_for, collapse=":")), ""), # allow job dependency on upstream setup
      paste0("#PBS -l walltime=", feat_time),
      sched_args_to_header(gpa), # analysis-level PBS directives
      "",
      "",
      gpa$parallel$fsl$compute_environment,
      "",
      "cd $PBS_O_WORKDIR"
    )
  }

  if (!file.exists(feat_output_directory)) {
    lg$debug("Creating l%d working directory: %s", level, feat_output_directory)
    dir.create(feat_output_directory, recursive = TRUE)
  }

  if (level == 3) {
    njobs <- length(to_run) # parallel across slices, one job per model
    feat_binary <- system.file("bin/feat_parallel", package = "fmri.pipeline")
    stopifnot(file.exists(feat_binary))
    df <- data.frame(fsf = to_run, job = 1:njobs, stringsAsFactors = FALSE)
  } else {
    njobs <- ceiling(length(to_run) / (feat_cpus * runsperproc))
    feat_binary <- "feat"
    # use length.out on rep to ensure that the vectors align even if chunks are uneven wrt files to run
    df <- data.frame(
      fsf = to_run,
      job = rep(1:njobs, each = feat_cpus * runsperproc, length.out = length(to_run)), stringsAsFactors = FALSE
    )
  }

  df <- df[order(df$job), ]

  submission_id <- basename(tempfile(pattern = "job"))
  joblist <- rep(NA_character_, njobs)
  for (j in seq_len(njobs)) {
    outfile <- paste0(feat_output_directory, "/featsep_l", level, "_", j, "_", submission_id, file_suffix)
    cat(preamble, file=outfile, sep="\n")
    thisrun <- with(df, fsf[job==j])
    cat(
      "function feat_runner() {",
      ifelse(level == 1L, "  local odir=\"${1/.fsf/.feat}\"", "  local odir=\"${1/.fsf/.gfeat}\""),
      "  [ -f \"${odir}/.feat_fail\" ] && rm -f \"${odir}/.feat_fail\"",
      "  start_time=$( date )",
      "  if [ $# -eq 2 ]; then",
      paste0("    ", feat_binary, " $1 -P $2"),
      "  else",
      paste0("    ", feat_binary, " $1"),
      "  fi",
      "  exit_code=$?",
      "  end_time=$( date )",
      "  if [ $exit_code -eq 0 ]; then",
      "    status_file=\"${odir}/.feat_complete\"",
      "  else",
      "    status_file=\"${odir}/.feat_fail\"",
      "  fi",
      "  echo $start_time > \"${status_file}\"",
      "  echo $end_time >> \"${status_file}\"",
      "}",
      "function feat_killed() {",
      "  kill_time=$( date )",
      "  echo $kill_time > \"${odir}/.feat_fail\"",
      "}",
      "trap feat_killed SIGTERM",
      sep = "\n", file = outfile, append = TRUE
    )
    if (level == 3L) {
      cat(paste("feat_runner", thisrun, feat_cpus, "&"), file = outfile, sep = "\n", append = TRUE)
    } else {
      cat(paste("feat_runner", thisrun, "&"), file=outfile, sep="\n", append=TRUE)
    }
    cat("wait\n\n", file=outfile, append=TRUE)
    if (level != 3L) {
      cat(paste(
        "bash", system.file("bash/gen_feat_reg_dir", package = "fmri.pipeline"),
        unique(dirname(thisrun))
      ), sep = "\n", file = outfile, append = TRUE)
    }
    joblist[j] <- cluster_job_submit(outfile, scheduler=gpa$scheduler)
    #joblist[j] <- "dummy"
  }

  # write the list of separate feat qsub jobs that are now queued (so that LVL2 can wait on these)
  # should also return this to the caller as a function?
  writeLines(joblist, con = file.path(feat_output_directory, paste0("sep_l", level, "_jobs.txt")))

  return(joblist)
}

