#' function to submit a set of jobs on a cluster to estimate many Feat level 1 models
#' 
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#' @param l1_model_names an optional level 1 model name used to subset the runs to submit to feat.
#'   If not provided, all level 1 models in \code{gpa} will be submitted for feat estimation.
#' @param rerun a logical indicating whether to re-run an existing directory. Default: FALSE
#' @param wait_for a parent job that should complete before these jobs commence
#' @return a vector of job ids for all files that were submitted to the cluster
run_feat_lvl1_sepjobs <- function(gpa, l1_model_names=NULL, rerun=FALSE, wait_for="") {
  ## this version of the FSL LVL1 feat estimation creates multiple qsub scripts in a temporary directory
  ## where each script has a number of feat calls that are forked, then the script waits for completion
  ## This circumvents the ICS limit on multiple nodes in R using a SNOW cluster.
  ## The primary calculation is how many files there are to run relative to processors and files per processor (chunking)


  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok=TRUE)
  checkmate::assert_logical(rerun)

  lg <- lgr::get_logger("glm_pipeline/l1_estimation")

  if (!checkmate::test_class(gpa$l1_model_setup, "l1_setup")) {
    lg$error("In run_feat_lvl1_sepjobs, did not find an l1_setup object in gpa$l1_model_setup.")
    lg$error("Make sure to run setup_l1_models before running run_feat_lvl1_sepjobs")
    return(NULL)
  }

  if (!is.null(l1_model_names[1L])) {
    l1_queue <- gpa$l1_model_setup$fsl %>% dplyr::filter(l1_model %in% !!l1_model_names)
  } else {
    l1_queue <- gpa$l1_model_setup$fsl
  }

  # need this to be cached somewhere...
  l1_working_directory <- file.path(gpa$working_directory, "feat_lvl1")

  cpusperjob <- 8 #number of cpus per qsub
  runsperproc <- 3 #number of feat calls per processor

  fsf_files <- l1_queue$l1_feat_fsf

  #figure out which fsf files have already been run
  dir_expect <- gsub("\\.fsf$", ".feat", fsf_files, perl=TRUE)

  to_run <- c()

  for (f in seq_along(fsf_files)) {
    if (file.exists(dir_expect[f])) {
      if (rerun) {
        cmd <- paste0("rm -rf \"", dir_expect[f], "\"")
        lg$info("Removing old directory: %s", cmd)
        system(cmd)
        to_run <- c(to_run, fsf_files[f]) #add to queue
      } else {
        lg$info("Skipping existing directory: %s", dir_expect[f])
      }
    } else {
      to_run <- c(to_run, fsf_files[f])
    }
  }

  if (length(to_run) == 0L) {
    lg$warn("No LVL1 .fsf files to execute.")
    return(NULL)
  }

  lg$info("About to run the following fsf files in parallel:")
  lg$info("  File: %s", to_run)

  if (gpa$scheduler == "slurm") {
    file_suffix <- ".sbatch"
    preamble <- c(
      "#!/bin/bash",
      "#SBATCH -p general",
      "#SBATCH -N 1", #always single node for now
      paste0("#SBATCH -n ", cpusperjob),
      paste0("#SBATCH --time=", gpa$parallel$fsl$l1_feat_time),
      paste0("#SBATCH --mem-per-cpu=", gpa$parallel$fsl$l1_feat_memgb, "G"),
      "",
      "",
      gpa$parallel$fsl$compute_environment
    )

  } else if (gpa$scheduler == "torque") {
    file_suffix <- ".pbs"
    preamble <- c(
      "#!/bin/bash",
      "#PBS -A mnh5174_c_g_sc_default",
      paste0("#PBS -l nodes=1:ppn=", cpusperjob),
      paste0("#PBS -l pmem=", gpa$parallel$fsl$l1_feat_memgb, "gb"),
      ifelse(wait_for != "", paste0("#PBS -W depend=afterok:", wait_for), ""), # allow job dependency on upstream setup
      paste0("#PBS -l walltime=", gpa$parallel$fsl$l1_feat_time),
      "#PBS -j oe",
      "#PBS -m n",
      "#PBS -W group_list=mnh5174_collab",
      "",
      "",
      gpa$parallel$fsl$compute_environment,
      "",
      "cd $PBS_O_WORKDIR"
    )
  }

  if (!file.exists(l1_working_directory)) {
    lg$debug("Creating l1 working directory: %s", l1_working_directory)
    dir.create(l1_working_directory, recursive = TRUE)
  }

  njobs <- ceiling(length(to_run)/(cpusperjob*runsperproc))

  #use length.out on rep to ensure that the vectors align even if chunks are uneven wrt files to run
  df <- data.frame(fsf=to_run, job=rep(1:njobs, each=cpusperjob*runsperproc, length.out=length(to_run)), stringsAsFactors=FALSE)
  df <- df[order(df$job), ]

  joblist <- rep(NA_character_, njobs)
  for (j in 1:njobs) {
    outfile <- paste0(l1_working_directory, "/featsep_lvl1_", j, "_", basename(tempfile()), file_suffix)
    cat(preamble, file=outfile, sep="\n")
    thisrun <- with(df, fsf[job==j])
    cat(
      "function feat_runner() {",
      "  local odir=\"${1/.fsf/.feat}\"",
      "  date > ${odir}/.feat_complete && feat $1 && date >> ${odir}/.feat_complete",
      "}",
      sep = "\n", file = outfile, append = TRUE
    )
    cat(paste("feat_runner", thisrun, "&"), file=outfile, sep="\n", append=TRUE)
    cat("wait\n\n", file=outfile, append=TRUE)
    cat(paste(
      "bash", file.path(gpa$pipeline_home, "inst", "bash", "gen_feat_reg_dir"),
      unique(dirname(thisrun))
    ), sep = "\n", file = outfile, append = TRUE)
    joblist[j] <- cluster_job_submit(outfile)
    #joblist[j] <- "dummy"
  }

  #write the list of separate feat qsub jobs that are now queued (so that LVL2 can wait on these)
  #should also return this to the caller as a function?
  writeLines(joblist, con=file.path(l1_working_directory, "sep_lvl1_jobs.txt"))

  return(joblist)
}
