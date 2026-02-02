#' Submit SPM jobs to the cluster for unattended MATLAB/Octave execution
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#' @param level level of analysis (1 or 3)
#' @param model_names optional model names used to subset SPM runs
#' @param rerun logical indicating whether to rerun existing directories
#' @param wait_for optional parent job ids that should complete before these jobs commence
#' @return a vector of job ids for all scripts that were submitted
#' @importFrom glue glue
#' @export
run_spm_sepjobs <- function(gpa, level = 1L, model_names = NULL, rerun = FALSE, wait_for = "") {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, len = 1L)
  checkmate::assert_character(model_names, null.ok = TRUE)
  checkmate::assert_logical(rerun)

  lg <- lgr::get_logger(paste0("glm_pipeline/spm_l", level, "_estimation"))
  lg$set_threshold(gpa$lgr_threshold)
  upd_job_status_path <- system.file("bin/upd_job_status.R", package = "fmri.pipeline")

  if (level == 1L) {
    if (!checkmate::test_class(gpa$l1_model_setup, "l1_setup")) {
      lg$error("In run_spm_sepjobs, did not find an l1_setup object in gpa$l1_model_setup.")
      lg$error("Make sure to run setup_l1_models before running run_spm_sepjobs")
      return(NULL)
    }

    spm_queue <- gpa$l1_model_setup$spm
    if (!is.null(model_names[1L])) {
      spm_queue <- spm_queue %>% dplyr::filter(l1_model %in% !!model_names)
    }

    spm_time <- gpa$parallel$spm$l1_spm_time
    spm_memgb <- gpa$parallel$spm$l1_spm_memgb
    spm_cpus <- gpa$parallel$spm$l1_spm_cpus_per_job
    runsperproc <- gpa$parallel$spm$l1_spm_runs_per_cpu

    run_setup <- isTRUE(gpa$glm_settings$spm$run_l1_setup)
    run_glm <- isTRUE(gpa$glm_settings$spm$run_l1_glm)
    run_contrasts <- isTRUE(gpa$glm_settings$spm$run_l1_contrasts)

    script_map <- c(
      setup = "setup_glm_design.m",
      glm = "run_glm.m",
      contrasts = "estimate_glm_contrasts.m"
    )

    spm_output_directory <- file.path(gpa$output_locations$scheduler_scripts, "spm_l1")
  } else if (level == 3L) {
    if (!checkmate::test_class(gpa$l3_model_setup, "l3_setup")) {
      lg$error("In run_spm_sepjobs, did not find a l3_setup object in gpa$l3_model_setup.")
      lg$error("Make sure to run setup_l3_models before running run_spm_sepjobs")
      return(NULL)
    }

    spm_queue <- gpa$l3_model_setup$spm
    if (!is.null(model_names[1L])) {
      spm_queue <- spm_queue %>% dplyr::filter(l3_model %in% !!model_names)
    }

    spm_time <- gpa$parallel$spm$l3_spm_time
    spm_memgb <- gpa$parallel$spm$l3_spm_memgb
    spm_cpus <- gpa$parallel$spm$l3_spm_cpus_per_job
    runsperproc <- gpa$parallel$spm$l3_spm_runs_per_cpu

    run_setup <- isTRUE(gpa$glm_settings$spm$run_l3_setup)
    run_glm <- isTRUE(gpa$glm_settings$spm$run_l3_glm)
    run_contrasts <- isTRUE(gpa$glm_settings$spm$run_l3_contrasts)

    script_map <- c(
      setup = "setup_l3_design.m",
      glm = "run_l3_glm.m",
      contrasts = "estimate_l3_contrasts.m"
    )

    spm_output_directory <- file.path(gpa$output_locations$scheduler_scripts, "spm_l3")
  }

  if (is.null(spm_queue) || nrow(spm_queue) == 0L) {
    lg$warn("No SPM jobs found for level %d.", level)
    return(NULL)
  }

  if (!"spm_dir" %in% names(spm_queue)) {
    stop("SPM queue is missing spm_dir column. Ensure setup_l1_models/setup_l3_models ran successfully.")
  }

  spm_job_df <- spm_queue %>%
    dplyr::select(spm_dir, spm_complete, to_run)

  if (isTRUE(rerun)) {
    lg$info("rerun = TRUE in run_spm_sepjobs. All SPM directories will be marked for job execution.")
    spm_job_df$to_run <- TRUE
  }

  to_run_dirs <- spm_job_df %>%
    dplyr::filter(to_run == TRUE) %>%
    dplyr::pull(spm_dir)

  if (length(to_run_dirs) == 0L) {
    lg$warn("No SPM level %d directories to execute.", level)
    return(NULL)
  }

  if (!file.exists(spm_output_directory)) {
    lg$debug("Creating spm_l%d working directory: %s", level, spm_output_directory)
    dir.create(spm_output_directory, recursive = TRUE)
  }

  if (is.null(gpa$parallel$compute_environment$spm) || length(gpa$parallel$compute_environment$spm) == 0L) {
    lg$warn("No SPM compute environment configured. MATLAB/Octave may not be available on compute nodes.")
  }

  matlab_cmd <- gpa$glm_settings$spm$matlab_cmd
  if (is.null(matlab_cmd) || !nzchar(matlab_cmd)) matlab_cmd <- "matlab"
  matlab_args <- gpa$glm_settings$spm$matlab_args
  if (is.null(matlab_args) || !nzchar(matlab_args)) matlab_args <- "-nodisplay -nosplash -r"
  matlab_exit <- gpa$glm_settings$spm$matlab_exit
  if (is.null(matlab_exit)) matlab_exit <- "exit;"

  build_matlab_call <- function(script_path) {
    cmd_str <- paste0("run('", script_path, "');", matlab_exit)
    paste(matlab_cmd, matlab_args, shQuote(cmd_str))
  }

  # Scheduler header
  if (gpa$scheduler == "slurm") {
    file_suffix <- ".sbatch"
    preamble <- c(
      "#!/bin/bash",
      "#SBATCH -N 1",
      paste0("#SBATCH -n ", spm_cpus),
      paste0("#SBATCH --time=", spm_time),
      paste0("#SBATCH --mem-per-cpu=", spm_memgb, "G"),
      ifelse(wait_for != "", paste0("#SBATCH --dependency=afterok:", paste(wait_for, collapse=":")), ""),
      sched_args_to_header(gpa),
      "",
      get_compute_environment(gpa, c("spm", "r")),
      "",
      "cd $SLURM_SUBMIT_DIR"
    )
  } else if (gpa$scheduler == "torque") {
    file_suffix <- ".pbs"
    preamble <- c(
      "#!/bin/bash",
      paste0("#PBS -l nodes=1:ppn=", spm_cpus),
      paste0("#PBS -l pmem=", spm_memgb, "gb"),
      ifelse(wait_for != "", paste0("#PBS -W depend=afterok:", paste(wait_for, collapse=":")), ""),
      paste0("#PBS -l walltime=", spm_time),
      sched_args_to_header(gpa),
      "",
      get_compute_environment(gpa, c("spm", "r")),
      "",
      "cd $PBS_O_WORKDIR"
    )
  } else {
    file_suffix <- ".sh"
    preamble <- c(
      "#!/bin/bash",
      get_compute_environment(gpa, c("spm", "r"))
    )
  }

  njobs <- ceiling(length(to_run_dirs) / (spm_cpus * runsperproc))
  df <- data.frame(
    spm_dir = to_run_dirs,
    job = rep(1:njobs, each = spm_cpus * runsperproc, length.out = length(to_run_dirs)),
    stringsAsFactors = FALSE
  )

  submission_id <- basename(tempfile(pattern = "job"))
  joblist <- rep(NA_character_, njobs)

  tracking_sqlite_db <- gpa$output_locations$sqlite_db
  for (j in seq_len(njobs)) {
    outfile <- file.path(spm_output_directory, paste0("spmsep_l", level, "_", j, "_", submission_id, file_suffix))
    cat(preamble, file = outfile, sep = "\n")

    thisrun <- with(df, spm_dir[job == j])
    cat(
      "",
      "job_failed=0",
      paste0("matlab_cmd=", shQuote(matlab_cmd)),
      paste0("matlab_args=", shQuote(matlab_args)),
      paste0("matlab_exit=", shQuote(matlab_exit)),
      paste0("spm_setup_script=", shQuote(script_map[["setup"]])),
      paste0("spm_glm_script=", shQuote(script_map[["glm"]])),
      paste0("spm_contrast_script=", shQuote(script_map[["contrasts"]])),
      "",
      "function spm_runner() {",
      "  odir=\"$1\"",
      "  [ -f \"${odir}/.spm_fail\" ] && rm -f \"${odir}/.spm_fail\"",
      "  if [ -f \"${odir}/.spm_complete\" ]; then",
      "    rm -f \"${odir}/.spm_complete\"",
      "  fi",
      "  start_time=$( date )",
      "  exit_code=0",
      "",
      "  run_matlab() {",
      "    script_path=\"$1\"",
      "    cmd_str=\"run('${script_path}');${matlab_exit}\"",
      "    ${matlab_cmd} ${matlab_args} \"${cmd_str}\"",
      "    cmd_exit=$?",
      "    if [ $cmd_exit -ne 0 ]; then",
      "      exit_code=$cmd_exit",
      "    fi",
      "  }",
      "",
      "  gunzip_script=\"${odir}/gunzip_commands.sh\"",
      "  if [ -f \"${gunzip_script}\" ]; then",
      "    bash \"${gunzip_script}\"",
      "    if [ $? -ne 0 ]; then exit_code=$?; fi",
      "  fi",
      "",
      if (isTRUE(run_setup)) {
        c(
          "  setup_script=\"${odir}/${spm_setup_script}\"",
          "  if [ -f \"${setup_script}\" ]; then",
          "    run_matlab \"${setup_script}\"",
          "  fi"
        )
      } else {
        character(0)
      },
      if (isTRUE(run_glm)) {
        c(
          "  glm_script=\"${odir}/${spm_glm_script}\"",
          "  if [ -f \"${glm_script}\" ]; then",
          "    run_matlab \"${glm_script}\"",
          "  fi"
        )
      } else {
        character(0)
      },
      "",
      "  contrast_setup_script=\"${odir}/setup_spm_contrasts.sh\"",
      "  if [ -f \"${contrast_setup_script}\" ]; then",
      "    bash \"${contrast_setup_script}\"",
      "    if [ $? -ne 0 ]; then exit_code=$?; fi",
      "  fi",
      "",
      if (isTRUE(run_contrasts)) {
        c(
          "  contrast_script=\"${odir}/${spm_contrast_script}\"",
          "  if [ -f \"${contrast_script}\" ]; then",
          "    run_matlab \"${contrast_script}\"",
          "  fi"
        )
      } else {
        character(0)
      },
      "",
      "  end_time=$( date )",
      "  if [ $exit_code -eq 0 ]; then",
      "    status_file=\"${odir}/.spm_complete\"",
      "  else",
      "    status_file=\"${odir}/.spm_fail\"",
      "  fi",
      "  echo $start_time > \"${status_file}\"",
      "  echo $end_time >> \"${status_file}\"",
      "  return $exit_code",
      "}",
      "",
      "function spm_killed() {",
      "  kill_time=$( date )",
      "  if [ -n \"${odir}\" ]; then",
      "    echo $kill_time > \"${odir}/.spm_fail\"",
      "  fi",
      paste("  Rscript", upd_job_status_path, "--job_id", "'$job_id'", "--sqlite_db", tracking_sqlite_db, "--status", "FAILED"),
      "  exit 1",
      "}",
      "trap spm_killed SIGTERM",
      paste("Rscript", upd_job_status_path, "--job_id", "'$job_id'", "--sqlite_db", tracking_sqlite_db, "--status", "STARTED"),
      "",
      sep = "\n",
      file = outfile,
      append = TRUE
    )

    for (spm_dir in thisrun) {
      if (!dir.exists(spm_dir)) {
        lg$warn("Skipping missing spm_dir: %s", spm_dir)
        next
      }
      cat(
        "spm_runner ", shQuote(spm_dir), "\n",
        "if [ $? -ne 0 ]; then job_failed=1; fi",
        sep = "",
        file = outfile,
        append = TRUE
      )
    }

    cat(
      "if [ $job_failed -ne 0 ]; then",
      paste("  Rscript", upd_job_status_path, "--job_id", "'$job_id'", "--sqlite_db", tracking_sqlite_db, "--status", "FAILED"),
      "  exit 1",
      "fi",
      paste("Rscript", upd_job_status_path, "--job_id", "'$job_id'", "--sqlite_db", tracking_sqlite_db, "--status", "COMPLETED"),
      "",
      sep = "\n",
      file = outfile,
      append = TRUE
    )

    tracking_args <- list(
      job_name = paste0("spmsep_l", level, "_", j),
      batch_directory = spm_output_directory,
      n_nodes = 1,
      n_cpus = spm_cpus,
      wall_time = spm_time,
      mem_per_cpu = spm_memgb,
      scheduler_options = gpa$parallel$sched_args
    )

    if (wait_for != "") tracking_args$parent_job_id <- wait_for

    joblist[j] <- cluster_job_submit(
      outfile,
      scheduler = gpa$scheduler,
      tracking_sqlite_db = tracking_sqlite_db,
      tracking_args = tracking_args
    )
  }

  writeLines(joblist, con = file.path(spm_output_directory, paste0("spmsep_l", level, "_jobs.txt")))

  return(joblist)
}
