#' function to submit 3dLMEr jobs on a cluster
#'
#' @param gpa a \code{glm_pipeline_arguments} object
#' @param model_names optional subset of L3 model names to run
#' @param rerun logical, whether to rerun existing models
#' @param wait_for job ID(s) to wait for
#' 
#' @return vector of submitted job IDs
#' @export
run_3dlmer_sepjobs <- function(gpa, level=3L, model_names=NULL, rerun=FALSE, wait_for="") {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower=3L, upper=3L) # only L3 for now
  
  lg <- lgr::get_logger("glm_pipeline/run_3dlmer_sepjobs")
  
  if (!checkmate::test_class(gpa$l3_model_setup$afni, "data.frame")) {
    lg$error("No AFNI 3dLMEr setup found in gpa$l3_model_setup$afni.")
    return(NULL)
  }
  
  setup_df <- gpa$l3_model_setup$afni
  if (!is.null(model_names)) {
    setup_df <- setup_df %>% dplyr::filter(l3_model %in% !!model_names)
  }
  
  if (nrow(setup_df) == 0L) {
    lg$warn("No 3dLMEr models found to run.")
    return(NULL)
  }
  
  # Filter based on rerun and status
  if (!isTRUE(rerun)) {
    setup_df <- setup_df %>% dplyr::filter(feat_complete == FALSE)
  }
  
  if (nrow(setup_df) == 0L) {
    lg$info("All requested 3dLMEr models are already complete. Use rerun=TRUE to force.")
    return(NULL)
  }
  
  # Resource settings
  lmer_time <- gpa$parallel$afni$l3_lmer_time %||% "48:00:00"
  lmer_memgb <- gpa$parallel$afni$l3_lmer_memgb %||% "32"
  lmer_cpus <- gpa$parallel$afni$l3_lmer_njobs %||% 8
  
  tracking_sqlite_db <- gpa$output_locations$sqlite_db
  upd_job_status_path <- system.file("bin/upd_job_status.R", package = "fmri.pipeline")
  
  job_ids <- c()
  
  for (i in seq_len(nrow(setup_df))) {
    script_to_run <- setup_df$afni_script[i]
    l3_name <- setup_df$l3_model[i]
    con_name <- paste0(setup_df$l1_cope_name[i], "_", setup_df$l2_cope_name[i])
    
    # We create a wrapper script for scheduling, similar to run_feat_sepjobs
    # but 3dLMEr scripts are already generated. We just need to add the scheduler header.
    
    submission_id <- basename(tempfile(pattern = paste0("3dlmer_", l3_name, "_")))
    sched_script <- file.path(dirname(script_to_run), paste0(submission_id, ".sbatch"))
    
    if (gpa$scheduler == "slurm") {
      header <- c(
        "#!/bin/bash",
        "#SBATCH -N 1",
        paste0("#SBATCH -n ", lmer_cpus),
        paste0("#SBATCH --time=", lmer_time),
        paste0("#SBATCH --mem=", lmer_memgb, "G"),
        ifelse(wait_for != "", paste0("#SBATCH --dependency=afterok:", paste(wait_for, collapse=":")), ""),
        sched_args_to_header(gpa),
        "",
        get_compute_environment(gpa, c("afni", "r")),
        "",
        "cd $SLURM_SUBMIT_DIR",
        "job_id=$SLURM_JOB_ID"
      )
    } else {
      # Torque fallback...
      header <- c(
        "#!/bin/bash",
        paste0("#PBS -l nodes=1:ppn=", lmer_cpus),
        paste0("#PBS -l mem=", lmer_memgb, "gb"),
        paste0("#PBS -l walltime=", lmer_time),
        ifelse(wait_for != "", paste0("#PBS -W depend=afterok:", paste(wait_for, collapse=":")), ""),
        sched_args_to_header(gpa),
        "",
        get_compute_environment(gpa, c("afni", "r")),
        "",
        "cd $PBS_O_WORKDIR",
        "job_id=$PBS_JOBID"
      )
    }
    
    # Add tracking and execution
    script_content <- c(
      header,
      "",
      paste("Rscript", upd_job_status_path, "--job_id" , "\"$job_id\"", "--sqlite_db", tracking_sqlite_db, "--status", "STARTED"),
      "",
      paste("bash", script_to_run),
      "exit_code=$?",
      "",
      "if [ $exit_code -eq 0 ]; then",
      paste("  Rscript", upd_job_status_path, "--job_id" , "\"$job_id\"", "--sqlite_db", tracking_sqlite_db, "--status", "COMPLETED"),
      "else",
      paste("  Rscript", upd_job_status_path, "--job_id" , "\"$job_id\"", "--sqlite_db", tracking_sqlite_db, "--status", "FAILED"),
      "fi",
      "",
      "exit $exit_code"
    )
    
    writeLines(script_content, sched_script)
    
    tracking_args <- list(
      job_name = paste0("3dlmer_", l3_name, "_", con_name),
      batch_directory = dirname(script_to_run),
      n_nodes = 1,
      n_cpus = lmer_cpus,
      wall_time = lmer_time,
      mem_total = lmer_memgb,
      scheduler_options = gpa$parallel$sched_args
    )
    if (wait_for != "") tracking_args$parent_job_id <- wait_for
    
    jid <- cluster_job_submit(sched_script, scheduler = gpa$scheduler, tracking_sqlite_db = tracking_sqlite_db, tracking_args = tracking_args)
    job_ids <- c(job_ids, jid)
  }
  
  return(job_ids)
}
