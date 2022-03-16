# TODO: Support ability not to run models at a certain level (e.g., level 3) if only lower-level models are needed, or also
# run just higher-level models if lower levels are already complete.

#' primary function for running a GLM analysis pipeline
#' 
#' @param gpa a glm_pipeline_arguments object containing a model specification (created by setup_glm_pipeline)
#' @param l1_model_names a character vector of level 1 model names (specified during build_l1_models) that should be executed
#' @param l2_model_names a character vector of level 2 model names (specified during build_l2_models) that should be executed
#' @param l3_model_names a character vector of level 3 model names (specified during build_l3_models) that should be executed
#' @param glm_software which glm software should be used for model estimation (not implemented yet)
#' @importFrom checkmate assert_string assert_class assert_subset assert_integerish
#' @export
run_glm_pipeline <- function(gpa, l1_model_names = "prompt", l2_model_names = "prompt",
l3_model_names = "prompt", glm_software = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l3_model_names, null.ok = TRUE)

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")
  lg$set_threshold(gpa$lgr_threshold)

  model_list <- choose_glm_set(gpa, l1_model_names, l2_model_names, l3_model_names, lg)
  if (is.null(model_list)) { return(invisible(NULL)) } # user canceled

  batch_id <- uuid::UUIDgenerate()
  if (is.null(gpa$batch_run)) gpa$batch_run <- list()

  # setup a unique directory for this GLM batch submission
  # this will contain all scripts for this run and contain a snapshot of the gpa object
  # that is amended as each batch job runs and completes.
  batch_directory <- file.path(gpa$output_locations$scheduler_scripts, paste0("batch_", batch_id))
  if (!dir.exists(batch_directory)) dir.create(batch_directory, recursive = TRUE)
  gpa_cache <- file.path(batch_directory, "run_pipeline_cache.RData")
  save(gpa, file=gpa_cache)

  # batch job to finalize pipeline configuration
  f_batch <- R_batch_job$new(
    job_name = "finalize_configuration", batch_directory = batch_directory, scheduler = gpa$scheduler,
    input_rdata_file = gpa_cache, output_rdata_file = gpa_cache,
    n_nodes = 1, n_cpus = 1, wall_time = gpa$parallel$finalize_time,
    mem_total = "16G",
    r_code = "gpa <- finalize_pipeline_configuration(gpa)", r_packages = "fmri.pipeline",
    batch_code = gpa$parallel$compute_environment,
    scheduler_options = gpa$parallel$sched_args
  )

  if (is.null(gpa$finalize_complete) || isFALSE(gpa$finalize_complete)) {
    lg$info("finalize_pipeline_configuration has not been run on this object. We will start with this step.")
    run_finalize <- TRUE
  } else {
    run_finalize <- FALSE
  }

  # Note: one tricky job dependency problem occurs when a job spawns multiple child jobs.
  # In this case, downstream jobs cannot know the job ids to wait for because these have not been generated yet.
  # This occurs in the case of L1 job submissions, where run_feat_sepjobs submit many jobs that need to complete before
  # L2 jobs should start. The new approach is to keep the parent job active until all children complete. This is
  # implemented by wait_for_children in R_batch_job. Here, we run the l1 model setup, the launch all feat runs in batch
  # jobs and wait for all of these to complete before the l1_setup_batch (parent) job completes.

  l1_setup_batch <- NULL
  l1_execute_batch <- NULL
  if (!is.null(model_list$l1_model_names)) {
    # batch job for setting up l1 models -- calls setup_l1_models to create relevant FSFs
      l1_setup_batch <- f_batch$copy(
        job_name = "setup_l1", n_cpus = gpa$parallel$l1_setup_cores,
        wall_time = gpa$parallel$l1_setup_time,
        r_code = sprintf(
          "gpa <- setup_l1_models(gpa, l1_model_names=%s)", paste(deparse(model_list$l1_model_name), collapse = "")
        )
      )
      l1_setup_batch$mem_total <- "24G"

      l1_setup_batch$depends_on_parents <- "finalize_configuration"

      # batch job for executing l1 jobs (and waiting) after setup
      l1_execute_batch <- f_batch$copy(
        job_name = "run_l1", n_cpus = 1,
        wall_time = gpa$parallel$fsl$l1_feat_alljobs_time,
        r_code = "child_job_ids <- run_feat_sepjobs(gpa, level = 1L)" # execute l1 jobs
      )

      l1_execute_batch$depends_on_parents <- "setup_l1"
      l1_execute_batch$wait_for_children <- TRUE # need to wait for l1 feat jobs to complete before moving to l2/l3
  }


  # todo
  # gpa <- verify_lv1_runs(gpa)

  l2_batch <- NULL

  # only run level 2 if this is a multi-run dataset and user requests l2 model estimation
  if (isTRUE(gpa$multi_run) && !is.null(model_list$l2_model_names)) {
    # setup of l2 models (should follow l1)
    l2_batch <- f_batch$copy(
      job_name = "setup_run_l2", n_cpus = gpa$parallel$l2_setup_cores,
      wall_time = gpa$parallel$l2_setup_run_time,
      r_code = c(
        "gpa <- setup_l2_models(gpa)",
        "child_job_ids <- run_feat_sepjobs(gpa, level = 2L)"
      )
    )

    l2_batch$depends_on_parents <- "run_l1"
    l2_batch$wait_for_children <- TRUE # need to wait for l2 feat jobs to complete before moving to l3
  }

  if (!is.null(model_list$l3_model_names)) {
    l3_batch <- f_batch$copy(
      job_name = "setup_run_l3", n_cpus = gpa$parallel$l2_setup_cores,
      wall_time = gpa$parallel$l3_setup_run_time,
      r_code = c(
        #"gpa <- setup_l3_models(gpa)", # doesn't pass forward subset
        sprintf("gpa <- setup_l3_models(gpa, l3_model_names=%s)", paste(deparse(model_list$l3_model_name), collapse = "")),
        "child_job_ids <- run_feat_sepjobs(gpa, level = 3L)"
      )
    )

    l3_batch$depends_on_parents <- ifelse(isTRUE(gpa$multi_run), "setup_run_l2", "run_l1")
    l3_batch$wait_for_children <- TRUE # need to wait for l3 feat jobs to complete before moving to cleanup
  }

  # cleanup step: refresh l3 feat status and copy gpa back to main directory
  cleanup_batch <- f_batch$copy(
    job_name = "cleanup_fsl", n_cpus = 1,
    wall_time = "30:00", # 30 minutes should be plenty
    r_code = c(
      "gpa <- cleanup_glm_pipeline(gpa)"
    )
  )

  cleanup_batch$depends_on_parents <- "setup_run_l3"

  if (isTRUE(run_finalize)) {
    glm_batch <- R_batch_sequence$new(f_batch, l1_setup_batch, l1_execute_batch, l2_batch, l3_batch, cleanup_batch)
  } else {
    glm_batch <- R_batch_sequence$new(l1_setup_batch, l1_execute_batch, l2_batch, l3_batch, cleanup_batch)
  }
  glm_batch$submit()
}

#' helper function to guide user through process of choosing which models to run in GLM pipeline
#'
#' @param l1_model_names a character vector of level 1 model names
#' @param l2_model_names a character vector of level 2 model names
#' @param l3_model_names a character vector of level 3 model names
#' 
#' @return a named list containing all models that were selected along with additional information
#'   about whether to rerun existing models
#' @keywords internal
choose_glm_set <- function(gpa, l1_model_names=NULL, l2_model_names=NULL, l3_model_names=NULL, lg=NULL) {
  checkmate::assert_class(lg, "Logger")

  if (is.null(l1_model_names)) {
    lg$debug("l1_model_names was NULL. Defaulting to running all l1 models")
    l1_model_names <- "all"
  }

  if (is.null(l2_model_names)) {
    lg$debug("l2_model_names was NULL. Defaulting to running all l2 models")
    l2_model_names <- "all"
  }

  if (is.null(l3_model_names)) {
    lg$debug("l3_model_names was NULL. Defaulting to running all l3 models")
    l3_model_names <- "all"
  }

  m_list <- named_list(l1_model_names, l2_model_names, l3_model_names)

  for (nn in names(m_list)) {
    # enforce that all, none, and prompt must be singleton arguments
    if (any(m_list[[nn]] %in% c("all", "none", "prompt")) && length(m_list[[nn]]) > 1L) {
      msg <- sprintf(
        "Argument %s has value 'all', 'none', or 'prompt'. These must be passed alone, not with other model names.", 
        m_list[[nn]]
      )
      lg$error(msg)
      stop(msg)
    }
  }

  m_string <- function(str) { if (is.null(str)) "none" else str }

  models_specified <- FALSE
  while (isFALSE(models_specified)) {
    l1_model_names <- choose_glm_models(gpa, l1_model_names, level = 1)
    if (isTRUE(gpa$multi_run)) l2_model_names <- choose_glm_models(gpa, l2_model_names, level = 2)
    l3_model_names <- choose_glm_models(gpa, l3_model_names, level = 3)

    cat("\nGLM models to run:\n------------------\n\n")
    cat("Level 1: ", paste(m_string(l1_model_names), collapse = ", "), "\n")
    if (isTRUE(gpa$multi_run)) cat("Level 2: ", paste(m_string(l2_model_names), collapse = ", "), "\n")
    cat("Level 3: ", paste(m_string(l3_model_names), collapse = ", "), "\n")
    cat("\n------------------\n\n")

    if (isTRUE(gpa$multi_run)) {
      menu_options <- c(
        "Yes (run)", "No, respecify level 1 models",
        "No, respecify level 2 models", "No, respecify level 3 models", "Cancel"
      )
    } else {
      menu_options <- c(
        "Yes (run)", "No, respecify level 1 models",
        "No, respecify level 3 models", "Cancel"
      )
    }

    cat("\n")
    respecify <- menu(menu_options, title = "Do you want to continue (this submits the models for execution)?")
    if (respecify == 0L) {
      respecify <- "Cancel"
    } else {
      respecify <- menu_options[respecify]
    }

    if (respecify == "Yes (run)") {
      models_specified <- TRUE
    } else if (respecify == "No, respecify level 1 models") {
      l1_model_names <- "prompt"
    } else if (respecify == "No, respecify level 2 models") {
      l2_model_names <- "prompt"
    } else if (respecify == "No, respecify level 3 models") {
      l3_model_names <- "prompt"
    } else if (respecify == "Cancel") {
      return(invisible(NULL))
    }
  }

  return(named_list(l1_model_names, l2_model_names, l3_model_names))

}