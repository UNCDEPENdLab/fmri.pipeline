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
  checkmate::assert_string(l1_model_names, null.ok = TRUE)
  checkmate::assert_string(l2_model_names, null.ok = TRUE)
  checkmate::assert_string(l3_model_names, null.ok = TRUE)

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")

  choose_models <- function(gpa, model_names, level) {
    checkmate::assert_integerish(level, min=1, max=3)
    all_m_names <- names(gpa[[paste0("l", level, "_models")]]$models)
    checkmate::assert_subset(model_names, c("prompt", all_m_names))

    if (is.null(model_names)) {
      # null indicates that we should run all models at this level
      chosen_models <- all_m_names
    } else if (model_names[1] == "prompt") {
      chosen_models <- select.list(all_m_names,
        multiple = TRUE,
        title = paste("Choose all level", level, "models to include in this pipeline run: ")
      )
      if (identical(chosen_models, character(0))) {
        lg$info(paste("No level", level, "models were selected."))
        chosen_models <- NULL
      }
    } else {
      chosen_models <- model_names # user-specified set
    }
    return(chosen_models)
  }

  models_specified <- FALSE
  while (isFALSE(models_specified)) {
    l1_model_names <- choose_models(gpa, l1_model_names, level = 1)
    if (isTRUE(gpa$multi_run)) l2_model_names <- choose_models(gpa, l2_model_names, level = 2)
    l3_model_names <- choose_models(gpa, l3_model_names, level = 3)

    cat("\nGLM models to run:\n------------------\n\n")
    cat("Level 1: ", paste(l1_model_names, collapse = ", "), "\n")
    if (isTRUE(gpa$multi_run)) cat("Level 2: ", paste(l2_model_names, collapse = ", "), "\n")
    cat("Level 3: ", paste(l3_model_names, collapse = ", "), "\n")

    respecify <- menu(c("Yes (run)", "No (re-choose models)"),
      title = "Do you want to continue (this submits the models for execution)?"
    )
    if (respecify == 1L) models_specified <- TRUE
  }

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")
  batch_id <- uuid::UUIDgenerate()
  if (is.null(gpa$batch_run)) gpa$batch_run <- list()

  # setup a unique directory for this GLM batch submission
  # this will contain all scripts for this run and contain a snapshot of the gpa object
  # that is amended as each batch job runs and completes.
  batch_directory <- file.path(gpa$output_locations$scheduler_scripts, paste0("batch_", batch_id))
  if (!dir.exists(batch_directory)) dir.create(batch_directory, recursive = TRUE)
  gpa_cache <- file.path(batch_directory, "run_pipeline_cache.RData")
  save(gpa, file=gpa_cache)

  if (is.null(gpa$finalize_complete) || isFALSE(gpa$finalize_complete)) {
    lg$info("finalize_pipeline_configuration has not been run on this object. We will start with this step.")
    f_batch <- R_batch_job$new(
      job_name = "finalize_configuration", batch_directory = batch_directory, scheduler = gpa$scheduler,
      input_environment = gpa_cache, output_environment = gpa_cache,
      n_nodes = 1, n_cpus = 1, cpu_time = gpa$parallel$finalize_time,
      r_code = "gpa <- finalize_pipeline_configuration(gpa)", r_packages = "fmri.pipeline",
      batch_code = gpa$parallel$compute_environment
    )

  }

  # Note: one tricky job dependency problem occurs when a job spawns multiple child jobs.
  # In this case, downstream jobs cannot know the job ids to wait for because these have not been generated yet.
  # This occurs in the case of L1 job submissions, where run_feat_sepjobs submit many jobs that need to complete before
  # L2 jobs should start. The new approach is to keep the parent job active until all children complete. This is
  # implemented by wait_for_children in R_batch_job. Here, we run the l1 model setup, the launch all feat runs in batch
  # jobs and wait for all of these to complete before the l1_batch (parent) job completes.

  l1_batch <- f_batch$copy(
    job_name = "setup_run_l1", n_cpus = gpa$parallel$l1_setup_cores,
    cpu_time = gpa$parallel$l1_setup_time,
    r_code = c(
      "gpa <- setup_l1_models(gpa)", # create all FSF files for level one runs
      "child_job_ids <- run_feat_sepjobs(gpa, level = 1L)" # execute l1 jobs
    )
  )

  l1_batch$depends_on_parents <- "finalize_configuration"
  l1_batch$wait_for_children <- TRUE # need to wait for l1 feat jobs to complete before moving to l2/l3

  # todo
  # gpa <- verify_lv1_runs(gpa)

  # only run level 2 if this is a multi-run dataset
  if (isTRUE(gpa$multi_run)) {
    # setup of l2 models (should follow l1)
    l2_batch <- f_batch$copy(
      job_name = "setup_run_l2", n_cpus = gpa$parallel$l2_setup_cores,
      cpu_time = gpa$parallel$l2_setup_time,
      r_code = c(
        "gpa <- setup_l2_models(gpa)",
        "child_job_ids <- run_feat_sepjobs(gpa, level = 2L)"
      )
    )

    l2_batch$depends_on_parents <- "l1_batch"
    l2_batch$wait_for_children <- TRUE # need to wait for l2 feat jobs to complete before moving to l3
  } else {
    l2_batch <- NULL
  }

  l3_batch <- f_batch$copy(
    job_name = "setup_run_l3", n_cpus = gpa$parallel$l2_setup_cores,
    cpu_time = gpa$parallel$l2_setup_time,
    r_code = c(
      "gpa <- setup_l3_models(gpa)",
      "jobs <- run_feat_sepjobs(gpa, level = 3L)"
    )
  )

  l3_batch$depends_on_parents <- ifelse(isTRUE(gpa$multi_run), "l2_batch", "l1_batch")

  glm_batch <- R_batch_sequence$new(f_batch, l1_batch, l2_batch, l3_batch)
  glm_batch$submit()
}
