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
    batch_code = get_compute_environment(gpa),
    scheduler_options = gpa$parallel$sched_args,
    sqlite_db = gpa$output_locations$sqlite_db
  )

  if (is.null(gpa$finalize_complete) || isFALSE(gpa$finalize_complete)) {
    lg$info("finalize_pipeline_configuration has not been run on this object. We will start with this step.")
    run_finalize <- TRUE
  } else {
    run_finalize <- FALSE
  }

  # initialize batch sequence elements as NULL so that they are ignored in R_batch_sequence if not used
  l1_setup_batch <- split_cache_batch <- l2_batch <- NULL
  sync_l2_backend_caches_batch <- NULL
  l1_execute_batches <- list()
  l3_execute_batches <- list()
  backend_specs <- gpa$glm_backend_specs
  if (is.null(backend_specs)) backend_specs <- default_glm_backend_specs()
  requested_backend_names <- unique(gpa$glm_software)
  resolved_backends <- resolve_glm_backends(backend_specs)
  glm_backends <- resolved_backends[intersect(requested_backend_names, names(resolved_backends))]
  backend_names <- names(glm_backends)
  l3_dependency_backends <- if (!is.null(model_list$l3_model_names) && isTRUE(gpa$multi_run)) {
    get_l3_dependency_backends(glm_backends)
  } else {
    character(0)
  }
  execution_backend_names <- unique(c(backend_names, l3_dependency_backends))
  missing_execution_backends <- setdiff(execution_backend_names, names(resolved_backends))
  if (length(missing_execution_backends) > 0L) {
    stop(
      sprintf(
        "Missing backend specs for required execution backends: %s",
        paste(missing_execution_backends, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  execution_glm_backends <- resolved_backends[execution_backend_names]
  if (!identical(sort(unique(gpa$glm_software)), sort(execution_backend_names))) {
    gpa$glm_software <- execution_backend_names
  }
  use_fsl <- "fsl" %in% execution_backend_names
  l3_requires_l2 <- isTRUE(gpa$multi_run) && !is.null(model_list$l3_model_names) &&
    any_l3_backend_requires_l2(glm_backends)
  backend_cache_map <- setNames(
    file.path(batch_directory, paste0("run_pipeline_cache_", execution_backend_names, ".RData")),
    execution_backend_names
  )

  # Note: one tricky job dependency problem occurs when a job spawns multiple child jobs.
  # In this case, downstream jobs cannot know the job ids to wait for because these have not been generated yet.
  # This occurs in the case of L1 job submissions, where run_feat_sepjobs submit many jobs that need to complete before
  # L2 jobs should start. The new approach is to keep the parent job active until all children complete. This is
  # implemented by wait_for_children in R_batch_job. Here, we run the l1 model setup, the launch all feat runs in batch
  # jobs and wait for all of these to complete before the l1_setup_batch (parent) job completes.

  if (!is.null(model_list$l1_model_names)) {
    # batch job for setting up l1 models -- calls setup_l1_models to create relevant FSFs
    l1_setup_batch <- f_batch$copy(
      job_name = "setup_l1", n_cpus = gpa$parallel$l1_setup_cores,
      wall_time = gpa$parallel$l1_setup_time, mem_total = gpa$parallel$l1_setup_memgb,
      r_code = sprintf(
        "gpa <- setup_l1_models(gpa, l1_model_names=%s)", paste(deparse(model_list$l1_model_name), collapse = "")
      )
    )

    l1_setup_batch$depends_on_parents <- "finalize_configuration"

    if (length(execution_backend_names) > 0L) {
      split_cache_batch <- f_batch$copy(
        job_name = "split_backend_caches", n_cpus = 1,
        wall_time = "0:10:00",
        r_code = c(
          sprintf("shared_cache <- '%s'", gpa_cache),
          "if (!file.exists(shared_cache)) {",
          "  stop('Shared GPA cache not found: ', shared_cache, '. The setup_l1 job may have failed.')",
          "}",
          sprintf(
            "backend_cache_paths <- c(%s)",
            paste(shQuote(unname(backend_cache_map)), collapse = ", ")
          ),
          "for (path in backend_cache_paths) {",
          "  success <- file.copy(shared_cache, path, overwrite = TRUE)",
          "  if (!isTRUE(success)) {",
          "    stop('Failed to copy shared cache to backend cache: ', path)",
          "  }",
          "}"
        )
      )

      split_cache_batch$depends_on_parents <- "setup_l1"
    }

    # batch jobs for executing l1 jobs (and waiting) after setup, per backend
    backend_l1_exec_time <- function(backend_name) {
      if (backend_name == "fsl") {
        return(gpa$parallel$fsl$l1_feat_alljobs_time)
      }
      if (backend_name == "spm") {
        return(gpa$parallel$spm$l1_spm_alljobs_time)
      }
      stop(sprintf("No L1 execution time configured for backend '%s'.", backend_name))
    }

    for (backend_name in execution_backend_names) {
      spec <- backend_specs[[backend_name]]
      run_fn_name <- if (is.null(spec)) NULL else spec$l1_run
      if (is.null(run_fn_name) || identical(run_fn_name, "__not_implemented__")) next
      if (!is.character(run_fn_name)) {
        lg$warn("Skipping backend '%s' L1 runner because l1_run is not a character function name.", backend_name)
        next
      }
      if (!is.null(execution_glm_backends[[backend_name]]$l1_run) &&
        isTRUE(attr(execution_glm_backends[[backend_name]]$l1_run, "glm_backend_not_implemented"))) {
        next
      }

      l1_exec <- f_batch$copy(
        job_name = paste0("run_l1_", backend_name), n_cpus = 1,
        wall_time = backend_l1_exec_time(backend_name),
        r_code = c(
          "child_job_ids <- c()",
          sprintf("child_job_ids <- c(child_job_ids, %s(gpa, level = 1L))", run_fn_name)
        ),
        post_children_r_code = c(
          sprintf("gpa <- refresh_glm_status(gpa, level = 1L, glm_software = '%s')", backend_name)
        ) # refresh l1 status after all jobs complete
      )

      l1_exec$depends_on_parents <- if (!is.null(split_cache_batch)) "split_backend_caches" else "setup_l1"
      l1_exec$wait_for_children <- TRUE # need to wait for l1 jobs to complete before moving to l2/l3
      if (!is.null(backend_cache_map[[backend_name]])) {
        l1_exec$input_rdata_file <- backend_cache_map[[backend_name]]
        l1_exec$output_rdata_file <- backend_cache_map[[backend_name]]
      }
      l1_execute_batches[[backend_name]] <- l1_exec
    }

    run_l1 <- TRUE
  } else {
    run_l1 <- FALSE
  }

  # If L1 is skipped but later levels will run, split the shared cache for each backend.
  if (isFALSE(run_l1) && is.null(split_cache_batch) && length(execution_backend_names) > 0L &&
      (!is.null(model_list$l2_model_names) || !is.null(model_list$l3_model_names))) {
    split_cache_batch <- f_batch$copy(
      job_name = "split_backend_caches", n_cpus = 1,
      wall_time = "0:10:00",
      r_code = c(
        sprintf("shared_cache <- '%s'", gpa_cache),
        "if (!file.exists(shared_cache)) {",
        "  stop('Shared GPA cache not found: ', shared_cache, '. The finalize job may have failed.')",
        "}",
        sprintf(
          "backend_cache_paths <- c(%s)",
          paste(shQuote(unname(backend_cache_map)), collapse = ", ")
        ),
        "for (path in backend_cache_paths) {",
        "  success <- file.copy(shared_cache, path, overwrite = TRUE)",
        "  if (!isTRUE(success)) {",
        "    stop('Failed to copy shared cache to backend cache: ', path)",
        "  }",
        "}"
      )
    )
    split_cache_batch$depends_on_parents <- if (isTRUE(run_finalize)) "finalize_configuration" else NULL
  }


  # todo
  # gpa <- verify_lv1_runs(gpa)

  # L2 is a prerequisite for any multi-run L3 backend that declares a dependency on it.
  # AFNI 3dLMEr consumes FSL-produced L2 COPEs even when AFNI is the only requested L3 backend.
  can_run_l2 <- isTRUE(gpa$multi_run) && isTRUE(use_fsl)
  if (!is.null(model_list$l2_model_names) && isFALSE(can_run_l2)) {
    reasons <- c()
    if (!isTRUE(gpa$multi_run)) reasons <- c(reasons, "dataset is not multi-run")
    if (!isTRUE(use_fsl)) reasons <- c(reasons, "FSL backend is not enabled")
    spm_note <- if ("spm" %in% backend_names) {
      "SPM does not run standalone L2 jobs; its L2 model spec is projected during SPM L1 setup."
    } else {
      NULL
    }
    lg$warn(
      "L2 models were requested but cannot be run as standalone jobs: %s%s",
      if (length(reasons) > 0L) paste(reasons, collapse = "; ") else "unknown reason",
      if (!is.null(spm_note)) paste0(". ", spm_note) else "."
    )
  }
  if (isTRUE(can_run_l2) && (!is.null(model_list$l2_model_names) || isTRUE(l3_requires_l2))) {
    # setup of l2 models (should follow l1)
    l2_batch <- f_batch$copy(
      job_name = "setup_run_l2", n_cpus = gpa$parallel$l2_setup_cores,
      wall_time = gpa$parallel$l2_setup_run_time,
      r_code = c(
        "gpa <- setup_l2_models(gpa, backend = 'fsl')",
        "child_job_ids <- run_feat_sepjobs(gpa, level = 2L)"
      )
    )

    if (isTRUE(run_l1)) {
      l2_batch$depends_on_parents <- "run_l1_fsl"
    } else if (!is.null(split_cache_batch)) {
      l2_batch$depends_on_parents <- "split_backend_caches"
    } else if (isTRUE(run_finalize)) {
      l2_batch$depends_on_parents <- "finalize_configuration"
    } else {
      l2_batch$depends_on_parents <- NULL
    }
    l2_batch$wait_for_children <- TRUE # need to wait for l2 feat jobs to complete before moving to l3
    if (!is.null(backend_cache_map[["fsl"]])) {
      l2_batch$input_rdata_file <- backend_cache_map[["fsl"]]
      l2_batch$output_rdata_file <- backend_cache_map[["fsl"]]
    }
    
    run_l2 <- TRUE
  } else {
    run_l2 <- FALSE
  }

  dependent_l3_backends <- if (!is.null(model_list$l3_model_names)) {
    backend_names[vapply(backend_names, function(backend_name) {
      source_backend <- backend_l3_l2_source_backend(glm_backends[[backend_name]])
      isTRUE(gpa$multi_run) &&
        backend_l3_requires_l2(glm_backends[[backend_name]]) &&
        !is.null(source_backend) &&
        !identical(source_backend, backend_name)
    }, logical(1))]
  } else {
    character(0)
  }

  if (isTRUE(run_l2) && length(dependent_l3_backends) > 0L) {
    sync_lines <- c()
    for (backend_name in dependent_l3_backends) {
      source_backend <- backend_l3_l2_source_backend(glm_backends[[backend_name]])
      sync_lines <- c(
        sync_lines,
        sprintf("src <- %s", shQuote(backend_cache_map[[source_backend]])),
        sprintf("dest <- %s", shQuote(backend_cache_map[[backend_name]])),
        "if (!file.exists(src)) stop('Required backend cache not found: ', src)",
        "success <- file.copy(src, dest, overwrite = TRUE)",
        "if (!isTRUE(success)) stop('Failed to sync backend cache from ', src, ' to ', dest)"
      )
    }
    sync_l2_backend_caches_batch <- f_batch$copy(
      job_name = "sync_l2_backend_caches",
      n_cpus = 1,
      wall_time = "0:10:00",
      r_code = sync_lines
    )
    sync_l2_backend_caches_batch$depends_on_parents <- "setup_run_l2"
  }

  
  if (!is.null(model_list$l3_model_names)) {
    for (backend_name in backend_names) {
      spec <- backend_specs[[backend_name]]
      run_fn_name <- if (is.null(spec)) NULL else spec$l3_run
      if (is.null(run_fn_name) || identical(run_fn_name, "__not_implemented__")) next
      if (!is.character(run_fn_name)) {
        lg$warn("Skipping backend '%s' L3 runner because l3_run is not a character function name.", backend_name)
        next
      }
      if (!is.null(glm_backends[[backend_name]]$l3_run) &&
        isTRUE(attr(glm_backends[[backend_name]]$l3_run, "glm_backend_not_implemented"))) {
        next
      }

      l1_parent <- if (!is.null(l1_execute_batches[[backend_name]])) {
        paste0("run_l1_", backend_name)
      } else if (!is.null(split_cache_batch)) {
        "split_backend_caches"
      } else if (!is.null(l1_setup_batch)) {
        "setup_l1"
      } else if (isTRUE(run_finalize)) {
        "finalize_configuration"
      } else {
        NULL
      }

      l3_exec <- f_batch$copy(
        job_name = paste0("setup_run_l3_", backend_name), n_cpus = gpa$parallel$l2_setup_cores,
        wall_time = gpa$parallel$l3_setup_run_time,
        r_code = c(
          sprintf(
            "gpa <- setup_l3_models(gpa, l3_model_names=%s, backend = '%s')",
            paste(deparse(model_list$l3_model_name), collapse = ""), backend_name
          ),
          "child_job_ids <- c()",
          sprintf("child_job_ids <- c(child_job_ids, %s(gpa, level = 3L))", run_fn_name)
        ),
        post_children_r_code = c(
          sprintf("gpa <- refresh_glm_status(gpa, level = 3L, glm_software = '%s')", backend_name)
        )
      )

      if (isTRUE(gpa$multi_run) && backend_l3_requires_l2(glm_backends[[backend_name]]) && isTRUE(run_l2)) {
        source_backend <- backend_l3_l2_source_backend(glm_backends[[backend_name]])
        if (!is.null(sync_l2_backend_caches_batch) && !identical(source_backend, backend_name)) {
          l3_exec$depends_on_parents <- "sync_l2_backend_caches"
        } else {
          l3_exec$depends_on_parents <- "setup_run_l2"
        }
      } else {
        l3_exec$depends_on_parents <- l1_parent
      }
      l3_exec$wait_for_children <- TRUE
      if (!is.null(backend_cache_map[[backend_name]])) {
        l3_exec$input_rdata_file <- backend_cache_map[[backend_name]]
        l3_exec$output_rdata_file <- backend_cache_map[[backend_name]]
      }
      l3_execute_batches[[backend_name]] <- l3_exec
    }

    run_l3 <- length(l3_execute_batches) > 0L
  } else {
    run_l3 <- FALSE
  }

  # cleanup step: refresh l3 feat status and copy gpa back to main directory
  backend_cache_vec <- NULL
  if (length(backend_cache_map) > 0L) {
    backend_cache_vec <- paste(
      sprintf("%s=%s", names(backend_cache_map), shQuote(unname(backend_cache_map))),
      collapse = ", "
    )
  }

  cleanup_code <- c()
  if (!is.null(backend_cache_vec) && nzchar(backend_cache_vec)) {
    cleanup_code <- c(cleanup_code, sprintf("backend_cache_paths <- c(%s)", backend_cache_vec))
    cleanup_code <- c(cleanup_code, "gpa <- cleanup_glm_pipeline(gpa, backend_cache_paths = backend_cache_paths)")
  } else {
    cleanup_code <- c(cleanup_code, "gpa <- cleanup_glm_pipeline(gpa)")
  }

  cleanup_batch <- f_batch$copy(
    job_name = "cleanup_glm", n_cpus = 1,
    wall_time = "30:00", # 30 minutes should be plenty
    r_code = cleanup_code
  )

  if (length(l3_execute_batches) > 0L) {
    cleanup_batch$depends_on_parents <- vapply(l3_execute_batches, function(x) x$job_name, character(1))
  } else if (!is.null(l2_batch)) {
    cleanup_batch$depends_on_parents <- "setup_run_l2"
  } else if (length(l1_execute_batches) > 0L) {
    cleanup_batch$depends_on_parents <- vapply(l1_execute_batches, function(x) x$job_name, character(1))
  } else {
    cleanup_batch$depends_on_parents <- if (isTRUE(run_finalize)) "finalize_configuration" else NULL
  }
  
  run_cleanup <- TRUE # always TRUE for now

  if (isTRUE(run_finalize)) {
    glm_batch <- R_batch_sequence$new(
      joblist = c(
        list(f_batch, l1_setup_batch, split_cache_batch),
        l1_execute_batches,
        list(l2_batch),
        list(sync_l2_backend_caches_batch),
        l3_execute_batches,
        list(cleanup_batch)
      ),
      sequence_id = batch_id
    )
  } else {
    glm_batch <- R_batch_sequence$new(
      joblist = c(
        list(l1_setup_batch, split_cache_batch),
        l1_execute_batches,
        list(l2_batch),
        list(sync_l2_backend_caches_batch),
        l3_execute_batches,
        list(cleanup_batch)
      ),
      sequence_id = batch_id
    )
  }
  glm_batch$submit()
  
  l1_exec_flags <- if (length(l1_execute_batches) > 0L) {
    setNames(lapply(l1_execute_batches, function(x) !is.null(x)), paste0("l1_", names(l1_execute_batches)))
  } else {
    list()
  }
  l3_exec_flags <- if (length(l3_execute_batches) > 0L) {
    setNames(lapply(l3_execute_batches, function(x) !is.null(x)), paste0("l3_", names(l3_execute_batches)))
  } else {
    list()
  }

  job_sequence <- c(
    list(
      finalize = run_finalize,
      l1 = run_l1,
      split_backend_caches = !is.null(split_cache_batch),
      l2 = run_l2,
      l3 = run_l3,
      cleanup = run_cleanup
    ),
    l1_exec_flags,
    l3_exec_flags
  ) # logicals of which jobs were submitted
  update_project_config(gpa = gpa, job_sequence = job_sequence, 
                        sequence_id = batch_id, batch_directory) # update config file with this run's details
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
  is_prompt <- function(x) {
    !is.null(x) && length(x) > 0L && identical(x[1L], "prompt")
  }

  prompt_selection_requested <- is_prompt(l1_model_names) ||
    is_prompt(l3_model_names) ||
    (isTRUE(gpa$multi_run) && is_prompt(l2_model_names))

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

    if (!interactive()) {
      lg$info("Non-interactive session detected. Proceeding with selected models without prompt.")
      models_specified <- TRUE
      next
    }

    if (isFALSE(prompt_selection_requested)) {
      lg$info("Model names were provided as input arguments. Proceeding without confirmation prompt.")
      models_specified <- TRUE
      next
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
