#' This function generates the inputs for an FSL level 2 analysis, where multiple runs for a subject are combined using
#' fixed effects estimation.
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing analysis speecification
#' @param l2_model_names a subset of L2 models to be setup by this function. If not specified,
#'   all models in gpa$l2_models will be included
#' @param l1_model_names a subset of L1 models to be passed to L2 by this function. If not
#'   specified, all models in gpa$l1_models will be included
#' @param backend optional backend filter (e.g., "fsl"). If supplied, only those backends
#'   will be processed.
#'
#' @details
#'   This function will setup FSL level 2 (subject) .fsf files for all combinations of
#'   \code{l2_model_names} and \code{l1_model_names}.
#'
#' @author Michael Hallquist
#' @importFrom checkmate assert_class assert_character assert_data_frame
#' @importFrom lgr get_logger
#' @importFrom iterators iter
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach registerDoSEQ
#' @export
setup_l2_models <- function(gpa, l1_model_names=NULL, l2_model_names=NULL, backend = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_class(gpa$l2_models, "hi_model_set")
  checkmate::assert_class(gpa$l1_models, "l1_model_set")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_subset(l1_model_names, names(gpa$l1_models$models))
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_subset(l2_model_names, names(gpa$l2_models$models))
  checkmate::assert_character(backend, null.ok = TRUE)

  glm_backends <- get_glm_backends(gpa)
  backend_names <- names(glm_backends)
  if (!is.null(backend)) {
    backend <- tolower(backend)
    checkmate::assert_subset(backend, c("fsl", "spm", "afni"))
    missing_backends <- setdiff(backend, backend_names)
    if (length(missing_backends) > 0L) {
      warning(
        sprintf("Backends not available for setup_l2_models: %s", paste(missing_backends, collapse = ", "))
      )
    }
    backend_names <- intersect(backend_names, backend)
  }

  # Filter to backends that support L2 (currently only FSL)
  # SPM concatenates runs at L1, AFNI is not implemented
 l2_supported_backends <- c("fsl")
  backends_without_l2 <- setdiff(backend_names, l2_supported_backends)
  if (length(backends_without_l2) > 0L) {
    # Log but don't warn - this is expected behavior for SPM
    backend_names <- intersect(backend_names, l2_supported_backends)
  }

  use_fsl <- "fsl" %in% backend_names

  # if no l2 model subset is requested, output all models
  if (is.null(l2_model_names)) l2_model_names <- names(gpa$l2_models$models)

  # if no l1 model subset is requested, output all models
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)

  lg <- lgr::get_logger("glm_pipeline/l2_setup")
  lg$set_threshold(gpa$lgr_threshold)

  add_log_suffix <- function(path, suffix) {
    if (is.null(path) || is.null(suffix) || !nzchar(suffix)) return(path)
    ext <- tools::file_ext(path)
    if (nzchar(ext)) {
      base <- sub(paste0("\\.", ext, "$"), "", path)
      return(paste0(base, "_", suffix, ".", ext))
    }
    paste0(path, "_", suffix)
  }
  log_suffix <- if (!is.null(backend)) paste(backend, collapse = "_") else NULL
  setup_l2_log_txt <- add_log_suffix(gpa$output_locations$setup_l2_log_txt, log_suffix)
  setup_l2_log_json <- add_log_suffix(gpa$output_locations$setup_l2_log_json, log_suffix)

  if (isTRUE(gpa$log_txt) && !"setup_l2_log_txt" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderFile$new(setup_l2_log_txt), name = "setup_l2_log_txt")
  }

  if (isTRUE(gpa$log_json) && !"setup_l2_log_json" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderJson$new(setup_l2_log_json), name = "setup_l2_log_json")
  }

  lg$debug("In setup_l2_models, setting up the following L2 models:")
  lg$debug("L2 model: %s", l2_model_names)
  lg$debug("In setup_l2_models, passing the following L1 models to L2:")
  lg$debug("L1 model: %s", l1_model_names)

  if (length(backend_names) == 0L) {
    lg$info("No L2-capable backends requested. Only FSL supports L2 (SPM concatenates runs at L1).")
    gpa$l2_setup_status <- list(
      success = TRUE,
      reason = "no_l2_backends",
      n_models = 0L,
      timestamp = Sys.time()
    )
    return(gpa)
  }

  # setup parallel worker pool, if requested
  if (!is.null(gpa$parallel$l2_setup_cores) && gpa$parallel$l2_setup_cores > 1L) {
    lg$info("Initializing l2 setup cluster with %d cores", gpa$parallel$l2_setup_cores)
    cl <- parallel::makeCluster(gpa$parallel$l2_setup_cores)
    doParallel::registerDoParallel(cl)
    on.exit(try(parallel::stopCluster(cl))) # cleanup pool upon exit of this function
  } else {
    lg$info("Initializing l2 setup with serial execution")
    foreach::registerDoSEQ() # formally register a sequential 'pool' so that dopar is okay
  }

  excluded_runs <- gpa$run_data %>%
    dplyr::select(id, session, run_number, exclude_run, exclude_subject) %>%
    dplyr::filter(exclude_run == TRUE | exclude_subject == TRUE)

  if (nrow(excluded_runs) > 0L) {
    lg$info("In setup_l2_models, the following runs will be excluded from L2 modeling: ")
    lg$info(
      "  subject: %s, session: %s, run_number: %s",
      excluded_runs$id, excluded_runs$session, excluded_runs$run_number
    )
  }

  # only retain good runs and subjects
  run_data <- gpa$run_data %>%
    dplyr::filter(exclude_run == FALSE & exclude_subject == FALSE)

  # subset basic metadata to merge against a given l1 model to enforce run/subject exclusions
  good_runs <- run_data %>%
    dplyr::select(id, session, run_number, exclude_run, exclude_subject)

  if (nrow(run_data) == 0L) {
    msg <- "In setup_l2_models, no runs survived the exclude_subject and exclude_run step."
    lg$warn(msg)
    warning(msg)
    gpa$l2_setup_status <- list(
      success = FALSE,
      reason = msg,
      n_excluded = nrow(excluded_runs),
      timestamp = Sys.time()
    )
    return(gpa)
  }

  if (is.null(gpa$l1_model_setup) || !inherits(gpa$l1_model_setup, "l1_setup")) {
    lg$error("No l1_model_setup found in the glm pipeline object.")
    lg$error("You must run setup_l1_models before running setup_l2_models.")
    stop("No l1_model_setup found in the glm pipeline object.",
    "You must run setup_l1_models before running setup_l2_models.")
  }

  # respecify L2 models for each subject based on available runs
  for (mname in l2_model_names) {
    lg$info("Recalculating per-subject L2 models based on available runs for model: %s", mname)
    gpa$l2_models$models[[mname]] <- respecify_l2_models_by_subject(gpa$l2_models$models[[mname]], run_data)
  }

  # refresh l1 model status in $l1_model_setup
  gpa <- refresh_glm_status(gpa, level=1L, lg=lg)

  # loop over requested backends and setup all requested combinations of L1/L2 models
  backend_results <- list()
  for (backend_name in backend_names) {
    helper <- get_l2_backend_helper(backend_name)
    if (is.null(helper)) {
      lg$warn("No L2 helper registered for backend '%s'. Skipping.", backend_name)
      next
    }

    res <- helper(
      gpa = gpa,
      backend = glm_backends[[backend_name]],
      lg = lg,
      l1_model_names = l1_model_names,
      l2_model_names = l2_model_names,
      good_runs = good_runs
    )
    backend_results[[backend_name]] <- res
  }

  all_subj_l2_combined <- setNames(lapply(backend_results, `[[`, "data"), names(backend_results))
  class(all_subj_l2_combined) <- c("l2_setup", "list")

  # combine with any existing l2 setup to avoid overwriting other backends
  existing_l2 <- if (!is.null(gpa$l2_model_setup) && inherits(gpa$l2_model_setup, "l2_setup")) gpa$l2_model_setup else NULL
  if (!is.null(existing_l2)) {
    backend_union <- union(names(existing_l2), names(backend_results))
    for (backend_name in backend_union) {
      new_df <- all_subj_l2_combined[[backend_name]]
      existing_df <- existing_l2[[backend_name]]
      res <- backend_results[[backend_name]]
      id_cols <- if (!is.null(res)) res$id_cols else NULL

      if (!is.null(existing_df) && !is.null(new_df) && is.data.frame(new_df) &&
        nrow(new_df) > 0L && !is.null(id_cols)) {
        all_subj_l2_combined[[backend_name]] <- update_df(
          current = existing_df, new = new_df, id_cols = id_cols
        )
      } else if (is.null(new_df) || (is.data.frame(new_df) && nrow(new_df) == 0L)) {
        all_subj_l2_combined[[backend_name]] <- existing_df
      }
    }
  }

  # append l2 setup to gpa
  gpa$l2_model_setup <- all_subj_l2_combined

  # refresh l2 model status in $l2_model_setup
  gpa <- refresh_glm_status(gpa, level = 2L, lg = lg)

  # record successful completion status
  n_models <- 0L
  for (backend_name in names(backend_results)) {
    df <- all_subj_l2_combined[[backend_name]]
    if (!is.null(df)) n_models <- n_models + nrow(df)
  }

  gpa$l2_setup_status <- list(
    success = TRUE,
    n_models = n_models,
    backend = backend,
    timestamp = Sys.time()
  )

  return(gpa)
}


get_l2_backend_helper <- function(backend_name) {
  backend_name <- tolower(backend_name)
  switch(
    backend_name,
    fsl = setup_l2_backend_fsl,
    spm = setup_l2_backend_spm,
    afni = setup_l2_backend_afni,
    NULL
  )
}

setup_l2_backend_fsl <- function(gpa, backend, lg, l1_model_names, l2_model_names, good_runs) {
  if (is.null(backend)) {
    return(list(data = NULL, id_cols = NULL))
  }

  model_set <- expand.grid(l1_model = l1_model_names, l2_model = l2_model_names, stringsAsFactors = FALSE)
  all_l2_list <- foreach(
    model_info = iter(model_set, by = "row"), .inorder = FALSE,
    .packages = c("fmri.pipeline", "dplyr", "data.table"),
    .export = c("backend")
  ) %dopar% {
    model_info <- model_info # avoid complaints about visible global binding in R CMD check
    this_l1_model <- model_info$l1_model
    this_l2_model <- model_info$l2_model

    l2_file_setup <- list(fsl = list(), spm = list(), afni = list())

    # get list of runs to examine/include
    to_run <- gpa$l1_model_setup$fsl %>%
      dplyr::filter(l1_model == !!this_l1_model) %>%
      dplyr::select(id, session, run_number, l1_model, feat_fsf, feat_dir)

    # handle run and subject exclusions by joining against good runs
    to_run <- dplyr::inner_join(good_runs, to_run, by = c("id", "session", "run_number"))
    data.table::setDT(to_run) # convert to data.table for split

    by_subj_session <- split(to_run, by = c("id", "session"))

    # setup Feat L2 files for each id and session
    for (l1_df in by_subj_session) {
      subj_id <- l1_df$id[1L]
      subj_session <- l1_df$session[1L]
      feat_l2_df <- tryCatch({
        backend$l2_setup(
          l1_df = l1_df,
          l2_model = this_l2_model, gpa = gpa
        )},
        error = function(e) {
          lg$error(
            "Problem with fsl_l2_model. L1 Model: %s, L2 Model: %s, Subject: %s, Session: %s",
            this_l1_model, this_l2_model, subj_id, subj_session
          )
          lg$error("Error message: %s", as.character(e))
          return(NULL)
        }
      )

      if (!is.null(feat_l2_df)) {
        # add to tracking data.frame
        l2_file_setup$fsl <- rbind(l2_file_setup$fsl, feat_l2_df)
      }
    }

    return(l2_file_setup)
  }

  fsl_df <- data.table::rbindlist(lapply(all_l2_list, "[[", "fsl"), fill = TRUE)
  list(
    data = fsl_df,
    id_cols = c("id", "session", "l1_model", "l2_model")
  )
}

setup_l2_backend_spm <- function(gpa, backend, lg, l1_model_names, l2_model_names, good_runs) {
  lg$warn("spm not supported in setup_l2_models")
  list(data = NULL, id_cols = NULL)
}

setup_l2_backend_afni <- function(gpa, backend, lg, l1_model_names, l2_model_names, good_runs) {
  lg$warn("afni not supported in setup_l2_models")
  list(data = NULL, id_cols = NULL)
}
