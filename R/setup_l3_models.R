#' This function generates the inputs for level 3 analyses, where multi-subject data are analyzed
#'   in group analyses.
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing analysis speecification
#' @param l3_model_names a subset of L3 models to be setup by this function. If not specified,
#'   all models in gpa$l2_models will be included
#' @param l2_model_names a subset of L2 models to be setup by this function. If not specified,
#'   all models in gpa$l2_models will be included
#' @param l1_model_names a subset of L1 models to be passed to L2 by this function. If not
#'   specified, all models in gpa$l1_models will be included
#' @param backend optional backend filter (e.g., "spm"). If supplied, only those backends
#'   will be processed.
#'
#' @details
#'   This function will setup FSL level 2 (subject) .fsf files for all combinations of
#'   \code{l2_model_names} and \code{l1_model_names}.
#'
#' FSL 2-level versus 3-level setup
#'
#' 2-level setup (one run per subject)
#'   - Pass L1 .feat folders as input to L3 .fsf setup
#'   - In this approach, the copes in the .fsf pertain to the L1 cope numbers
#'   - Requires one .fsf per L3 model
#'
#' 3-level setup (multiple runs per subject, combined at L2)
#'   - Pass individual cope*.feat folders within subject .gfeat folders
#'   - The folder cope numbers pertain to L1 copes
#'   - The cope*.nii.gz in the cope*.feat subfolders pertain to the L2 contrasts
#'   - Requires one .fsf per L1 cope x L3 model combination
#'   - Example: FSL_L2.gfeat/cope3.feat/stats/cope1.nii.gz
#'      ==> cope3 is the third contrast in the L1 feat model
#'      ==> cope1 is the first contrast in the L2 feat model
#'
#' @author Michael Hallquist
#' @importFrom checkmate assert_class assert_character assert_data_frame
#' @importFrom data.table is.data.table
#' @importFrom lgr get_logger
#' @export
setup_l3_models <- function(gpa, l3_model_names = NULL, l2_model_names = NULL, l1_model_names = NULL, backend = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_character(backend, null.ok = TRUE)
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_data_frame(gpa$subject_data)
  checkmate::assert_class(gpa$l3_models, "hi_model_set")
  checkmate::assert_class(gpa$l1_models, "l1_model_set")

  glm_backends <- get_glm_backends(gpa)
  backend_names <- names(glm_backends)
  if (!is.null(backend)) {
    backend <- tolower(backend)
    checkmate::assert_subset(backend, c("fsl", "spm", "afni"))
    missing_backends <- setdiff(backend, backend_names)
    if (length(missing_backends) > 0L) {
      warning(
        sprintf("Backends not available for setup_l3_models: %s", paste(missing_backends, collapse = ", "))
      )
    }
    backend_names <- intersect(backend_names, backend)
  }

  glm_software_orig <- gpa$glm_software
  if (!is.null(backend)) {
    gpa$glm_software <- backend_names
    on.exit({ gpa$glm_software <- glm_software_orig }, add = TRUE)
  }

  use_fsl <- "fsl" %in% backend_names
  use_spm <- "spm" %in% backend_names
  use_afni <- "afni" %in% backend_names
  requires_l2 <- isTRUE(gpa$multi_run) && isTRUE(use_fsl)

  # Validate model name subsets
  checkmate::assert_subset(l3_model_names, names(gpa$l3_models$models))
  checkmate::assert_subset(l1_model_names, names(gpa$l1_models$models))

  # full 3-level analysis (runs, subject, sample)
  if (isTRUE(gpa$multi_run) && isTRUE(requires_l2)) {
    checkmate::assert_class(gpa$l2_models, "hi_model_set")
    checkmate::assert_subset(l2_model_names, names(gpa$l2_models$models))

    # if no l2 model subset is requested, output all models
    if (is.null(l2_model_names)) l2_model_names <- names(gpa$l2_models$models)
  } else if (isTRUE(gpa$multi_run)) {
    l2_model_names <- NULL
  }

  # if no l3 model subset is requested, output all models
  if (is.null(l3_model_names)) l3_model_names <- names(gpa$l3_models$models)

  # if no l1 model subset is requested, output all models
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)

  lg <- lgr::get_logger("glm_pipeline/l3_setup")
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
  setup_l3_log_txt <- add_log_suffix(gpa$output_locations$setup_l3_log_txt, log_suffix)
  setup_l3_log_json <- add_log_suffix(gpa$output_locations$setup_l3_log_json, log_suffix)

  if (isTRUE(gpa$log_txt) && !"setup_l3_log_txt" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderFile$new(setup_l3_log_txt), name = "setup_l3_log_txt")
  }

  if (isTRUE(gpa$log_json) && !"setup_l3_log_json" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderJson$new(setup_l3_log_json), name = "setup_l3_log_json")
  }

  lg$debug("In setup_l3_models, setting up the following L3 models:")
  lg$debug("L3 model: %s", l3_model_names)
  if (isTRUE(requires_l2)) {
    lg$debug("In setup_l3_models, passing the following L2 models to L3:")
    lg$debug("L2 model: %s", l2_model_names)
  }
  lg$debug("In setup_l3_models, passing the following L1 models to L3:")
  lg$debug("L1 model: %s", l1_model_names)

  if (!isTRUE(use_fsl) && !isTRUE(use_spm) && !isTRUE(use_afni)) {
    lg$info("No supported L3 backends requested in setup_l3_models. Skipping.")
    return(gpa)
  }

  # handle refresh of feat status for lower-level models
  # For multi-run data, L2 inputs must be complete for entry into L3
  # For single-run data, L1 inputs must be complete for entry into L3
  if (isTRUE(requires_l2)) {
    # refresh l2 model status in $l2_model_setup
    gpa <- refresh_glm_status(gpa, level = 2L, lg = lg)
  } else {
    # refresh l1 model status in $l1_model_setup
    gpa <- refresh_glm_status(gpa, level = 1L, lg = lg)
  }

  excluded_runs <- gpa$run_data %>%
    dplyr::select(id, session, run_number, exclude_run, exclude_subject) %>%
    dplyr::filter(exclude_run == TRUE | exclude_subject == TRUE)

  if (nrow(excluded_runs) > 1L) {
    lg$info("In setup_l3_models, the following runs will be excluded from L3 modeling: ")
    lg$info(
      "  subject: %s, session: %s, run_number: %s",
      excluded_runs$id, excluded_runs$session, excluded_runs$run_number
    )
  }

  # make sure l1 models have already been generated
  enforce_glms_complete(gpa, level = 1L, lg)

  if (isTRUE(requires_l2)) {
    lg$info("In setup_l3_models, using a multi-run 3-level setup with runs (l1), subjects (l2), sample (l3)")

    # in multi-run setup, an l2_model_setup must be present
    enforce_glms_complete(gpa, level = 2L, lg)

  } else {
    lg$info("In setup_l3_models, using a single run 2-level setup with subjects (l1), sample (l3)")
    # only retain good runs and subjects
    run_data <- gpa$run_data %>%
      dplyr::filter(exclude_run == FALSE & exclude_subject == FALSE)

    if (nrow(run_data) == 0L) {
      msg <- "In setup_l3_models, no runs survived the exclude_subject and exclude_run step."
      lg$warn(msg)
      warning(msg)
      gpa$l3_setup_status <- list(
        success = FALSE,
        reason = msg,
        n_excluded = nrow(excluded_runs),
        timestamp = Sys.time()
      )
      return(gpa)
    }
  }

  # subjects and sessions to run at l3
  subj_df <- gpa$subject_data %>%
    dplyr::filter(exclude_subject==FALSE) %>%
    select(id, session)

  # loop over requested backends and setup all requested combinations of L1/L2/L3 models
  backend_results <- list()
  metadata <- NULL
  for (backend_name in backend_names) {
    helper <- get_l3_backend_helper(backend_name)
    if (is.null(helper)) {
      lg$warn("No L3 helper registered for backend '%s'. Skipping.", backend_name)
      next
    }

    res <- helper(
      gpa = gpa,
      backend = glm_backends[[backend_name]],
      lg = lg,
      l1_model_names = l1_model_names,
      l2_model_names = l2_model_names,
      l3_model_names = l3_model_names,
      subj_df = subj_df,
      requires_l2 = requires_l2
    )
    backend_results[[backend_name]] <- res
    if (is.null(metadata) && !is.null(res$metadata)) metadata <- res$metadata
  }

  all_subj_l3_combined <- c(
    list(metadata = metadata),
    setNames(lapply(backend_results, `[[`, "data"), names(backend_results))
  )

  class(all_subj_l3_combined) <- c("l3_setup", "list")

  # Combine the generated l3_model_setup with any existing l3_model_setup information
  # This is important so that if a user requests a model subset (using input arguments, and from run_glm_pipeline),
  # we don't clear out information about many other models that may have already completed
  existing_l3 <- if (!is.null(gpa$l3_model_setup) && inherits(gpa$l3_model_setup, "l3_setup")) gpa$l3_model_setup else NULL

  if (!is.null(existing_l3)) {
    if (is.null(all_subj_l3_combined$metadata) && !is.null(existing_l3$metadata)) {
      all_subj_l3_combined$metadata <- existing_l3$metadata
    }

    backend_union <- union(setdiff(names(existing_l3), "metadata"), names(backend_results))
    for (backend_name in backend_union) {
      new_df <- all_subj_l3_combined[[backend_name]]
      existing_df <- existing_l3[[backend_name]]
      res <- backend_results[[backend_name]]
      id_cols <- if (!is.null(res)) res$id_cols else NULL

      if (!is.null(existing_df) && !is.null(new_df) && is.data.frame(new_df) &&
        nrow(new_df) > 0L && !is.null(id_cols)) {
        all_subj_l3_combined[[backend_name]] <- update_df(
          current = existing_df, new = new_df, id_cols = id_cols
        )
      } else if (is.null(new_df) || (is.data.frame(new_df) && nrow(new_df) == 0L)) {
        all_subj_l3_combined[[backend_name]] <- existing_df
      }
    }
  }

  # append l3 setup to gpa
  gpa$l3_model_setup <- all_subj_l3_combined

  # refresh l3 model status in $l3_model_setup
  gpa <- refresh_glm_status(gpa, level = 3L, lg = lg)

  # the expressions for l3 output locations should not generate any duplicate fsfs
  if (isTRUE(use_fsl) && !is.null(gpa$l3_model_setup$fsl)) {
    dupe_fsfs <- duplicated(gpa$l3_model_setup$fsl$feat_fsf)
    if (any(dupe_fsfs, na.rm = TRUE)) {
      lg$warn("There are duplicate fsfs in gpa$l3_model_setup$fsl.")
      lg$warn("This suggests a problem with gpa$output_locations$feat_l3_fsf and gpa$output_locations$feat_l3_directory settings.")
      lg$warn("gpa$output_locations$feat_l3_fsf: %s", gpa$output_locations$feat_l3_fsf)
      lg$warn("gpa$output_locations$feat_l3_directory: %s", gpa$output_locations$feat_l3_directory)
      lg$debug("%s", capture.output(print(gpa$l3_model_setup$fsl[dupe_fsfs, ])))
    }
  }

  n_models <- 0L
  for (backend_name in names(backend_results)) {
    df <- all_subj_l3_combined[[backend_name]]
    if (!is.null(df)) n_models <- n_models + nrow(df)
  }

  # record successful completion status
  gpa$l3_setup_status <- list(
    success = TRUE,
    n_models = n_models,
    backend = backend,
    timestamp = Sys.time()
  )

  return(gpa)
}


get_l3_backend_helper <- function(backend_name) {
  backend_name <- tolower(backend_name)
  switch(
    backend_name,
    fsl = setup_l3_backend_fsl,
    spm = setup_l3_backend_spm,
    afni = setup_l3_backend_afni,
    NULL
  )
}

setup_l3_backend_fsl <- function(gpa, backend, lg, l1_model_names, l2_model_names, l3_model_names, subj_df, requires_l2) {
  if (is.null(backend)) {
    return(list(metadata = NULL, data = NULL, id_cols = NULL))
  }

  if (isTRUE(requires_l2)) {
    model_set <- expand.grid(
      l1_model = l1_model_names, l2_model = l2_model_names,
      l3_model = l3_model_names, stringsAsFactors = FALSE
    )
  } else {
    model_set <- expand.grid(l1_model = l1_model_names, l3_model = l3_model_names, stringsAsFactors = FALSE)
  }

  l3_cope_config <- get_fsl_l3_model_df(gpa, model_set, subj_df)

  l3_cope_input_df <- l3_cope_config %>%
    dplyr::filter(l3_cope_number == 1L) %>%
    dplyr::select(-l3_cope_number, -l3_cope_name)

  to_run <- get_feat_l3_inputs(gpa, l3_cope_input_df, lg)

  if (length(to_run) == 0L) {
    lg$warn("No feat l3 inputs returned from get_feat_l3_inputs. Is it possible no lower-level analyses are complete?")
    return(list(metadata = l3_cope_config, data = NULL,
                id_cols = c("l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model")))
  }

  all_l3_list <- foreach(
    model_info = iter(to_run), .inorder = FALSE,
    .packages = c("fmri.pipeline", "dplyr", "data.table"), .export = c("lg", "gpa", "backend")
  ) %dopar% {
    model_info <- model_info # to avoid complaints about global variable binding in R CMD check

    l3_file_setup <- list(fsl = list(), spm = list(), afni = list())

    if (nrow(model_info) <= 3) {
      lg$warn(
        "Fewer than 4 complete feat input directories for l1 model %s, l2 model %s, l3 model %s",
        model_info$l1_model[1L], model_info$l2_model[1L], model_info$l3_model[1L]
      )
      l3_file_setup$fsl <- NULL
    }

    l3_file_setup$fsl <- tryCatch(backend$l3_setup(model_info, gpa = gpa),
      error = function(e) {
        lg$error(
          "Problem with fsl_l3_model. L1 Model: %s, L2 Model: %s, L3 model %s",
          model_info$l1_model[1L], model_info$l2_model[1L], model_info$l3_model[1L]
        )
        lg$error("Error message: %s", as.character(e))
        return(NULL)
      }
    )

    return(l3_file_setup)
  }

  fsl_df <- dplyr::bind_rows(lapply(all_l3_list, "[[", "fsl"))
  list(
    metadata = l3_cope_config,
    data = fsl_df,
    id_cols = c("l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model")
  )
}

setup_l3_backend_spm <- function(gpa, backend, lg, l1_model_names, l2_model_names, l3_model_names, subj_df, requires_l2) {
  if (is.null(backend)) {
    return(list(metadata = NULL, data = NULL, id_cols = NULL))
  }

  if (isTRUE(gpa$multi_run)) {
    lg$info("In setup_l3_models, using SPM L1->L3 setup with concatenated runs (skipping L2)")
  }

  spm_model_set <- expand.grid(
    l1_model = l1_model_names, l3_model = l3_model_names,
    stringsAsFactors = FALSE
  )

  spm_cope_config <- get_spm_l3_model_df(gpa, spm_model_set, subj_df)

  spm_cope_input_df <- spm_cope_config %>%
    dplyr::filter(l3_cope_number == 1L) %>%
    dplyr::select(-l3_cope_number, -l3_cope_name)

  spm_to_run <- get_spm_l3_inputs(gpa, spm_cope_input_df, lg)

  if (length(spm_to_run) == 0L) {
    lg$warn("No SPM L3 inputs returned. Is it possible no lower-level SPM analyses are complete?")
    return(list(metadata = spm_cope_config, data = NULL,
                id_cols = c("l1_model", "l1_cope_name", "l3_model")))
  }

  all_spm_list <- foreach(
    model_info = iter(spm_to_run), .inorder = FALSE,
    .packages = c("fmri.pipeline", "dplyr", "data.table"), .export = c("lg", "gpa", "backend")
  ) %dopar% {
    model_info <- model_info # to avoid complaints about global variable binding in R CMD check

    l3_file_setup <- list(fsl = list(), spm = list(), afni = list())

    # Validate minimum number of inputs (same check as FSL backend)
    if (nrow(model_info) <= 3) {
      lg$warn(
        "Fewer than 4 complete SPM input directories for l1 model %s, l3 model %s",
        model_info$l1_model[1L], model_info$l3_model[1L]
      )
      l3_file_setup$spm <- NULL
      return(l3_file_setup)
    }

    l3_file_setup$spm <- tryCatch(backend$l3_setup(model_info, gpa = gpa),
      error = function(e) {
        lg$error(
          "Problem with spm_l3_model. L1 Model: %s, L3 model %s",
          model_info$l1_model[1L], model_info$l3_model[1L]
        )
        lg$error("Error message: %s", as.character(e))
        return(NULL)
      }
    )

    return(l3_file_setup)
  }

  spm_df <- dplyr::bind_rows(lapply(all_spm_list, "[[", "spm"))
  list(
    metadata = spm_cope_config,
    data = spm_df,
    id_cols = c("l1_model", "l1_cope_name", "l3_model")
  )
}

setup_l3_backend_afni <- function(gpa, backend, lg, l1_model_names, l2_model_names, l3_model_names, subj_df, requires_l2) {
  lg$warn("AFNI L3 setup is not implemented.")
  list(metadata = NULL, data = NULL, id_cols = NULL)
}


############
# helper function to get a cope data.frame for all level 1 models
get_l1_cope_df <- function(gpa, model_set, subj_df=NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(exclude_subject == FALSE) %>%
      dplyr::select(id, session)
  }
  checkmate::assert_data_frame(model_set)
  dt <- dplyr::bind_rows(
    lapply(unique(model_set$l1_model), function(mm) {
      data.frame(
        l1_model = mm,
        l1_cope_number = seq_along(gpa$l1_cope_names[[mm]]),
        l1_cope_name = gpa$l1_cope_names[[mm]]
      )
    })
  )

  #return subject-specific rows for each cope
  dt <- dt %>% tidyr::crossing(subj_df)
  return(dt)
}

# helper function to get a cope data.frame for all level 3 models
# handles the per-subject cope numbering problem
get_l2_cope_df <- function(gpa, model_set, subj_df=NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(exclude_subject == FALSE) %>%
      dplyr::select(id, session)
  }
  checkmate::assert_data_frame(model_set)

  dplyr::bind_rows(
    lapply(unique(model_set$l2_model), function(mm) {
      if (!is.null(gpa$l2_models$models[[mm]]$by_subject)) {
        # combine as single data frame from nested list columns
        l2_df <- dplyr::bind_rows(gpa$l2_models$models[[mm]]$by_subject$cope_list)
        l2_df$l2_model <- mm #retain model name
      } else {
        cope_names <- rownames(gpa$l2_models$models[[mm]]$contrasts)
        l2_df <- data.frame(
          l2_model = mm, l2_cope_number = seq_along(cope_names),
          l2_cope_name = cope_names
        )
        l2_df <- l2_df %>% tidyr::crossing(subj_df)
      }
      return(l2_df)

    })
  )
}

# helper function to get a cope data.frame for all level 3 models
get_l3_cope_df <- function(gpa, model_set, subj_df=NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(exclude_subject == FALSE) %>%
      dplyr::select(id, session)
  }
  checkmate::assert_data_frame(model_set)
  dt <- dplyr::bind_rows(
    lapply(unique(model_set$l3_model), function(mm) {
      cope_names <- rownames(gpa$l3_models$models[[mm]]$contrasts)
      data.frame(
        l3_model = mm,
        l3_cope_number = seq_along(cope_names),
        l3_cope_name = cope_names
      )
    })
  )

  # return subject-specific rows for each cope
  dt <- dt %>% tidyr::crossing(subj_df)
  return(dt)

}


# data.table cross-join (tidyr::crossing is a bit faster and already exists)
# https://stackoverflow.com/questions/10600060/how-to-do-cross-join-in-r
# CJ.table <- function(X, Y) {
#   setkey(X[, c(k = 1, .SD)], k)[Y[, c(k = 1, .SD)], allow.cartesian = TRUE][, k := NULL]
# }

get_fsl_l3_model_df <- function(gpa, model_df, subj_df=NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(exclude_subject == FALSE) %>%
      dplyr::select(id, session)
  }
  model_df$model_id <- seq_len(nrow(model_df))

  l1_df <- get_l1_cope_df(gpa, model_df, subj_df)

  l3_df <- get_l3_cope_df(gpa, model_df, subj_df)

  if (isTRUE(gpa$multi_run)) {
    # model_df has l1_model, l2_model, l3_model
    l2_df <- get_l2_cope_df(gpa, model_df)

    combined <- model_df %>%
      tidyr::crossing(subj_df) %>%
      left_join(l1_df, by = c("id", "session", "l1_model")) %>%
      left_join(l2_df, by = c("id", "session", "l2_model")) %>%
      left_join(l3_df, by = c("id", "session", "l3_model"))
  } else {
    combined <- model_df %>%
      tidyr::crossing(subj_df) %>%
      left_join(l1_df, by = c("id", "session", "l1_model")) %>%
      left_join(l3_df, by = c("id", "session", "l3_model"))
  }

  return(combined)
}

get_feat_l3_inputs <- function(gpa, l3_cope_config, lg=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(l3_cope_config)

  if (is.null(lg)) { lg <- lgr::get_logger() }
  if (isTRUE(gpa$multi_run)) {
    # feat directories in $l2_model_setup
    feat_inputs <- gpa$l2_model_setup$fsl %>%
      dplyr::filter(feat_complete == TRUE)

    #join up combination of all models with cope directories
    feat_inputs <- feat_inputs %>%
      dplyr::inner_join(l3_cope_config, by=c("id", "session", "l1_model", "l2_model"))

    # sort out expected cope files for each model combination
    feat_inputs <- feat_inputs %>%
      dplyr::mutate(cope_file = file.path(
        feat_dir,
        paste0("cope", l1_cope_number, ".feat"),
        "stats",
        paste0("cope", l2_cope_number, ".nii.gz")
      )) %>%
      dplyr::select(
        id, session, l1_model, l2_model, l3_model, l1_cope_name, l2_cope_name, feat_dir, cope_file
      )

    split_on <- c("l1_cope_name", "l2_cope_name", "l1_model", "l2_model", "l3_model")
  } else {
    # feat directories in $l1_model_setup
    feat_inputs <- gpa$l1_model_setup$fsl %>%
      dplyr::filter(feat_complete == TRUE)

    feat_inputs <- feat_inputs %>%
      dplyr::inner_join(l3_cope_config, by = c("id", "session", "l1_model"))

    feat_inputs <- feat_inputs %>%
      dplyr::mutate(cope_file = file.path(
        feat_dir,
        "stats",
        paste0("cope", l1_cope_number, ".nii.gz")
      )) %>%
      dplyr::select(
        id, session, l1_model, l3_model, l1_cope_name, feat_dir, cope_file
      )

    split_on <- c("l1_cope_name", "l1_model", "l3_model")

  }

  if (!is.data.table(feat_inputs)) data.table::setDT(feat_inputs)
  feat_inputs <- split(feat_inputs, by = split_on)

  return(feat_inputs)
}

get_spm_contrast_file <- function(spm_dir, cope_number) {
  base <- sprintf("con_%04d.nii", cope_number)
  cand <- file.path(spm_dir, base)
  if (file.exists(cand)) return(cand)
  cand_gz <- paste0(cand, ".gz")
  if (file.exists(cand_gz)) return(cand_gz)
  cand_gz2 <- sub("\\.nii$", ".nii.gz", cand)
  if (file.exists(cand_gz2)) return(cand_gz2)
  return(cand)
}

get_spm_l3_model_df <- function(gpa, model_df, subj_df = NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(exclude_subject == FALSE) %>%
      dplyr::select(id, session)
  }
  model_df$model_id <- seq_len(nrow(model_df))

  l1_df <- get_l1_cope_df(gpa, model_df, subj_df)
  l3_df <- get_l3_cope_df(gpa, model_df, subj_df)

  combined <- model_df %>%
    tidyr::crossing(subj_df) %>%
    left_join(l1_df, by = c("id", "session", "l1_model")) %>%
    left_join(l3_df, by = c("id", "session", "l3_model"))

  return(combined)
}

get_spm_l3_inputs <- function(gpa, l3_cope_config, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(l3_cope_config)

  if (is.null(lg)) { lg <- lgr::get_logger() }

  spm_inputs <- gpa$l1_model_setup$spm
  if (is.null(spm_inputs) || (is.data.frame(spm_inputs) && nrow(spm_inputs) == 0L)) {
    lg$warn("No SPM L1 model setup data found. Cannot setup SPM L3 models.")
    return(list())
  }

  if (!"spm_dir" %in% names(spm_inputs)) {
    lg$warn("No spm_dir column in gpa$l1_model_setup$spm. Cannot setup SPM L3 models.")
    return(list())
  }

  if ("spm_complete" %in% names(spm_inputs)) {
    spm_inputs <- spm_inputs %>%
      dplyr::filter(spm_complete == TRUE)
  } else {
    lg$warn("No spm_complete column in gpa$l1_model_setup$spm. Proceeding without completion filter for L3 inputs.")
  }

  if ("spm_contrast_exists" %in% names(spm_inputs)) {
    spm_inputs <- spm_inputs %>%
      dplyr::filter(spm_contrast_exists == TRUE)
  }

  spm_inputs <- spm_inputs %>%
    dplyr::inner_join(l3_cope_config, by = c("id", "session", "l1_model"))

  spm_inputs <- spm_inputs %>%
    dplyr::mutate(con_file = mapply(get_spm_contrast_file, spm_dir, l1_cope_number, USE.NAMES = FALSE)) %>%
    dplyr::select(id, session, l1_model, l3_model, l1_cope_name, spm_dir, con_file)

  spm_inputs <- spm_inputs %>%
    dplyr::filter(file.exists(con_file))

  if (!is.data.table(spm_inputs)) data.table::setDT(spm_inputs)
  spm_inputs <- split(spm_inputs, by = c("l1_cope_name", "l1_model", "l3_model"))

  return(spm_inputs)
}
