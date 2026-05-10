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

  glm_backends <- get_glm_backends(gpa, level = 2L)
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

  # Filter to backends that support standalone L2 estimation (currently only FSL)
  # SPM uses L2 model specs during L1 setup for run/session projection
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

  # normalize model signatures so l2_scope is always available and validated
  gpa <- normalize_longitudinal_model_signatures(gpa, lg)
  gpa <- refresh_l1_cope_names(gpa, lg = lg)

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
    lg$info(
      "No backend with standalone L2 estimation requested. FSL estimates L2 directly; SPM consumes L2 model specs during L1 setup."
    )
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

  # compose L2 modeling data once so builder/respecification use the same columns
  l2_data <- compose_l2_model_data(gpa, lg = lg)

  # only retain good runs and subjects
  run_data <- l2_data %>%
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
    this_scope <- gpa$l2_models$models[[mname]]$l2_scope
    if (is.null(this_scope) || !nzchar(this_scope)) this_scope <- "id_session"
    split_on <- if (identical(this_scope, "id")) "id" else c("id", "session")
    lg$info(
      "Recalculating per-subject L2 models based on available runs for model: %s (l2_scope=%s)",
      mname, this_scope
    )
    gpa$l2_models$models[[mname]] <- respecify_l2_models_by_subject(
      gpa$l2_models$models[[mname]],
      run_data,
      split_on = split_on,
      aggregated_session = 0L
    )
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

  l1_cope_validity <- setNames(lapply(backend_results, `[[`, "l1_cope_validity"), names(backend_results))
  l1_cope_validity <- l1_cope_validity[!vapply(l1_cope_validity, is.null, logical(1))]
  if (length(l1_cope_validity) > 0L) {
    class(l1_cope_validity) <- c("l1_cope_validity", "list")
    gpa$l1_cope_validity <- l1_cope_validity

    validity_tsv <- gpa$output_locations$setup_l2_l1_cope_validity_tsv
    if (is.null(validity_tsv) || !nzchar(validity_tsv)) {
      validity_tsv <- file.path(gpa$output_directory, "setup_l2_l1_cope_validity.tsv")
    }
    validity_out <- dplyr::bind_rows(l1_cope_validity, .id = "backend")
    if (nrow(validity_out) > 0L) {
      dir.create(dirname(validity_tsv), recursive = TRUE, showWarnings = FALSE)
      utils::write.table(validity_out, file = validity_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }

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

is_all_zero_numeric <- function(x, tol = sqrt(.Machine$double.eps)) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(TRUE)
  all(abs(x) <= tol)
}

parse_fsl_lvl1_ev_map <- function(feat_fsf) {
  checkmate::assert_file_exists(feat_fsf)
  fsf_lines <- readLines(feat_fsf, warn = FALSE)

  extract_field <- function(pattern) {
    hits <- grep(pattern, fsf_lines, value = TRUE, perl = TRUE)
    if (length(hits) == 0L) {
      return(data.frame(index = integer(0), value = character(0), stringsAsFactors = FALSE))
    }

    parsed <- regmatches(hits, regexec(pattern, hits, perl = TRUE))
    parsed <- parsed[lengths(parsed) > 0L]
    data.frame(
      index = as.integer(vapply(parsed, `[[`, character(1), 2L)),
      value = vapply(parsed, `[[`, character(1), 3L),
      stringsAsFactors = FALSE
    )
  }

  titles <- extract_field("^set fmri\\(evtitle([0-9]+)\\) \"(.*)\"$")
  shapes <- extract_field("^set fmri\\(shape([0-9]+)\\) ([0-9]+)$")
  customs <- extract_field("^set fmri\\(custom([0-9]+)\\) \"(.*)\"$")

  out <- dplyr::full_join(titles, shapes, by = "index", suffix = c("_title", "_shape")) %>%
    dplyr::full_join(customs, by = "index")
  if (nrow(out) == 0L) {
    return(data.frame(index = integer(0), name = character(0), shape = integer(0), timing_file = character(0)))
  }

  out <- out %>%
    dplyr::transmute(
      index = .data$index,
      name = dplyr::coalesce(.data$value_title, ""),
      shape = suppressWarnings(as.integer(.data$value_shape)),
      timing_file = dplyr::coalesce(.data$value, "")
    ) %>%
    dplyr::arrange(.data$index)

  out
}

detect_fsl_custom_ev_empty <- function(timing_file, shape_code) {
  if (is.na(shape_code)) return(FALSE)
  if (identical(shape_code, 10L)) return(TRUE)
  if (!nzchar(timing_file) || !file.exists(timing_file)) return(FALSE)

  if (identical(shape_code, 2L)) {
    vals <- tryCatch(scan(timing_file, what = numeric(), quiet = TRUE), error = function(e) numeric(0))
    if (length(vals) == 0L) {
      raw_vals <- tryCatch(scan(timing_file, what = character(), quiet = TRUE), error = function(e) character(0))
      vals <- suppressWarnings(as.numeric(raw_vals))
    }
    return(is_all_zero_numeric(vals))
  }

  if (identical(shape_code, 3L)) {
    ev_df <- tryCatch(utils::read.table(timing_file, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(ev_df) || nrow(ev_df) == 0L) return(TRUE)
    value_col <- if (ncol(ev_df) >= 3L) 3L else ncol(ev_df)
    vals <- suppressWarnings(as.numeric(ev_df[[value_col]]))
    return(is_all_zero_numeric(vals))
  }

  FALSE
}

infer_empty_ev_reason <- function(run_number, ev_name, timing_file, shape_code) {
  if (is.na(shape_code)) return("unknown")
  if (identical(shape_code, 10L)) return("empty_ev_shape")
  if (!nzchar(timing_file) || !file.exists(timing_file)) return("missing_timing_file")

  raw_3col <- file.path(dirname(timing_file), paste0("run", run_number, "_", ev_name, "_FSL3col.txt"))
  if (file.exists(raw_3col)) {
    raw_df <- tryCatch(utils::read.table(raw_3col, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(raw_df)) {
      return("unreadable_raw_fsl3col")
    }
    if (nrow(raw_df) == 0L) {
      return("empty_raw_events")
    }
    if (ncol(raw_df) >= 3L) {
      vals_chr <- trimws(as.character(raw_df[[3L]]))
      vals_num <- suppressWarnings(as.numeric(vals_chr))
      finite_vals <- vals_num[is.finite(vals_num)]
      if (length(vals_chr) > 0L && length(finite_vals) == 0L) {
        return("all_missing_parametric_modulator")
      }
      if (length(finite_vals) > 0L && is_all_zero_numeric(finite_vals)) {
        return("all_zero_raw_values")
      }
    }
  }

  "all_zero_timing_file"
}

empty_ev_reason_is_expected <- function(reason) {
  reason %in% c("empty_ev_shape", "empty_raw_events")
}

infer_l1_empty_ev_info_for_row <- function(l1_row, gpa, lg = NULL) {
  if (is.null(lg)) {
    lg <- lgr::get_logger("glm_pipeline/l2_setup")
  }

  l1_model <- l1_row$l1_model[[1L]]
  feat_fsf <- l1_row$feat_fsf[[1L]]
  run_number <- l1_row$run_number[[1L]]
  contrasts <- gpa$l1_models$models[[l1_model]]$contrasts

  if (is.null(contrasts)) {
    return(list(
      empty_evs = character(0),
      empty_ev_reasons = character(0)
    ))
  }
  if (checkmate::test_data_frame(contrasts)) contrasts <- as.matrix(contrasts)
  if (!is.matrix(contrasts) || nrow(contrasts) == 0L) {
    return(list(
      empty_evs = character(0),
      empty_ev_reasons = character(0)
    ))
  }

  if (!checkmate::test_file_exists(feat_fsf)) {
    lg$warn("Cannot inspect L1 FEAT FSF for contrast validity because it is missing: %s", feat_fsf)
    return(list(
      empty_evs = character(0),
      empty_ev_reasons = character(0)
    ))
  }

  ev_map <- parse_fsl_lvl1_ev_map(feat_fsf)
  if (nrow(ev_map) == 0L) {
    return(list(
      empty_evs = character(0),
      empty_ev_reasons = character(0)
    ))
  }

  ev_is_empty <- mapply(
    detect_fsl_custom_ev_empty,
    timing_file = ev_map$timing_file,
    shape_code = ev_map$shape,
    USE.NAMES = FALSE
  )
  empty_evs <- ev_map$name[ev_is_empty]
  if (length(empty_evs) == 0L) {
    return(list(
      empty_evs = character(0),
      empty_ev_reasons = character(0)
    ))
  }

  ev_reasons <- setNames(
    vapply(which(ev_is_empty), function(ii) {
      infer_empty_ev_reason(run_number, ev_map$name[ii], ev_map$timing_file[ii], ev_map$shape[ii])
    }, character(1)),
    ev_map$name[ev_is_empty]
  )

  list(
    empty_evs = empty_evs,
    empty_ev_reasons = ev_reasons
  )
}

build_l1_cope_validity <- function(l1_setup_fsl, gpa, lg = NULL) {
  checkmate::assert_data_frame(l1_setup_fsl)
  if (is.null(lg)) {
    lg <- lgr::get_logger("glm_pipeline/l2_setup")
  }

  checkmate::assert_subset(c("id", "session", "run_number", "l1_model", "feat_fsf", "feat_dir"), names(l1_setup_fsl))
  valid_states <- c("valid", "invalid_expected_empty", "invalid_unexpected_empty", "invalid_missing_file")

  validity_rows <- lapply(seq_len(nrow(l1_setup_fsl)), function(ii) {
    l1_row <- l1_setup_fsl[ii, , drop = FALSE]
    l1_model <- l1_row$l1_model[[1L]]
    cope_names <- gpa$l1_cope_names[[l1_model]]
    if (is.null(cope_names) || length(cope_names) == 0L) {
      cope_names <- extract_contrast_names(
        gpa$l1_models$models[[l1_model]]$contrasts,
        model_name = l1_model,
        level_label = "L1",
        allow_empty = TRUE
      )
    }
    if (length(cope_names) == 0L) {
      return(NULL)
    }

    contrasts <- gpa$l1_models$models[[l1_model]]$contrasts
    if (checkmate::test_data_frame(contrasts)) contrasts <- as.matrix(contrasts)
    empty_info <- infer_l1_empty_ev_info_for_row(l1_row, gpa, lg = lg)

    do_check_files <- "feat_complete" %in% names(l1_row) && isTRUE(l1_row$feat_complete[[1L]])

    dplyr::bind_rows(lapply(seq_along(cope_names), function(cc) {
      l1_cope_name <- cope_names[cc]
      cope_file <- file.path(l1_row$feat_dir[[1L]], "stats", paste0("cope", cc, ".nii.gz"))
      varcope_file <- file.path(l1_row$feat_dir[[1L]], "stats", paste0("varcope", cc, ".nii.gz"))
      missing_files <- character(0)
      if (isTRUE(do_check_files)) {
        if (!file.exists(cope_file)) missing_files <- c(missing_files, "cope_file")
        if (!file.exists(varcope_file)) missing_files <- c(missing_files, "varcope_file")
      }

      contrast_evs <- character(0)
      if (is.matrix(contrasts) && nrow(contrasts) >= cc) {
        contrast_evs <- colnames(contrasts)[which(abs(contrasts[cc, ]) > 0)]
      }
      empty_for_cope <- intersect(contrast_evs, empty_info$empty_evs)
      empty_reasons <- empty_info$empty_ev_reasons[empty_for_cope]
      empty_reasons <- empty_reasons[!is.na(empty_reasons)]

      validity_state <- "valid"
      validity_reason <- NA_character_
      if (length(missing_files) > 0L) {
        validity_state <- "invalid_missing_file"
        validity_reason <- paste0("missing_", paste(missing_files, collapse = ";missing_"))
      } else if (length(empty_for_cope) > 0L) {
        if (length(empty_reasons) > 0L && all(empty_ev_reason_is_expected(unique(empty_reasons)))) {
          validity_state <- "invalid_expected_empty"
        } else {
          validity_state <- "invalid_unexpected_empty"
        }
        validity_reason <- paste(unique(empty_reasons), collapse = ";")
      }

      stopifnot(validity_state %in% valid_states)
      data.frame(
        id = l1_row$id[[1L]],
        session = l1_row$session[[1L]],
        run_number = l1_row$run_number[[1L]],
        l1_model = l1_model,
        l1_cope_number = cc,
        l1_cope_name = l1_cope_name,
        validity_state = validity_state,
        valid_for_higher_level = identical(validity_state, "valid"),
        validity_reason = validity_reason,
        empty_evs = paste(empty_for_cope, collapse = ";"),
        missing_files = paste(missing_files, collapse = ";"),
        cope_file = cope_file,
        varcope_file = varcope_file,
        feat_fsf = l1_row$feat_fsf[[1L]],
        feat_dir = l1_row$feat_dir[[1L]],
        stringsAsFactors = FALSE
      )
    }))
  })

  out <- dplyr::bind_rows(validity_rows)
  if (nrow(out) == 0L) {
    return(data.frame(
      id = character(0), session = integer(0), run_number = integer(0),
      l1_model = character(0), l1_cope_number = integer(0), l1_cope_name = character(0),
      validity_state = character(0), valid_for_higher_level = logical(0),
      validity_reason = character(0), empty_evs = character(0), missing_files = character(0),
      cope_file = character(0), varcope_file = character(0),
      feat_fsf = character(0), feat_dir = character(0),
      stringsAsFactors = FALSE
    ))
  }

  out
}

is_strict_intercept_only_l2_model <- function(mobj) {
  if (is.null(mobj)) return(FALSE)

  mm <- mobj$model_matrix
  intercept_only_matrix <- is.matrix(mm) &&
    ncol(mm) == 1L &&
    identical(colnames(mm), "(Intercept)")

  if (!isTRUE(intercept_only_matrix) && !is.null(mobj$lmfit)) {
    model_terms <- tryCatch(stats::terms(mobj$lmfit), error = function(e) NULL)
    intercept_only_matrix <- !is.null(model_terms) &&
      identical(attr(model_terms, "intercept"), 1L) &&
      length(attr(model_terms, "term.labels")) == 0L
  }

  if (!isTRUE(intercept_only_matrix)) return(FALSE)

  cmat <- mobj$contrasts
  if (checkmate::test_data_frame(cmat)) cmat <- as.matrix(cmat)
  is.matrix(cmat) &&
    nrow(cmat) == 1L &&
    ncol(cmat) == 1L &&
    identical(colnames(cmat), "(Intercept)") &&
    isTRUE(all.equal(as.numeric(cmat[1L, 1L]), 1))
}

make_l2_passthrough_row <- function(l1_df, l2_model, gpa, lg = NULL) {
  if (is.null(lg)) {
    lg <- lgr::get_logger("glm_pipeline/l2_setup")
  }

  checkmate::assert_data_frame(l1_df)
  stopifnot(nrow(l1_df) == 1L)
  l2_scope <- gpa$l2_models$models[[l2_model]]$l2_scope
  if (is.null(l2_scope) || !is.character(l2_scope) || length(l2_scope) != 1L || !nzchar(l2_scope)) {
    l2_scope <- "id_session"
  }

  id <- l1_df$id[1L]
  session <- if (identical(l2_scope, "id")) 0L else l1_df$session[1L]
  l1_model <- l1_df$l1_model[1L]
  l1_cope_number <- l1_df$l1_cope_number[1L]
  l1_cope_name <- l1_df$l1_cope_name[1L]
  passthrough_cope_file <- l1_df$cope_file[1L]

  contrast_names <- extract_contrast_names(
    gpa$l2_models$models[[l2_model]]$contrasts,
    model_name = l2_model,
    level_label = "L2",
    allow_empty = FALSE
  )
  cope_list <- data.frame(
    id = id,
    session = session,
    l2_cope_number = 1L,
    l2_cope_name = contrast_names[1L],
    stringsAsFactors = FALSE
  )

  feat_l2_df <- data.frame(
    id = id, session = session,
    l1_model = l1_model,
    l1_cope_number = l1_cope_number,
    l1_cope_name = l1_cope_name,
    l2_model = l2_model,
    l2_scope = l2_scope,
    l2_input_mode = "l1_cope_file_passthrough",
    l2_passthrough = TRUE,
    n_l2_copes = 1L,
    n_input_files = 1L,
    passthrough_cope_file = passthrough_cope_file,
    stringsAsFactors = FALSE
  )
  feat_l2_df$cope_list <- list(cope_list)

  status <- get_feat_status(
    feat_dir = l1_df$feat_dir[1L],
    feat_fsf = l1_df$feat_fsf[1L],
    lg = lg
  )
  status$feat_complete <- isTRUE(status$feat_complete) && file.exists(passthrough_cope_file)
  status$feat_failed <- if (isTRUE(status$feat_complete)) FALSE else status$feat_failed

  feat_l2_df <- dplyr::bind_cols(feat_l2_df, status)
  feat_l2_df$to_run <- FALSE

  lg$warn(
    paste(
      "Only one valid run remains for subject %s session %s, L1 model %s, L1 cope %s, L2 model %s.",
      "Using direct L1 cope pass-through instead of creating a one-input L2 FEAT."
    ),
    id, session, l1_model, l1_cope_name, l2_model
  )

  feat_l2_df
}

setup_l2_backend_fsl <- function(gpa, backend, lg, l1_model_names, l2_model_names, good_runs) {
  if (is.null(backend)) {
    return(list(data = NULL, id_cols = NULL))
  }

  model_set <- expand.grid(l1_model = l1_model_names, l2_model = l2_model_names, stringsAsFactors = FALSE)
  all_l2_list <- foreach(
    model_info = iter(model_set, by = "row"), .inorder = FALSE,
    .packages = c("fmri.pipeline", "dplyr", "data.table")
  ) %dopar% {
    model_info <- model_info # avoid complaints about visible global binding in R CMD check
    this_l1_model <- model_info$l1_model
    this_l2_model <- model_info$l2_model
    this_scope <- gpa$l2_models$models[[this_l2_model]]$l2_scope
    if (is.null(this_scope) || !nzchar(this_scope)) this_scope <- "id_session"
    split_cols <- if (identical(this_scope, "id")) "id" else c("id", "session")

    l2_file_setup <- list(fsl = list(), l1_cope_validity = list())

    # get list of runs to examine/include
    to_run <- gpa$l1_model_setup$fsl %>%
      dplyr::filter(l1_model == !!this_l1_model) %>%
      dplyr::select(
        id, session, run_number, l1_model, feat_fsf, feat_dir,
        dplyr::any_of(c("feat_complete", "feat_failed"))
      )

    # handle run and subject exclusions by joining against good runs
    to_run <- dplyr::inner_join(good_runs, to_run, by = c("id", "session", "run_number"))
    if (nrow(to_run) == 0L) {
      lg$warn(
        "No eligible L1 runs remain before L2 setup for L1 model %s and L2 model %s after run/subject exclusions.",
        this_l1_model, this_l2_model
      )
      return(l2_file_setup)
    }

    cope_names <- gpa$l1_cope_names[[this_l1_model]]
    if (is.null(cope_names) || length(cope_names) == 0L) {
      lg$warn("No L1 cope names found for model '%s'; skipping L2 setup.", this_l1_model)
      return(l2_file_setup)
    }

    l1_cope_validity <- build_l1_cope_validity(to_run, gpa, lg = lg) %>%
      dplyr::mutate(l2_model = this_l2_model)
    l2_file_setup$l1_cope_validity <- l1_cope_validity

    dropped_inputs <- l1_cope_validity %>% dplyr::filter(.data$valid_for_higher_level == FALSE)
    if (nrow(dropped_inputs) > 0L) {
      for (ii in seq_len(nrow(dropped_inputs))) {
        rr <- dropped_inputs[ii, , drop = FALSE]
        lg$warn(
          "Dropping run %s for subject %s session %s, L1 model %s, L1 cope %s before L2 setup because the lower-level cope is %s (%s).",
          rr$run_number[1], rr$id[1], rr$session[1], rr$l1_model[1], rr$l1_cope_name[1],
          rr$validity_state[1],
          ifelse(is.na(rr$validity_reason[1]), "unknown", rr$validity_reason[1])
        )
      }
    }

    all_groups <- split(
      data.table::as.data.table(l1_cope_validity),
      by = c(split_cols, "l1_cope_name", "l1_cope_number")
    )
    for (cand in all_groups) {
      if (!all(cand$valid_for_higher_level == FALSE)) next
      subj_id <- cand$id[1L]
      subj_session <- if ("session" %in% names(cand) && length(unique(cand$session)) == 1L) cand$session[1L] else 0L
      reason_txt <- paste(unique(stats::na.omit(cand$validity_reason)), collapse = ";")
      if (!nzchar(reason_txt)) reason_txt <- "unknown"
      state_txt <- paste(unique(cand$validity_state), collapse = ";")
      lg$warn(
        "No valid runs remain for subject %s session %s, L1 model %s, L1 cope %s, L2 model %s. Skipping this L2 analysis. Dropped runs: %s. States: %s. Reasons: %s.",
        subj_id, subj_session, this_l1_model, cand$l1_cope_name[1L], this_l2_model,
        paste(cand$run_number, collapse = ","),
        state_txt,
        reason_txt
      )
    }

    to_run <- l1_cope_validity %>%
      dplyr::filter(.data$valid_for_higher_level == TRUE) %>%
      dplyr::mutate(
        l2_model = this_l2_model
      ) %>%
      dplyr::select(
        id, session, run_number, l1_model, l1_cope_number, l1_cope_name, l2_model, feat_fsf, feat_dir, cope_file
      )

    if (nrow(to_run) == 0L) {
      return(l2_file_setup)
    }

    data.table::setDT(to_run) # convert to data.table for split
    by_subj_session <- split(to_run, by = c(split_cols, "l1_cope_name", "l1_cope_number"))

    # setup Feat L2 files for each id and session
    for (l1_df in by_subj_session) {
      subj_id <- l1_df$id[1L]
      subj_session <- if (length(unique(l1_df$session)) == 1L) l1_df$session[1L] else 0L
      feat_l2_df <- tryCatch({
        backend$l2_setup(
          l1_df = l1_df,
          l2_model = this_l2_model, gpa = gpa
        )},
        error = function(e) {
          lg$error(
            "Problem with fsl_l2_model. L1 Model: %s, L2 Model: %s, Subject: %s, Session: %s, Scope: %s",
            this_l1_model, this_l2_model, subj_id, subj_session, this_scope
          )
          lg$error("Error message: %s", as.character(e))
          return(NULL)
        }
      )

      if (!is.null(feat_l2_df)) {
        # add to tracking data.frame
        l2_file_setup$fsl <- dplyr::bind_rows(l2_file_setup$fsl, feat_l2_df)
      }
    }

    return(l2_file_setup)
  }

  fsl_df <- dplyr::bind_rows(lapply(all_l2_list, "[[", "fsl"))
  validity_df <- dplyr::bind_rows(lapply(all_l2_list, "[[", "l1_cope_validity"))
  list(
    data = fsl_df,
    l1_cope_validity = validity_df,
    id_cols = c("id", "session", "l1_model", "l1_cope_name", "l2_model")
  )
}

setup_l2_backend_spm <- function(gpa, backend, lg, l1_model_names, l2_model_names, good_runs) {
  lg$info("SPM has no standalone L2 setup step; L2 model specs are applied during SPM L1 setup.")
  list(data = NULL, id_cols = NULL)
}

setup_l2_backend_afni <- function(gpa, backend, lg, l1_model_names, l2_model_names, good_runs) {
  lg$warn("afni not supported in setup_l2_models")
  list(data = NULL, id_cols = NULL)
}
