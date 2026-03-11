#' Function for setting up level 1 model specifications and corresponding files
#'
#' @param gpa A glm_model_arguments function setup by the \code{setup_glm_pipeline} function
#' @param l1_model_names A character vector of model names within \code{gpa$l1_models} whose l1
#'    inputs should be generated. If omitted, all models within \code{gpa$l1_models} will be generated.
#'
#' @details The \code{l1_model_names} argument allows the creation of multiple l1 models to be parallelized at
#'   a superordinate level or, even better, to be spread across independent jobs on a cluster. This function
#'   already provides the option to parallelize over subjects for a single model if \code{gpa$l1_setup_cpus}
#'   is greater than 1.
#'
#' @author Michael Hallquist
#'
#' @importFrom checkmate assert_class assert_character
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach registerDoSEQ
#' @importFrom iterators iter
#' @importFrom parallel makeCluster stopCluster
#' @importFrom RNifti niftiHeader
#' @importFrom dplyr filter pull
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @export
setup_l1_models <- function(gpa, l1_model_names=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_subset(l1_model_names, names(gpa$l1_models$models))
  checkmate::assert_data_frame(gpa$subject_data)
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_names(names(gpa$run_data), must.include = c("id", "session", "run_number", "run_nifti", "exclude_run")) # required columns

  #if no model subset is requested, output all models
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)

  lg <- lgr::get_logger("glm_pipeline/l1_setup")
  lg$set_threshold(gpa$lgr_threshold)

  if (isTRUE(gpa$log_json) && !"setup_l1_log_json" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderJson$new(gpa$output_locations$setup_l1_log_json), name = "setup_l1_log_json")
  }

  if (isTRUE(gpa$log_txt) && !"setup_l1_log_txt" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderFile$new(gpa$output_locations$setup_l1_log_txt), name = "setup_l1_log_txt")
  }

  gpa <- refresh_l1_cope_names(gpa, lg = lg)

  # TODO: This is a mess if use the l1_model_names since we will always be overwriting what's in the l1_model_setup
  # data.frame. We need more of an append/update approach, perhaps like the sqlite setup for the overall pipeline.

  # refresh status of feat inputs and outputs at level 1 before continuing
  gpa <- refresh_glm_status(gpa, level=1L, lg=lg)

  spm_l1_session_mode <- "separate"
  spm_anchor_by_id <- NULL
  if ("spm" %in% gpa$glm_software) {
    spm_mode_raw <- NULL
    if (!is.null(gpa$glm_settings) && !is.null(gpa$glm_settings$spm) &&
        "l1_session_mode" %in% names(gpa$glm_settings$spm)) {
      spm_mode_raw <- gpa$glm_settings$spm$l1_session_mode
    }
    spm_l1_session_mode <- resolve_spm_l1_session_mode(spm_mode_raw, lg = lg)
    if (identical(spm_l1_session_mode, "pooled")) {
      subj_ids <- unique(as.character(gpa$subject_data$id))
      spm_anchor_by_id <- setNames(
        vapply(subj_ids, function(sid) {
          min(as.integer(gpa$subject_data$session[gpa$subject_data$id == sid]), na.rm = TRUE)
        }, integer(1)),
        subj_ids
      )
      lg$info("SPM L1 session mode is 'pooled'; SPM L1 models will be created once per id across all sessions.")
    } else {
      lg$info("SPM L1 session mode is 'separate'; SPM L1 models will be created per id/session.")
    }
  }

  #setup parallel worker pool, if requested
  if (!is.null(gpa$parallel$l1_setup_cores) && gpa$parallel$l1_setup_cores[1L] > 1L) {
    lg$info("Initializing l1 setup cluster with %d cores", gpa$parallel$l1_setup_cores[1L])
    cl <- parallel::makeCluster(gpa$parallel$l1_setup_cores[1L])
    doParallel::registerDoParallel(cl)
    on.exit(try(parallel::stopCluster(cl))) #cleanup pool upon exit of this function
  } else {
    lg$info("Initializing l1 setup with serial execution")
    foreach::registerDoSEQ() #formally register a sequential 'pool' so that dopar is okay
  }

  # loop over each subject, identify relevant fMRI data, and setup level 1 analysis files
  # all_subj_l1_list <- lapply(seq_len(nrow(gpa$subject_data)), function(ii) {
  #   subj_df <- gpa$subject_data[ii, ]

  all_subj_l1_list <- foreach(
    subj_df = iter(gpa$subject_data, by = "row"), .inorder = FALSE, .packages = c("dplyr", "fmri.pipeline"), .errorhandling = "remove",
    .export = c(
      "truncate_runs", "fsl_l1_model", "spm_l1_model", "get_spm_status",
      "get_mr_abspath", "get_output_directory", "run_fsl_command",
      "get_feat_status", "add_custom_feat_syntax"
    )
  ) %dopar% {
    subj_df <- subj_df # avoid complaints about visible global binding in R CMD check
    subj_id <- subj_df$id
    subj_session <- subj_df$session
    glm_backends <- get_glm_backends(gpa)
    backend_fsl <- glm_backends[["fsl"]]
    backend_spm <- glm_backends[["spm"]]
    is_spm_anchor <- !identical(spm_l1_session_mode, "pooled") ||
      (!is.null(spm_anchor_by_id) &&
       as.character(subj_id) %in% names(spm_anchor_by_id) &&
       as.integer(subj_session) == as.integer(spm_anchor_by_id[[as.character(subj_id)]]))

    # find the run data for analysis
    rdata <- gpa$run_data %>% dplyr::filter(id == !!subj_id & session == !!subj_session & run_nifti_present == TRUE)
    run_nifti <- get_mr_abspath(rdata, "run_nifti")
    mr_run_nums <- rdata$run_number
    run_volumes <- rdata$run_volumes
    nvoxels <- rdata$nvoxels
    exclude_run <- rdata$exclude_run

    # ensure that we have a TR column and that TR does not vary by run
    if (!"tr" %in% names(rdata)) {
      lg$error("No tr column in run data for subject: %s, session: %d", subj_id, subj_session)
      return(NULL)
    }

    if (length(unique(rdata$tr)) > 1L) {
      lg$error("More than one TR value for runs within a subject. This is not currently supported! subject: %s, session: %d", subj_id, subj_session)
      return(NULL)
    }

    if ("l1_confound_file" %in% names(rdata)) {
      l1_confound_files <- rdata$l1_confound_file
    } else {
      l1_confound_files <- rep("", length(run_nifti))
    }

    subj_mr_dir <- subj_df$mr_dir

    # process files
    if (length(run_nifti) == 0L) {
      lg$warn("Unable to find any preprocessed fMRI files in dir: %s", subj_mr_dir)
      return(NULL)
    } else if (any(!file.exists(run_nifti))) {
      nii_present <- file.exists(run_nifti)
      lg$warn("Could not find some of the expected preprocessed fMRI files. These will be dropped.")
      lg$warn("Missing: %s", run_nifti[!nii_present])
      run_nifti <- run_nifti[nii_present]
      mr_run_nums <- mr_run_nums[nii_present]
      l1_confound_files <- l1_confound_files[nii_present]
      run_volumes <- run_volumes[nii_present]
      exclude_run <- exclude_run[nii_present]
    } else {
      lg$debug(paste("MR files to analyze:", run_nifti)) # log files that were found
    }

    lg$debug("Volumes in run_nifti: %s", paste(run_volumes, collapse = ", "))

    # get all events that pertain to this participant
    m_events <- data.table::rbindlist(
      lapply(gpa$l1_models$events, function(this_event) {
        this_event$data %>% dplyr::filter(id == !!subj_id & session == !!subj_session)
      })
    )

    # lookup subject directory for placing truncated files (now handled upstream)
    # subj_output_directory <- get_output_directory(id = subj_id, session = subj_session, gpa = gpa, create_if_missing = FALSE, what = "sub")

    # initialize mr data frame with elements of $run_data
    mr_df <- data.frame(
      id = subj_id, session = subj_session, run_number = mr_run_nums, run_nifti, l1_confound_file = l1_confound_files,
      run_volumes, exclude_run, row.names = NULL
    )

    # Precompute concatenated confounds for SPM/AFNI if requested
    l1_confound_files_spm <- l1_confound_files
    if ("spm" %in% gpa$glm_software && length(run_nifti) > 1L) {
      spm_defaults <- list(concatenate_runs = FALSE)
      spm_settings <- populate_defaults(gpa$glm_settings$spm, spm_defaults)
      if (is.null(spm_settings$concatenate_runs)) {
        spm_settings$concatenate_runs <- FALSE
      }

      if (isTRUE(spm_settings$concatenate_runs)) {
        has_nonempty_confounds <- length(l1_confound_files) > 0L &&
          any(!is.na(l1_confound_files) & nzchar(l1_confound_files))
        valid_confounds <- has_nonempty_confounds &&
          all(!is.na(l1_confound_files) & nzchar(l1_confound_files) & file.exists(l1_confound_files))

        if (isTRUE(valid_confounds)) {
          confound_outdir <- NULL
          if (length(l1_confound_files) > 0L && file.exists(l1_confound_files[1L])) {
            confound_outdir <- dirname(l1_confound_files[1L])
          } else {
            confound_outdir <- get_output_directory(
              id = subj_id, session = subj_session, gpa = gpa,
              what = "sub", create_if_missing = TRUE
            )
          }
          concat_file <- concat_l1_confounds(
            gpa = gpa, id = subj_id, session = subj_session,
            run_numbers = mr_run_nums, confound_files = l1_confound_files,
            output_dir = confound_outdir, lg = lg
          )
          if (!is.null(concat_file)) {
            l1_confound_files_spm <- concat_file
            lg$info("Using concatenated confounds for SPM: %s", concat_file)
          } else {
            l1_confound_files_spm <- NULL
            lg$warn(
              "Unable to concatenate confounds for SPM (id=%s session=%s); proceeding without confounds.",
              subj_id, subj_session
            )
          }
        } else if (!isTRUE(has_nonempty_confounds)) {
          l1_confound_files_spm <- NULL
          lg$info(
            "SPM concatenate_runs=TRUE with no confounds configured for id=%s session=%s; proceeding without confounds.",
            subj_id, subj_session
          )
        } else {
          l1_confound_files_spm <- NULL
          invalid <- l1_confound_files[is.na(l1_confound_files) | !nzchar(l1_confound_files) | !file.exists(l1_confound_files)]
          invalid_fmt <- if (length(invalid) == 0L) "<none>" else paste(ifelse(is.na(invalid), "NA", invalid), collapse = ", ")
          lg$warn(
            "SPM concatenate_runs=TRUE but confounds are missing/invalid for id=%s session=%s; proceeding without confounds. Missing/invalid: %s",
            subj_id, subj_session, invalid_fmt
          )
        }
      }
    }

    # Tracking list containing data.frames for each software, where we expect one row per run-level model (FSL)
    # or subject-level model (AFNI). The structure varies because FSL estimates per-run GLMs, while AFNI concatenates.
    l1_file_setup <- list(fsl = list(), spm = list(), afni = list(), metadata = mr_df)

    # loop over models to output
    for (ii in seq_along(l1_model_names)) {
      this_model <- l1_model_names[ii]

      # look for all ts_multiplier columns used in model signals
      ts_multiplier_cols <- unlist(lapply(gpa$l1_models$signals[gpa$l1_models$models[[this_model]]$signals], function(this_signal) {
        if (isFALSE(this_signal$ts_multiplier)) {
          return(NULL)
        } else {
          return(this_signal$ts_multiplier)
        }
      }))

      # pull corresponding ppi data for these columns, pass to BDM
      ts_multiplier_data <- if (!is.null(ts_multiplier_cols)) {
        gpa$ppi_data %>%
          dplyr::filter(id == !!subj_id & session == !!subj_session)
      } else {
        NULL
      }
        
      # setup design matrix for any given software package
      m_signals <- lapply(gpa$l1_models$signals[gpa$l1_models$models[[this_model]]$signals], function(this_signal) {
        # filter down to this id if the signal is a data.frame
        if (inherits(this_signal$value, "data.frame")) {
          this_signal$value <- this_signal$value %>%
            dplyr::filter(id == !!subj_id & session == !!subj_session)

          if (nrow(this_signal$value) == 0L) {
            msg <- glue(
              "In L1 model setup, failed to find any rows in in gpa$l1_models$signals${this_signal$name}$value",
              " for id: {subj_id}, session: {subj_session}, model: {this_model}.\n  We will create an empty regressor.",
              .trim = FALSE
            )
            lg$warn(msg)
            warning(msg)
          }

          # refit wi model if needed
          if (!is.null(this_signal$wi_model) && nrow(this_signal$value) > 0L) {
            this_signal <- fit_wi_model(this_signal)
          }
        } else {
          msg <- "Unable to sort out how to refit wi_model with signal that doesn't have a data.frame in the $value slot"
          lg$error(msg)
          stop(msg)
        }
        return(this_signal)
      })

      l1_output_dir <- get_output_directory(
        id = subj_id, session = subj_session,
        l1_model = this_model, gpa = gpa, what = "l1"
      )

      if (!dir.exists(l1_output_dir)) {
        lg$info("Creating subject output directory: %s", l1_output_dir)
        dir.create(l1_output_dir, showWarnings = FALSE, recursive = TRUE)
      }

      bdm_out_file <- file.path(l1_output_dir, paste0(gpa$l1_models$models[[this_model]]$name, "_bdm_setup.RData"))
      run_bdm <- TRUE
      if (file.exists(bdm_out_file)) {
        lg$info("Loading BDM info from extant file: %s", bdm_out_file)
        # Attempt to load BDM. If it fails, re-run build_design matrix
        run_bdm <- tryCatch(
          {
            load(bdm_out_file)
            if (!exists("d_obj") || is.null(d_obj)) {
              lg$warn("Regenerating design because d_obj is missing/NULL in extant BDM file: %s.", bdm_out_file)
              TRUE # regenerate
            } else {
              FALSE # use cache
            }
          },
          error = function(e) {
            lg$error("Failed to load BDM file: %s with error: %s. I will regenerate the design matrix.", bdm_out_file, e)
            return(TRUE)
          }
        )

        if ("run_4d_files" %in% names(d_obj)) { # older nomenclature ca. mid 2021
          d_obj$run_nifti <- d_obj$run_4d_files
        }
      }

      if (run_bdm) {
        t_out <- gpa$glm_software
        if (isTRUE(gpa$use_preconvolve)) t_out <- c("convolved", t_out) # compute preconvolved regressors
        bdm_args <- gpa$additional$bdm_args
        bdm_args$events <- m_events
        bdm_args$signals <- m_signals
        bdm_args$tr <- rdata$tr[1L]
        bdm_args$write_timing_files <- t_out
        bdm_args$drop_volumes <- gpa$drop_volumes
        bdm_args$run_data <- mr_df
        bdm_args$runs_to_output <- mr_run_nums
        bdm_args$output_directory <- file.path(l1_output_dir, "timing_files")
        bdm_args$lg <- lg
        bdm_args$ts_multipliers <- ts_multiplier_data

        d_obj <- tryCatch(do.call(build_design_matrix, bdm_args), error = function(e) {
          lg$error("Failed build_design_matrix for id: %s, session: %s, model: %s", subj_id, subj_session, this_model)
          lg$error("Error message: %s", as.character(e))
          return(NULL)
        })

        save(d_obj, bdm_args, mr_df, mr_run_nums, subj_mr_dir, run_nifti, run_volumes, exclude_run, nvoxels,
          l1_confound_files, subj_id, subj_session, this_model,
          file = bdm_out_file
        )

        if (is.null(d_obj)) next # skip to next iteration on error
      }

      if (!is.null(backend_fsl)) {
        # Setup FSL run-level models for each combination of signals
        # Returns a data.frame of feat l1 inputs and the fsf file
        feat_l1_df <- tryCatch(
          {
            backend_fsl$l1_setup(
              id = subj_id, session = subj_session, l1_confound_files = l1_confound_files, d_obj = d_obj,
              gpa = gpa, this_model, nvoxels = nvoxels
            )
          },
          error = function(e) {
            lg$error("Problem with fsl_l1_model. Model: %s, Subject: %s, Session: %s", this_model, subj_id, subj_session)
            lg$error("Error message: %s", as.character(e))
            return(NULL)
          }
        )

        if (!is.null(feat_l1_df)) {
          # add to tracking data.frame (simple, inefficient append)
          l1_file_setup$fsl <- rbind(l1_file_setup$fsl, feat_l1_df)
        }
      }

      if (is.null(backend_spm)) {
        lg$warn("SPM backend is not configured; skipping SPM L1 setup for subject %s session %s.", subj_id, subj_session)
      } else if (identical(spm_l1_session_mode, "pooled") && !isTRUE(is_spm_anchor)) {
        lg$debug(
          "Skipping SPM pooled L1 setup for non-anchor row id=%s session=%s model=%s (anchor session=%s).",
          subj_id, subj_session, this_model, spm_anchor_by_id[[as.character(subj_id)]]
        )
      } else if (!identical(spm_l1_session_mode, "pooled")) {
        lg$info("Setting up SPM L1 model %s for subject %s session %s.", this_model, subj_id, subj_session)
        spm_l1_df <- tryCatch(
          {
            backend_spm$l1_setup(
              id = subj_id, session = subj_session, l1_confound_files = l1_confound_files_spm,
              d_obj = d_obj, gpa = gpa, model_name = this_model, run_nifti = run_nifti,
              run_numbers = mr_run_nums, run_sessions = rep.int(as.integer(subj_session), length(mr_run_nums)),
              source_run_numbers = mr_run_nums, l1_session_mode = "separate"
            )
          },
          error = function(e) {
            lg$error("Problem running spm_l1_model. Model: %s, Subject: %s, Session: %s", this_model, subj_id, subj_session)
            lg$error("Error message: %s", as.character(e))
            return(NULL)
          }
        )

        if (!is.null(spm_l1_df)) {
          l1_file_setup$spm <- rbind(l1_file_setup$spm, spm_l1_df)
        } else {
          lg$warn("SPM L1 model setup returned NULL for model %s, subject %s, session %s.", this_model, subj_id, subj_session)
        }
      } else {
        lg$info("Setting up pooled-session SPM L1 model %s for subject %s.", this_model, subj_id)

        pooled_rdata <- gpa$run_data %>% dplyr::filter(id == !!subj_id)
        if ("run_nifti_present" %in% names(pooled_rdata)) {
          pooled_rdata <- pooled_rdata %>% dplyr::filter(run_nifti_present == TRUE)
        }
        pooled_rdata <- pooled_rdata[order(pooled_rdata$session, pooled_rdata$run_number), , drop = FALSE]

        if (nrow(pooled_rdata) == 0L) {
          lg$warn("No run data found for pooled SPM L1 setup: id=%s model=%s", subj_id, this_model)
          next
        }
        if (!"tr" %in% names(pooled_rdata) || length(unique(pooled_rdata$tr)) > 1L) {
          lg$error("SPM pooled L1 requires one TR per subject across sessions. id=%s model=%s", subj_id, this_model)
          next
        }

        pooled_run_nifti <- get_mr_abspath(pooled_rdata, "run_nifti")
        pooled_run_sessions <- as.integer(pooled_rdata$session)
        pooled_source_runs <- as.integer(pooled_rdata$run_number)
        pooled_run_numbers <- seq_along(pooled_source_runs)
        pooled_run_volumes <- pooled_rdata$run_volumes
        pooled_exclude_run <- pooled_rdata$exclude_run
        pooled_l1_confound_files <- if ("l1_confound_file" %in% names(pooled_rdata)) pooled_rdata$l1_confound_file else rep("", length(pooled_run_numbers))

        run_key <- paste(pooled_run_sessions, pooled_source_runs, sep = "::")
        remap_run_numbers <- function(df, what = "data") {
          if (!is.data.frame(df) || nrow(df) == 0L) return(df)
          if (!all(c("session", "run_number") %in% names(df))) return(df)
          idx <- match(paste(df$session, df$run_number, sep = "::"), run_key)
          if (anyNA(idx)) {
            lg$warn(
              "Dropping %d row(s) with unmatched session/run_number while remapping pooled SPM %s for id=%s model=%s.",
              sum(is.na(idx)), what, subj_id, this_model
            )
            df <- df[!is.na(idx), , drop = FALSE]
            idx <- idx[!is.na(idx)]
          }
          df$run_number <- pooled_run_numbers[idx]
          df
        }

        m_events_spm <- data.table::rbindlist(
          lapply(gpa$l1_models$events, function(this_event) {
            this_event$data %>% dplyr::filter(id == !!subj_id)
          })
        )
        m_events_spm <- remap_run_numbers(as.data.frame(m_events_spm), what = "events")

        m_signals_spm <- lapply(gpa$l1_models$signals[gpa$l1_models$models[[this_model]]$signals], function(this_signal) {
          if (inherits(this_signal$value, "data.frame")) {
            this_signal$value <- this_signal$value %>% dplyr::filter(id == !!subj_id)
            this_signal$value <- remap_run_numbers(this_signal$value, what = paste0("signal ", this_signal$name))
            if (!is.null(this_signal$wi_model) && nrow(this_signal$value) > 0L) {
              this_signal <- fit_wi_model(this_signal)
            }
          } else {
            msg <- "Unable to sort out how to refit wi_model with signal that doesn't have a data.frame in the $value slot"
            lg$error(msg)
            stop(msg)
          }
          this_signal
        })

        ts_multiplier_data_spm <- if (!is.null(ts_multiplier_cols)) {
          gpa$ppi_data %>% dplyr::filter(id == !!subj_id)
        } else {
          NULL
        }
        if (is.data.frame(ts_multiplier_data_spm) && nrow(ts_multiplier_data_spm) > 0L) {
          ts_multiplier_data_spm <- remap_run_numbers(ts_multiplier_data_spm, what = "ts_multipliers")
        }

        mr_df_spm <- data.frame(
          id = subj_id,
          session = pooled_run_sessions,
          run_number = pooled_run_numbers,
          run_nifti = pooled_run_nifti,
          l1_confound_file = pooled_l1_confound_files,
          run_volumes = pooled_run_volumes,
          exclude_run = pooled_exclude_run,
          row.names = NULL
        )

        spm_outdir <- get_output_directory(
          id = subj_id, session = 0L, l1_model = this_model,
          gpa = gpa, glm_software = "spm", what = "l1", create_if_missing = TRUE
        )
        bdm_args_spm <- gpa$additional$bdm_args
        t_out_spm <- "spm"
        if (isTRUE(gpa$use_preconvolve)) t_out_spm <- c("convolved", t_out_spm)
        bdm_args_spm$events <- m_events_spm
        bdm_args_spm$signals <- m_signals_spm
        bdm_args_spm$tr <- pooled_rdata$tr[1L]
        bdm_args_spm$write_timing_files <- t_out_spm
        bdm_args_spm$drop_volumes <- gpa$drop_volumes
        bdm_args_spm$run_data <- mr_df_spm
        bdm_args_spm$runs_to_output <- pooled_run_numbers
        bdm_args_spm$output_directory <- file.path(spm_outdir, "timing_files")
        bdm_args_spm$lg <- lg
        bdm_args_spm$ts_multipliers <- ts_multiplier_data_spm

        d_obj_spm <- tryCatch(do.call(build_design_matrix, bdm_args_spm), error = function(e) {
          lg$error("Failed pooled SPM build_design_matrix for id: %s, model: %s", subj_id, this_model)
          lg$error("Error message: %s", as.character(e))
          return(NULL)
        })
        if (is.null(d_obj_spm)) next

        spm_settings <- populate_defaults(gpa$glm_settings$spm, list(concatenate_runs = FALSE))
        l1_confound_files_spm_use <- pooled_l1_confound_files
        if (isTRUE(spm_settings$concatenate_runs)) {
          has_nonempty_confounds <- length(pooled_l1_confound_files) > 0L &&
            any(!is.na(pooled_l1_confound_files) & nzchar(pooled_l1_confound_files))
          valid_confounds <- has_nonempty_confounds &&
            all(!is.na(pooled_l1_confound_files) & nzchar(pooled_l1_confound_files) & file.exists(pooled_l1_confound_files))

          if (isTRUE(valid_confounds)) {
            concat_file <- concat_l1_confounds(
              gpa = gpa,
              id = subj_id,
              session = 0L,
              run_numbers = pooled_source_runs,
              run_sessions = pooled_run_sessions,
              confound_files = pooled_l1_confound_files,
              output_dir = spm_outdir,
              file_name = paste0("l1_confounds_concat_id", subj_id, ".txt"),
              lg = lg
            )
            if (is.null(concat_file)) {
              l1_confound_files_spm_use <- NULL
              lg$warn(
                "Unable to concatenate confounds for pooled SPM setup (id=%s model=%s); proceeding without confounds.",
                subj_id, this_model
              )
            } else {
              l1_confound_files_spm_use <- concat_file
            }
          } else if (!isTRUE(has_nonempty_confounds)) {
            l1_confound_files_spm_use <- NULL
            lg$info(
              "SPM concatenate_runs=TRUE with no pooled confounds configured for id=%s model=%s; proceeding without confounds.",
              subj_id, this_model
            )
          } else {
            l1_confound_files_spm_use <- NULL
            invalid <- pooled_l1_confound_files[
              is.na(pooled_l1_confound_files) |
                !nzchar(pooled_l1_confound_files) |
                !file.exists(pooled_l1_confound_files)
            ]
            invalid_fmt <- if (length(invalid) == 0L) "<none>" else paste(ifelse(is.na(invalid), "NA", invalid), collapse = ", ")
            lg$warn(
              "SPM concatenate_runs=TRUE but pooled confounds are missing/invalid for id=%s model=%s; proceeding without confounds. Missing/invalid: %s",
              subj_id, this_model, invalid_fmt
            )
          }
        }

        spm_l1_df <- tryCatch(
          {
            backend_spm$l1_setup(
              id = subj_id, session = 0L, l1_confound_files = l1_confound_files_spm_use,
              d_obj = d_obj_spm, gpa = gpa, model_name = this_model, run_nifti = pooled_run_nifti,
              run_numbers = pooled_run_numbers, run_sessions = pooled_run_sessions,
              source_run_numbers = pooled_source_runs, l1_session_mode = "pooled"
            )
          },
          error = function(e) {
            lg$error("Problem running pooled spm_l1_model. Model: %s, Subject: %s", this_model, subj_id)
            lg$error("Error message: %s", as.character(e))
            return(NULL)
          }
        )

        if (!is.null(spm_l1_df)) {
          l1_file_setup$spm <- rbind(l1_file_setup$spm, spm_l1_df)
        } else {
          lg$warn("SPM pooled L1 model setup returned NULL for model %s, subject %s.", this_model, subj_id)
        }
      }
    }

    lg$info("Completed processing of subject: %s", subj_id)
    return(l1_file_setup)
  }

  # N.B. rbindlist converts the bound elements into a single data.table object
  spm_list <- lapply(all_subj_l1_list, "[[", "spm")
  spm_combined <- NULL
  if (!all(vapply(spm_list, is.null, logical(1)))) {
    spm_combined <- rbindlist(spm_list, fill = TRUE)
  }

  all_subj_l1_combined <- list(
    fsl = rbindlist(lapply(all_subj_l1_list, "[[", "fsl"), fill = TRUE),
    spm = spm_combined,
    metadata = rbindlist(lapply(all_subj_l1_list, "[[", "metadata"), fill = TRUE)
  )

  class(all_subj_l1_combined) <- c("l1_setup", "list")

  #append l1 setup to gpa
  gpa$l1_model_setup <- all_subj_l1_combined
  
  # refresh l1 model status in $l1_model_setup
  gpa <- refresh_glm_status(gpa, level=1L, lg=lg)

  return(gpa)

}
