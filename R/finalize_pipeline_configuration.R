#' This is a small helper function to validate the glm_model_arguments list structure.
#' It adds a few details such as the output directory to make it less burdensome for to setup a pipeline
#' N.B. gpa is a shorthand abbreviation for glm_model_arguments, to save typing
#'
#' @param gpa A \code{glm_pipeline_arguments} object setup by \code{setup_glm_pipeline}
#' @importFrom stringr str_count fixed
#' @importFrom magrittr %>%
#' @importFrom lgr get_logger
#' @importFrom rhdf h5save
finalize_pipeline_configuration <- function(gpa) {

  lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")

  #sort out file locations
  if (is.null(gpa$group_output_directory) || gpa$group_output_directory == "default") {
    gpa$group_output_directory <- file.path(getwd(), "group_analyses", gpa$analysis_name)
  }

  gpa$sqlite_db <- file.path(gpa$working_directory, paste0(gpa$analysis_name, ".sqlite"))
  gpa$object_cache <- file.path(gpa$working_directory, paste0(gpa$analysis_name, ".rds"))

  if (is.null(gpa$sqlite_con) || !DBI::dbIsValid(gpa$sqlite_con)) {
    lg$info("Opening SQLite connection to: %s", gpa$sqlite_db)
    gpa$sqlite_con <- DBI::dbConnect(RSQLite::SQLite(), gpa$sqlite_db)
  }

  # new approach: use internal model names for creating output directories at subject level
  # default to <analysis_name>/<l1_model_name>
  # add suffix if using preconvolution approach
  gpa$l1_models$models <- lapply(gpa$l1_models$models, function(mm) {
    mm$outdir <- file.path(gpa$analysis_name, paste0(mm$name, ifelse(gpa$use_preconvolve, "_preconvolve", "")))
    return(mm)
  })

  if (!is.null(gpa$run_number_regex)) {
    if (stringr::str_count(gpa$run_number_regex, stringr::fixed("(")) != 1L) {
      stop(
        "run_number_regex: ", gpa$run_number_regex,
        " must have exactly one opening parenthesis, denoting start of run number capture"
      )
    }

    if (stringr::str_count(gpa$run_number_regex, stringr::fixed(")")) != 1L) {
      stop(
        "run_number_regex: ", gpa$run_number_regex,
        " must have exactly one closing parenthesis, denoting end of run number capture"
      )
    }
  }

  # setup l1 copes, cope names, and contrasts.
  gpa$n_l1_copes <- sapply(gpa$l1_models$models, function(mm) {
    nrow(mm$contrasts)
  }) # number of level 1 copes per model
  gpa$l1_cope_names <- lapply(gpa$l1_models$models, function(mm) {
    rownames(mm$contrasts)
  }) # names of level 1 copes for each model

  #TODO: not currently supported
  #gpa$l1_working_directory <- file.path(gpa$working_directory, gpa$outdir) # temp folder for each analysis variant
  if (is.null(gpa$force_l1_creation)) {
    # whether to overwrite existing level 1 setup files (e.g., .fsf)
    gpa$force_l1_creation <- FALSE
  }

  # ---- PARALLELISM SETUP
  # pipeline_cores: number of cores used in push_pipeline when looping over l1 model variants
  if (is.null(gpa$parallel$pipeline_cores) || gpa$parallel$pipeline_cores == "default") {
    # number of workers to setup at the pipeline level (i.e., over l1 model variants)
    gpa$pipeline_cpus <- length(gpa$l1_models$models)
  }

  # l1_setup_cores defines how many cores to use when looping over subjects within a single l1 model setup
  if (is.null(gpa$parallel$l1_setup_cores) || gpa$parallel$l1_setup_cores == "default") {
    # default to serial execution within a single l1 model variant in setup_lvl1_models
    gpa$parallel$l1_setup_cores <- 1L
  } else {
    checkmate::assert_integerish(gpa$parallel$l1_setup_cores, lower = 1)
  }

  # number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)
  if (is.null(gpa$parallel$slurm)) gpa$parallel$slurm <- list()
  if (is.null(gpa$parallel$torque)) gpa$parallel$torque <- list()
  if (gpa$scheduler == "slurm") {
    if (is.null(gpa$parallel$sched_args)) {
      gpa$parallel$sched_args <- c("-p general")
      lg$info("Using default SLURM scheduler arguments: ")
      lg$info("Argument: %s", gpa$parallel$sched_args)
    }
  } else if (gpa$scheduler == "torque") {
    # gpa$parallel$sched_args <- c("-A mnh5174_c_g_sc_default", "-W group_list=mnh5174_collab")
    if (is.null(gpa$parallel$sched_args)) {
      gpa$parallel$sched_args <- c("-j oe", "-m n")
      lg$info("Using default PBS scheduler arguments: ")
      lg$info("Argument: %s", gpa$parallel$sched_args)
    }
  }

  if (is.null(gpa$parallel$fsl$l2_cores)) gpa$parallel$fsl$l2_cores <- 20
  if (is.null(gpa$parallel$fsl$l1_feat_time)) gpa$parallel$fsl$l1_feat_time <- "6:00:00" # 6 hours
  if (is.null(gpa$parallel$fsl$l1_feat_memgb)) gpa$parallel$fsl$l1_feat_memgb <- "12" # 12 GB by default
  if (is.null(gpa$parallel$fsl$l2_feat_time)) gpa$parallel$fsl$l2_feat_time <- "1:00:00" # 1 hour
  if (is.null(gpa$parallel$fsl$l2_feat_memgb)) gpa$parallel$fsl$l2_feat_memgb <- "12" # 12 GB by default
  if (is.null(gpa$parallel$fsl$l3_feat_time)) gpa$parallel$fsl$l3_feat_time <- "24:00:00" # 24 hours
  if (is.null(gpa$parallel$fsl$l3_feat_memgb)) gpa$parallel$fsl$l3_feat_memgb <- "32" # 32 GB by default

  if (is.null(gpa$parallel$fsl$compute_environment)) {
    lg$info("Using default compute environment for UNC Longleaf")
    gpa$parallel$fsl$compute_environment <- c(
      "module unload fsl", # remove any current fsl module
      "module load fsl/6.0.4" # load latest version (2021)
    )
  }

  #old ICS-ACI settings
  # "source /gpfs/group/mnh5174/default/lab_resources/ni_path.bash",
  # "module unload fsl", #make sure that the ni_path version of FSL is unloaded
  # "#module load \"openblas/0.2.20\" >/dev/null 2>&1",
  # "module load \"fsl/6.0.1\" >/dev/null 2>&1",
  # "module load gsl/2.5", #for dependlab R package to work (some new dependency)

  if (is.null(gpa$parallel$fsl$slurm_l1_array)) {

  }

  # TODO: deprecate this -- should not be required when executing as an R package
  if (is.null(gpa$pipeline_home)) gpa$pipeline_home <- "/proj/mnhallqlab/users/michael/fmri.pipeline"

  if (is.null(gpa$center_l3_predictors)) gpa$center_l3_predictors <- TRUE
  if (is.null(gpa$bad_ids)) gpa$bad_ids <- c()
  if (is.null(gpa$scheduler)) gpa$scheduler <- "slurm" # HPC batch system

  if (is.null(gpa$zthresh)) gpa$zthresh <- 3.09 # 1-tailed p=.001 for z stat
  if (is.null(gpa$clustsize)) gpa$clustsize <- 50 # arbitrary reasonable lower bound on clusters
  if (is.null(gpa$glm_software)) gpa$glm_software <- "fsl" # default to FSL FEAT

  if (is.null(gpa$log_json)) gpa$log_json <- TRUE # whether to write JSON log files
  if (is.null(gpa$log_txt)) gpa$log_txt <- TRUE # whether to write text log files

  if (is.null(gpa$l1_setup_log)) {
    l1_setup_log <- paste0(names(gpa$l1_models$models), "_l1setup") %>% setNames(names(gpa$l1_models$models))
    lg$debug("l1_setup_log is %s", l1_setup_log)
    gpa$l1_setup_log <- l1_setup_log
  }
  if (is.null(gpa$l1_execution_log)) {
    l1_execution_log <- paste0(names(gpa$l1_models$models), "_l1execution") %>% setNames(names(gpa$l1_models$models))
    lg$debug("l1_execution_log is %s", l1_execution_log)
    gpa$l1_execution_log <- l1_execution_log
  }

  if (is.null(gpa$n_expected_runs)) gpa$n_expected_runs <- 1 # assume single run case

  # remove bad ids before running anything further
  if (!is.null(gpa$bad_ids) && length(gpa$bad_ids) > 0L) {
    lg$info("Removing the following IDs from data structure before beginning analysis: ", paste(gpa$bad_ids, collapse=", "))
    gpa$subject_data <- gpa$subject_data %>% filter(!id %in% gpa$bad_ids) # remove bad ids
    gpa$run_data <- gpa$run_data %>% filter(!id %in% gpa$bad_ids) # remove bad ids
    gpa$trial_data <- gpa$trial_data %>% filter(!id %in% gpa$bad_ids) # remove bad ids
  }

  # build design matrix default arguments
  if (is.null(gpa$additional$bdm_args)) {
    gpa$additional$bdm_args <- list(
      baseline_coef_order = 2, center_values = TRUE,
      plot = FALSE, convolve_wi_run = TRUE, output_directory = "run_timing"
    )
  }

  # default settings for feat l1
  if (is.null(gpa$additional$feat_l1_args)) {
    gpa$additional$feat_l1_args <- list(z_thresh = 1.96, prob_thresh = .05, paradigm_hp = 120)
  }

  # identify and validate niftis for each run
  gpa <- lookup_nifti_inputs(gpa)

  # populate confounds
  if (!is.null(gpa$confound_settings$motion_params_file)) {
    checkmate::assert_string(gpa$confound_settings$motion_params_file)
    if ("motion_params" %in% names(gpa$run_data)) {
      message("motion_params column already in run_data. Not using motion_params_file specification.")
    } else {
      gpa$run_data$motion_params <- sapply(seq_len(nrow(gpa$run_data)), function(ii) {
        if (isTRUE(gpa$run_data$run_nifti_present[ii])) {
          normalizePath(file.path(dirname(gpa$run_data$run_nifti[ii]), gpa$confound_settings$motion_params_file))
        } else {
          NA_character_
        }
      })
    }
    gpa$run_data$motion_params_present <- file.exists(get_mr_abspath(gpa$run_data, "motion_params"))
  } else {
    if (!"motion_params" %in% names(gpa$run_data)) {
      gpa$run_data$motion_params <- gpa$run_data$motion_params_present <- NA_character_
    }
  }

  if (!is.null(gpa$confound_settings$confound_input_file)) {
    checkmate::assert_string(gpa$confound_settings$confound_input_file)
    if ("confound_input_file" %in% names(gpa$run_data)) {
      message("confound_input_file column already in run_data. Not using confound_input_file specification.")
    } else {
      gpa$run_data$confound_input_file <- sapply(seq_len(nrow(gpa$run_data)), function(ii) {
        if (isTRUE(gpa$run_data$run_nifti_present[ii])) {
          normalizePath(file.path(dirname(gpa$run_data$run_nifti[ii]), gpa$confound_settings$confound_input_file))
        } else {
          NA_character_
        }
      })
    }
    gpa$run_data$confound_input_file_present <- file.exists(get_mr_abspath(gpa$run_data, "confound_input_file"))
  } else {
    if (!"confound_input_file" %in% names(gpa$run_data)) {
      gpa$run_data$confound_input_file <- gpa$run_data$confound_input_file_present <- NA_character_
    }
  }

  #validate confound settings
  if (is.null(gpa$confound_settings)) {
    lg$info("Using default settings for confounds and exclusions")
    lg$info("Look for confounds in confounds.tsv")
    lg$info("Exclude run if mean(FD) > 0.5 or max(FD) > 6")
    gpa$confound_settings <- list(
      motion_params_file = NULL, 
      motion_params_colnames = NULL,
      confound_input_file = "confounds.tsv",
      l1_confound_regressors = NULL, # column names in motion_params and/or confound_input_file
      exclude_run = "mean(FD) > 0.5 || max(FD) > 6)",
      spike_volume = "FD > 0.9"
    )
  } else {
    checkmate::assert_string(gpa$confound_settings$exclude_run, null.ok=TRUE)
    checkmate::assert_string(gpa$confound_settings$exclude_subject, null.ok = TRUE)
  }

  #figure out all confound columns that will be used in the pipeline
  gpa$confound_settings$run_exclusion_columns <- if (is.null(gpa$confound_settings$exclude_run)) {
    NULL
  } else {
    all.vars(as.formula(paste("~", gpa$confound_settings$exclude_run)))
  }

  #TODO: Should this become 'id_exclusion_columns' and should we support session versus subject exclusion
  # (E.g., in longitudinal analysis)
  gpa$confound_settings$subject_exclusion_columns <- if (is.null(gpa$confound_settings$exclude_subject)) {
    NULL
  } else {
    all.vars(as.formula(paste("~", gpa$confound_settings$exclude_subject)))
  }

  gpa$confound_settings$all_confound_columns <- unique(c(
    gpa$confound_settings$l1_confound_regressors,
    gpa$confound_settings$run_exclusion_columns,
    gpa$confound_settings$subject_exclusion_columns
  ))

  # handle lookup and creation of all confound files


  # numeric row number of each input to aid in tracking
  gpa$run_data$input_number <- seq_len(nrow(gpa$run_data))

  lg$debug("Setting pipeline_finalized to TRUE")
  gpa$pipeline_finalized <- TRUE

  #cache gpa object to file
  saveRDS(gpa, file=gpa$object_cache)

  #save subject, run, and trial data to the database, too
  lg$info("Writing run_data to sqlite db: %s", gpa$sqlite_db)
  DBI::dbWriteTable(conn=gpa$sqlite_con, name="run_data", value=gpa$run_data, overwrite=TRUE)

  lg$info("Writing subject_data to sqlite db: %s", gpa$sqlite_db)
  DBI::dbWriteTable(conn=gpa$sqlite_con, name="subject_data", value=gpa$subject_data, overwrite=TRUE)

  lg$info("Writing trial_data to sqlite db: %s", gpa$sqlite_db)
  DBI::dbWriteTable(conn = gpa$sqlite_con, name = "trial_data", value = gpa$trial_data, overwrite = TRUE)

  # populate confounds in SQLite database
  # TODO: allow external $exclude_run from user, add internal calculated $calc_exclude run

  gpa$run_data$exclude_run <- FALSE
  xx <- sapply(seq_len(nrow(gpa$run_data)), function(ii) {
    browser()
    # this should add rows to the SQLite data for a subject if not yet present, or just return those rows if they exist
    l1_info <- get_l1_confounds(id = gpa$run_data$id[ii], session = gpa$run_data$session[ii], run_number = gpa$run_data$run_number[ii], gpa = gpa, drop_volumes = gpa$drop_volumes)
    insert_df_sqlite(gpa, id = gpa$run_data$id[ii], session = gpa$run_data$session[ii], run_number = gpa$run_data$run_number[ii], data=l1_info$confounds_df, table="test")

    #
    
  })

  #get_l1_confounds <- function(id = NULL, session = NULL, run_number = NULL, gpa, drop_volumes=0L, last_volume=NULL, demean=TRUE) {

  # determine whether to include each run
  mrdf$exclude_run <- sapply(seq_len(nrow(mrdf)), function(rr) {
    ll <- as.list(mrdf[rr, , drop = FALSE]) # rrth row of mrdf
    ll[["gpa"]] <- gpa
    ll[["run_nifti"]] <- NULL
    ex <- do.call(get_l1_confounds, ll)
    return(ex$exclude_run)
  })

  return(gpa)
}
