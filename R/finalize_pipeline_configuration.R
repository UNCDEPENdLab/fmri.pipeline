#' This is a small helper function to validate the glm_model_arguments list structure.
#' It adds a few details such as the output directory to make it less burdensome for to setup a pipeline
#' N.B. gpa is a shorthand abbreviation for glm_model_arguments, to save typing
#'
#' @param gpa A \code{glm_pipeline_arguments} object setup by \code{setup_glm_pipeline}
#' @param refinalize A logical indicating whether to force checks and finalize steps on an object that
#'   was previously finalized.
#' @importFrom stringr str_count fixed
#' @importFrom magrittr %>%
#' @importFrom lgr get_logger
#' @export
finalize_pipeline_configuration <- function(gpa, refinalize = FALSE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_logical(refinalize, len = 1L)
  lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")

  if (isTRUE(gpa$finalize_complete) && isFALSE(refinalize)) {
    lg$debug("In finalize_pipeline_configuration, finalization of gpa already complete. Returning object unchanged.")
    return(gpa)
  }

  if (is.null(gpa$sqlite_con) || !DBI::dbIsValid(gpa$sqlite_con)) {
    lg$info("Opening SQLite connection to: %s", gpa$output_locations$sqlite_db)
    gpa$sqlite_con <- DBI::dbConnect(RSQLite::SQLite(), gpa$output_locations$sqlite_db)
  }

  # new approach: use internal model names for creating output directories at subject level
  # default to <analysis_name>/<l1_model>
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
  gpa$l1_cope_names <- lapply(gpa$l1_models$models, function(mm) {
    rownames(mm$contrasts)
  }) # names of level 1 copes for each model


  if (is.null(gpa$center_l3_predictors)) gpa$center_l3_predictors <- TRUE
  if (is.null(gpa$bad_ids)) gpa$bad_ids <- c()
  if (is.null(gpa$scheduler)) gpa$scheduler <- "slurm" # HPC batch system
  gpa$scheduler <- tolower(gpa$scheduler)
  checkmate::assert_subset(gpa$scheduler, c("slurm", "sbatch", "torque", "qsub", "local", "sh"))

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
    lg$info("Removing the following IDs from data structure before beginning analysis: %s", paste(gpa$bad_ids, collapse = ", "))
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

  #####
  # handle GLM settings and defaults
  if (is.null(gpa$glm_settings) || gpa$glm_settings[1L] == "default") {
    lg$info("Using default settings for GLM implementation")
    gpa$glm_settings <- list(
      fsl = list(),
      afni = list(),
      spm = list()
    )
  }

  if (is.null(gpa$glm_settings$fsl$force_l1_creation)) {
    # whether to overwrite existing level 1 setup files (e.g., .fsf)
    gpa$glm_settings$fsl$force_l1_creation <- FALSE
  }

  if (is.null(gpa$glm_settings$fsl$force_l2_creation)) {
    # whether to overwrite existing level 2 setup files (e.g., .fsf)
    gpa$glm_settings$fsl$force_l2_creation <- FALSE
  }

  if (is.null(gpa$glm_settings$fsl$force_l3_creation)) {
    # whether to overwrite existing level 3 setup files (e.g., .fsf)
    gpa$glm_settings$fsl$force_l3_creation <- FALSE
  }

  # process confound settings
  gpa <- finalize_confound_settings(gpa, lg)

  # populate subject exclusions
  gpa <- calculate_subject_exclusions(gpa)

  # cache gpa object to file
  res <- tryCatch(saveRDS(gpa, file = gpa$output_locations$object_cache), error = function(e) {
    lg$error("Could not save gpa object to file: %s", gpa$output_locations$object_cache)
    return(NULL)
  })

  # save subject, run, and trial data to the database, too
  lg$info("Writing run_data to sqlite db: %s", gpa$output_locations$sqlite_db)
  DBI::dbWriteTable(conn = gpa$sqlite_con, name = "run_data", value = gpa$run_data, overwrite = TRUE)

  lg$info("Writing subject_data to sqlite db: %s", gpa$output_locations$sqlite_db)
  DBI::dbWriteTable(conn = gpa$sqlite_con, name = "subject_data", value = gpa$subject_data, overwrite = TRUE)

  lg$info("Writing trial_data to sqlite db: %s", gpa$output_locations$sqlite_db)
  DBI::dbWriteTable(conn = gpa$sqlite_con, name = "trial_data", value = gpa$trial_data, overwrite = TRUE)

  lg$debug("Setting finalize_complete to TRUE")
  gpa$finalize_complete <- TRUE

  return(gpa)
}

setup_parallel_settings <- function(gpa, lg = NULL) {
  checkmate::assert_class(lg, "Logger")

  # ---- PARALLELISM SETUP
  specify_cores <- function(gpa, field_name, default=1L) {
    # l1_setup_cores defines how many cores to use when looping over subjects within setup_l1_models
    if (is.null(gpa$parallel[[field_name]]) || gpa$parallel[[field_name]] == "default") {
      gpa$parallel[[field_name]] <- default
    }

    checkmate::assert_integerish(gpa$parallel[[field_name]], lower = 1)
    return(gpa)
  }

  # pipeline_cores: number of cores used in push_pipeline when looping over l1 model variants
  gpa <- specify_cores(gpa, "pipeline_cores", ifelse(is.null(gpa$l1_models), 4, length(gpa$l1_models$models)))

  # l1_setup_cores defines how many cores to use when looping over models in setup_l1_models
  # default to 4 cores in setup_lvl1_models
  gpa <- specify_cores(gpa, "l1_setup_cores", 4)
  if (is.null(gpa$parallel$l1_setup_memgb)) gpa$parallel$fsl$l1_setup_memgb <- "16G"

  # l2_setup_cores defines how many cores to use when looping over models in setup_l2_models
  gpa <- specify_cores(gpa, "l2_setup_cores", 4)

  # finalize_cores defines how many cores to use when finalizing the pipeline before setting up and running models
  gpa <- specify_cores("finalize_cores", 4)

  if (is.null(gpa$parallel$slurm)) gpa$parallel$slurm <- list()
  if (is.null(gpa$parallel$torque)) gpa$parallel$torque <- list()
  if (is.null(gpa$parallel$batch_code)) {

  }
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

  # time for finalize_pipeline_configuration in run_glm_pipeline
  if (is.null(gpa$parallel$finalize_time)) gpa$parallel$finalize_time <- "1:30:00" # 1.5 hours
  if (is.null(gpa$parallel$l1_setup_time)) gpa$parallel$l1_setup_time <- "1:30:00" # 1.5 hours
  if (is.null(gpa$parallel$l2_setup_time)) gpa$parallel$l2_setup_time <- "1:30:00" # 1.5 hours
  if (is.null(gpa$parallel$compute_environment)) {
    lg$info("Using default R compute environment for UNC Longleaf")
    gpa$parallel$compute_environment <- c(
      "module use /proj/mnhallqlab/sw/modules",
      "module load r/4.0.3_depend"
    )
  }

  # number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)
  if (is.null(gpa$parallel$fsl$l2_cores)) gpa$parallel$fsl$l2_cores <- 20
  if (is.null(gpa$parallel$fsl$l1_feat_time)) gpa$parallel$fsl$l1_feat_time <- "6:00:00" # 6 hours
  if (is.null(gpa$parallel$fsl$l1_feat_memgb)) gpa$parallel$fsl$l1_feat_memgb <- "12" # 12 GB by default
  if (is.null(gpa$parallel$fsl$l2_feat_time)) gpa$parallel$fsl$l2_feat_time <- "1:00:00" # 1 hour
  if (is.null(gpa$parallel$fsl$l2_feat_memgb)) gpa$parallel$fsl$l2_feat_memgb <- "12" # 12 GB by default
  if (is.null(gpa$parallel$fsl$l3_feat_time)) gpa$parallel$fsl$l3_feat_time <- "24:00:00" # 24 hours
  if (is.null(gpa$parallel$fsl$l3_feat_memgb)) gpa$parallel$fsl$l3_feat_memgb <- "32" # 32 GB by default
  if (is.null(gpa$parallel$fsl$l3_feat_cpusperjob)) gpa$parallel$fsl$l3_feat_cpusperjob <- 16 # cpus used to process all slices

  if (is.null(gpa$parallel$fsl$compute_environment)) {
    lg$info("Using default FSL compute environment for UNC Longleaf")
    gpa$parallel$fsl$compute_environment <- c(
      "module unload fsl", # remove any current fsl module
      "module load fsl/6.0.4" # load latest version (2021)
    )
  }

  # old ICS-ACI settings
  # "source /gpfs/group/mnh5174/default/lab_resources/ni_path.bash",
  # "module unload fsl", #make sure that the ni_path version of FSL is unloaded
  # "#module load \"openblas/0.2.20\" >/dev/null 2>&1",
  # "module load \"fsl/6.0.1\" >/dev/null 2>&1",
  # "module load gsl/2.5", #for dependlab R package to work (some new dependency)

  if (is.null(gpa$parallel$fsl$slurm_l1_array)) {

  }

  return(gpa)
}

# helper function for settings up $output_directory and $output_locations
setup_output_locations <- function(gpa, lg = NULL) {
  checkmate::assert_class(lg, "Logger")

  # sort out file locations
  if (is.null(gpa$output_directory) || gpa$output_directory == "default") {
    gpa$output_directory <- file.path(getwd(), gpa$analysis_name)
  }

  if (!basename(gpa$output_directory) == gpa$analysis_name) {
    lg$info("Appending analysis_name %s to output_directory %s", gpa$analysis_name, gpa$output_directory)
  }

  lg$info("Output directory for this analysis will be: %s", gpa$output_directory)

  # see quickstart.Rmd > output settings

  # build out ability to consolidate outputs in one folder, to use specific paths for some outputs, etc.
  # if user specifies gpa$output_directory that matches gpa$analysis_name, don't at this as subfolder

  if (length(unique(gpa$session)) == 1L) {
    feat_sub_directory <- file.path("{gpa$output_directory}", "feat_l1", "sub-{id}")
  } else {
    feat_sub_directory <- file.path("{gpa$output_directory}", "feat_l1", "sub-{id}", "ses-{session}")
  }

  output_defaults <- list(
    # default to BIDS-style consolidated output
    consolidated = TRUE,
    feat_sub_directory = feat_sub_directory,
    feat_ses_directory = feat_sub_directory, # no difference in defaults
    feat_l1_directory = file.path(feat_sub_directory, "{l1_model}"),
    feat_l2_directory = feat_sub_directory,
    feat_l3_directory = file.path(gpa$output_directory, "feat_l3", "{l1_contrast}", "{l1_model}", "{l2_contrast}"),
    scheduler_scripts = file.path(gpa$output_directory, "scheduler_scripts"),
    sqlite_db = file.path(gpa$output_directory, paste0(gpa$analysis_name, ".sqlite")),
    object_cache = file.path(gpa$output_directory, paste0(gpa$analysis_name, ".rds")),
    setup_l1_log_txt = file.path(gpa$output_directory, "setup_l1_models.txt"),
    setup_l1_log_json = file.path(gpa$output_directory, "setup_l1_models.json"),
    setup_l2_log_txt = file.path(gpa$output_directory, "setup_l2_models.txt"),
    setup_l2_log_json = file.path(gpa$output_directory, "setup_l2_models.json"),
    setup_l3_log_txt = file.path(gpa$output_directory, "setup_l3_models.txt"),
    setup_l3_log_json = file.path(gpa$output_directory, "setup_l3_models.json")
  )

  if (checkmate::test_string(gpa$output_locations) && gpa$output_locations[1L] == "default") {
    gpa$output_locations <- output_defaults
  }

  miss_fields <- setdiff(names(output_defaults), names(gpa$output_locations))
  if (length(miss_fields) > 0L) {
    for (mm in miss_fields) {
      lg$info("Populating missing $output_locations field: %s with default: %s", mm, output_defaults[[mm]])
      gpa$output_locations[[mm]] <- output_defaults[[mm]]
    }
  }

  return(gpa)
}

#' Helper function to populate confound information for pipeline files
#' 
#' @param gpa a glm_pipeline_arguments object for population
#' @param lg a Logger object for logging results of confound processing
#' @keywords internal
#' @importFrom parallel mclapply
finalize_confound_settings <- function(gpa, lg) {
  checkmate::assert_class(lg, "Logger")
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  #####
  # populate confounds
  if (!is.null(gpa$confound_settings$motion_params_file)) {
    checkmate::assert_string(gpa$confound_settings$motion_params_file)
    if ("motion_params_file" %in% names(gpa$run_data)) {
      lg$info("motion_params_file column already in run_data. Not using $confound_settings$motion_params_file specification.")
    } else {
      gpa$run_data$motion_params_file <- sapply(seq_len(nrow(gpa$run_data)), function(ii) {
        if (isTRUE(gpa$run_data$run_nifti_present[ii])) {
          normalizePath(file.path(dirname(gpa$run_data$run_nifti[ii]), gpa$confound_settings$motion_params_file))
        } else {
          NA_character_
        }
      })
    }
    gpa$run_data$motion_params_present <- file.exists(get_mr_abspath(gpa$run_data, "motion_params_file"))
  } else {
    if (!"motion_params_file" %in% names(gpa$run_data)) {
      gpa$run_data$motion_params_file <- gpa$run_data$motion_params_present <- NA_character_
    }
  }

  if ("confound_input_file" %in% names(gpa$run_data)) {
    lg$info("confound_input_file column already in run_data. Not using confound_input_file specification.")
  } else if (!is.null(gpa$confound_settings$confound_input_file)) {
    checkmate::assert_string(gpa$confound_settings$confound_input_file)
    gpa$run_data$confound_input_file <- sapply(seq_len(nrow(gpa$run_data)), function(ii) {
      if (isTRUE(gpa$run_data$run_nifti_present[ii])) {
        normalizePath(file.path(dirname(gpa$run_data$run_nifti[ii]), gpa$confound_settings$confound_input_file))
      } else {
        NA_character_
      }
    })
  }

  # determine whether confound input files are present
  gpa$run_data$confound_input_file_present <- if ("confound_input_file" %in% names(gpa$run_data)) {
    file.exists(get_mr_abspath(gpa$run_data, "confound_input_file"))
  } else {
    NA_character_
  }

  # validate confound settings
  if (is.null(gpa$confound_settings)) {
    lg$info("Using default settings for confounds and exclusions")
    lg$info("Look for confounds in confounds.tsv")
    lg$info("Exclude run if mean(FD) > 0.5 or max(FD) > 6")
    gpa$confound_settings <- list(
      motion_params_file = NULL,
      motion_params_colnames = NULL,
      confound_input_file = "confounds.tsv",
      l1_confound_regressors = NULL, # column names in motion_params_file and/or confound_input_file
      exclude_run = "mean(FD) > 0.5 || max(FD) > 6)",
      spike_volumes = "FD > 0.9"
    )
  } else {
    checkmate::assert_string(gpa$confound_settings$exclude_run, null.ok = TRUE)
    checkmate::assert_string(gpa$confound_settings$exclude_subject, null.ok = TRUE)
  }

  # default na.strings argument for data.table::fread calls to motion parameters and confounds
  if (is.null(gpa$confound_settings$na_strings)) {
    gpa$confound_settings$na_strings <- getOption("datatable.na.strings", "NA")
  }

  # figure out all confound columns that will be used in the pipeline
  gpa$confound_settings$run_exclusion_columns <- if (is.null(gpa$confound_settings$exclude_run)) {
    NULL
  } else {
    all.vars(as.formula(paste("~", gpa$confound_settings$exclude_run)))
  }

  # TODO: Should this become 'id_exclusion_columns' and should we support session versus subject exclusion
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

  # populate confounds in SQLite database and calculate run exclusions
  # TODO: allow external $exclude_run from user, add internal calculated $calc_exclude run

  # confound_info <- parallel::mclapply(seq_len(nrow(gpa$run_data)), function(ii) {
  #   # this should add rows to the SQLite data for a subject if not yet present, or just return those rows if they exist
  #   l1_info <- get_l1_confounds(
  #     id = gpa$run_data$id[ii], session = gpa$run_data$session[ii], run_number = gpa$run_data$run_number[ii],
  #     gpa = gpa, drop_volumes = gpa$drop_volumes
  #   )[c("l1_confound_file", "exclude_run")]
  #   return(l1_info)
  # }, mc.cores = gpa$parallel$finalize_cores)

  confound_info <- lapply(seq_len(nrow(gpa$run_data)), function(ii) {
    # this should add rows to the SQLite data for a subject if not yet present, or just return those rows if they exist
    l1_info <- get_l1_confounds(
      id = gpa$run_data$id[ii], session = gpa$run_data$session[ii], run_number = gpa$run_data$run_number[ii],
      gpa = gpa, drop_volumes = gpa$drop_volumes
    )[c("l1_confound_file", "exclude_run")]
    return(l1_info)
  })

  gpa$run_data$exclude_run <- sapply(confound_info, "[[", "exclude_run")
  gpa$run_data$l1_confound_file <- sapply(confound_info, "[[", "l1_confound_file")

  return(gpa)
}
