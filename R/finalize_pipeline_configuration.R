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
#' @importFrom DBI dbConnect dbIsValid dbWriteTable
#' @importFrom RSQLite SQLite
#' @export
finalize_pipeline_configuration <- function(gpa, refinalize = FALSE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_logical(refinalize, len = 1L)
  lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")
  lg$set_threshold(gpa$lgr_threshold)

  if (isTRUE(gpa$finalize_complete) && isFALSE(refinalize)) {
    lg$debug("In finalize_pipeline_configuration, finalization of gpa already complete. Returning object unchanged.")
    return(gpa)
  }

  if (is.null(gpa$sqlite_con) || !DBI::dbIsValid(gpa$sqlite_con)) {
    lg$info("Opening SQLite connection to: %s", gpa$output_locations$sqlite_db)
    gpa$sqlite_con <- DBI::dbConnect(RSQLite::SQLite(), gpa$output_locations$sqlite_db)
  }

  # final checks on compute environment now that we're running inside the compute environment
  test_compute_environment(gpa, stop_on_fail=TRUE)

  # l1 models must be specified to get started (hard enforcement)
  if (!checkmate::test_class(gpa$l1_models, "l1_model_set")) {
    msg <- "Could not find valid $l1_models specification in gpa. Be sure to run build_l1_models() before proceeding!"
    lg$error(msg)
    stop(msg)
  }

  # l2 models are optional and depend on whether it's a multi-run setup
  if (isTRUE(gpa$multi_run) && !checkmate::test_class(gpa$l2_models, "hi_model_set")) {
    msg <- "Could not find valid $l2_models specification in gpa. We will proceed, but weird things may happen! I suggest running build_l2_models()."
    lg$warn(msg)
    warning(msg)
  }

  # l3 models are optional but almost always should be in place
  if (!checkmate::test_class(gpa$l3_models, "hi_model_set")) {
    msg <- "Could not find valid $l3_models specification in gpa. We will proceed, but weird things may happen! I suggest running build_l3_models()."
    lg$warn(msg)
    warning(msg)
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
  if (is.null(gpa$clustsize)) gpa$clustsize <- 50 # arbitrary reasonable lower bound on cluster size
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
    gpa$subject_data <- gpa$subject_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
    gpa$run_data <- gpa$run_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
    gpa$trial_data <- gpa$trial_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
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

  fsl_defaults <- list(
    force_l1_creation = FALSE, # whether to overwrite existing level 1 setup files (e.g., .fsf)
    failed_l1_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l1_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    force_l2_creation = FALSE, # whether to overwrite existing level 2 setup files (e.g., .fsf)
    failed_l2_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l2_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    force_l3_creation = FALSE, # whether to overwrite existing level 3 setup files (e.g., .fsf)
    failed_l3_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l3_folder_action = "archive", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    replace_l1_nifti_symlink = TRUE # whether to replace filtered_func_data with symlink to the (same) input data
  )

  gpa$glm_settings$fsl <- populate_defaults(gpa$glm_settings$fsl, fsl_defaults)

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
  gpa <- specify_cores(gpa, "pipeline_cores", ifelse(is.null(gpa$l1_models$models), 4, length(gpa$l1_models$models)))

  # l1_setup_cores defines how many cores to use when looping over models in setup_l1_models
  # default to 4 cores in setup_lvl1_models
  gpa <- specify_cores(gpa, "l1_setup_cores", 4)

  # default to 4GB per core for l1 setup
  if (is.null(gpa$parallel$l1_setup_memgb)) gpa$parallel$l1_setup_memgb <- paste0(4 * gpa$parallel$l1_setup_cores, "G")

  # l2_setup_cores defines how many cores to use when looping over models in setup_l2_models
  gpa <- specify_cores(gpa, "l2_setup_cores", 4)

  # finalize_cores defines how many cores to use when finalizing the pipeline before setting up and running models
  gpa <- specify_cores(gpa, "finalize_cores", 4)

  if (is.null(gpa$parallel$slurm)) gpa$parallel$slurm <- list()
  if (is.null(gpa$parallel$torque)) gpa$parallel$torque <- list()

  if (gpa$scheduler == "slurm") {
    if (is.null(gpa$parallel$sched_args)) {
      # Jan 2022: turns out Longleaf doesn't want us to use a partition by default!
      # gpa$parallel$sched_args <- c("-p general")
      # lg$info("Using default SLURM scheduler arguments: ")
      # lg$info("Argument: %s", gpa$parallel$sched_args)
    } else {
      checkmate::assert_character(gpa$parallel$sched_args)
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
  if (is.null(gpa$parallel$finalize_time)) gpa$parallel$finalize_time <- "6:00:00" # 6.0 hours (includes run truncation, which is slow)
  if (is.null(gpa$parallel$l1_setup_time)) gpa$parallel$l1_setup_time <- "4:00:00" # 4.0 hours

  # 14 hours for l2 setup and execution to clear scheduler (all jobs)
  if (is.null(gpa$parallel$l2_setup_run_time)) gpa$parallel$l2_setup_run_time <- "14:00:00"

  # 80 hours for all L3 analyses to clear scheduler (all jobs)
  if (is.null(gpa$parallel$l3_setup_run_time)) gpa$parallel$l3_setup_run_time <- "80:00:00"
  if (is.null(gpa$parallel$compute_environment$global) && isTRUE(grepl("(longleaf|ll\\.unc\\.edu)", gpa$nodename))) {
    lg$info("Using default global compute environment for UNC Longleaf")
    gpa$parallel$compute_environment$global <- "module use /proj/mnhallqlab/sw/modules"
  }

  if (is.null(gpa$parallel$compute_environment$r) && isTRUE(grepl("(longleaf|ll\\.unc\\.edu)", gpa$nodename))) {
    lg$info("Using default R compute environment for UNC Longleaf")
    gpa$parallel$compute_environment$r <- c(
      "module unload r",
      "module load r/4.2.1"
    )
  }

  # fsl_parallel_defaults <- list(
  #   l1_feat_time = "8:00:00", # 8 hours
  # )

  # number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)
  if (is.null(gpa$parallel$fsl$l1_feat_alljobs_time)) gpa$parallel$fsl$l1_feat_alljobs_time <- "72:00:00" # 3 days for all jobs
  if (is.null(gpa$parallel$fsl$l1_feat_time)) gpa$parallel$fsl$l1_feat_time <- "10:00:00" # 10 hours
  if (is.null(gpa$parallel$fsl$l1_feat_memgb)) gpa$parallel$fsl$l1_feat_memgb <- "12" # 12 GB by default
  if (is.null(gpa$parallel$fsl$l1_feat_cpus_per_job)) gpa$parallel$fsl$l1_feat_cpus_per_job <- 8
  if (is.null(gpa$parallel$fsl$l1_feat_runs_per_cpu)) gpa$parallel$fsl$l1_feat_runs_per_cpu <- 2 # 2 runs per cpu in a job
  if (is.null(gpa$parallel$fsl$l2_feat_time)) gpa$parallel$fsl$l2_feat_time <- "2:00:00" # 2 hours
  if (is.null(gpa$parallel$fsl$l2_feat_memgb)) gpa$parallel$fsl$l2_feat_memgb <- "12" # 12 GB by default
  if (is.null(gpa$parallel$fsl$l3_feat_time)) gpa$parallel$fsl$l3_feat_time <- "24:00:00" # 24 hours
  if (is.null(gpa$parallel$fsl$l3_feat_memgb)) gpa$parallel$fsl$l3_feat_memgb <- "32" # 32 GB by default
  if (is.null(gpa$parallel$fsl$l3_feat_cpusperjob)) gpa$parallel$fsl$l3_feat_cpusperjob <- 16 # cpus used to process all slices

  if (is.null(gpa$parallel$compute_environment$fsl) && isTRUE(grepl("(longleaf|ll\\.unc\\.edu)", gpa$nodename))) {
    lg$info("Using default FSL compute environment for UNC Longleaf")
    gpa$parallel$compute_environment$fsl <- c(
      "module unload fsl", # remove any current fsl module
      "module load fsl/6.0.4" # load latest version (2021)
    )
  }

  if (is.null(gpa$parallel$compute_environment$afni) && isTRUE(grepl("(longleaf|ll\\.unc\\.edu)", gpa$nodename))) {
    lg$info("Using default AFNI compute environment for UNC Longleaf")
    gpa$parallel$compute_environment$afni <- c(
      "module unload afni", # remove any current afni module
      "module load afni/23.0.07" # load latest version (2023)
    )
  }

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

  lg$info("Output directory for this analysis will be: %s", gpa$output_directory)

  # see quickstart.Rmd > output settings

  # build out ability to consolidate outputs in one folder, to use specific paths for some outputs, etc.
  # if user specifies gpa$output_directory that matches gpa$analysis_name, don't at this as subfolder

  if (length(unique(gpa$run_data$session)) == 1L) {
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
    # default structure is like: L1m-abspexrew/L2m-modl2_l2c-EV_overall/L3m-int_only/FEAT_l1c-{l1_cope_name}.fsf
    feat_l3_directory = ifelse(isTRUE(gpa$multi_run),
      file.path("{gpa$output_directory}", "feat_l3", "L1m-{l1_model}", "L2m-{l2_model}_l2c-{l2_contrast}", "L3m-{l3_model}"),
      file.path("{gpa$output_directory}", "feat_l3", "L1m_{l1_model}", "L3m-{l3_model}")
    ),
    feat_l3_fsf = "FEAT_l1c-{l1_contrast}.fsf",
    feat_l3_combined_filename = ifelse(isTRUE(gpa$multi_run), # settings for combining FEAT L3 models into a smaller set of AFNI files for visualization
      file.path("{gpa$output_directory}", "feat_l3_combined", "L1m-{l1_model}", "l1c-{l1_cope_name}", "L2m-{l2_model}_L3m-{l3_model}_stats"),
      file.path("{gpa$output_directory}", "feat_l3_combined", "L1m-{l1_model}", "l1c-{l1_cope_name}", "L3m-{l3_model}_stats")
    ),
    feat_l3_combined_briknames = ifelse(isTRUE(gpa$multi_run),
      "l2c-{l2_cope_name}_l3c-{l3_cope_name}",
      "l3c-{l3_cope_name}"
    ),
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

  # validate confound settings
  confound_defaults <- list(
    motion_params_file = NULL,
    motion_params_colnames = NULL,
    confound_input_file = "confounds.tsv",
    l1_confound_regressors = NULL, # column names in motion_params_file and/or confound_input_file
    exclude_run = "mean(framewise_displacement) > 0.5 | max(framewise_displacement) > 6",
    truncate_run = NULL, # example: framewise_displacement > 1 & time > last_onset
    spike_volumes = "framewise_displacement > 0.9",
    na_strings = getOption("datatable.na.strings", "NA") # default na.strings argument for data.table::fread calls
  )

  if (is.null(gpa$confound_settings)) {
    lg$info("Using default settings for confounds and exclusions")
    lg$info("Look for confounds in confounds.tsv")
    lg$info("Exclude run if mean(framewise_displacement) > 0.5 or max(framewise_displacement) > 6")
  } else {
    checkmate::assert_string(gpa$confound_settings$exclude_run, null.ok = TRUE)
    checkmate::assert_string(gpa$confound_settings$exclude_subject, null.ok = TRUE)
  }

  gpa$confound_settings <- populate_defaults(gpa$confound_settings, confound_defaults)

  #####
  # populate confounds
  if ("motion_params_file" %in% names(gpa$run_data)) {
    lg$info("motion_params_file column already in run_data. Not using motion_params_file confounds specification.")
  } else if (!is.null(gpa$confound_settings$motion_params_file)) {
    checkmate::assert_string(gpa$confound_settings$motion_params_file)
    gpa$run_data$motion_params_file <- gpa$confound_settings$motion_params_file # this gets expanded by get_mr_abspath in get_l1_confounds
  }

  if (!"motion_params_file" %in% names(gpa$run_data)) {
    gpa$run_data$motion_params_file <- gpa$run_data$motion_params_present <- NA_character_
  } else {
    gpa$run_data$motion_params_present <- file.exists(get_mr_abspath(gpa$run_data, "motion_params_file"))
  }

  if ("confound_input_file" %in% names(gpa$run_data)) {
    lg$info("confound_input_file column already in run_data. Not using confound_input_file confounds specification.")
  } else if (!is.null(gpa$confound_settings$confound_input_file)) {
    checkmate::assert_string(gpa$confound_settings$confound_input_file)
    gpa$run_data$confound_input_file <- gpa$confound_settings$confound_input_file # this gets expanded by get_mr_abspath in get_l1_confounds
  }

  # determine whether confound input files are present
  if (!"confound_input_file" %in% names(gpa$run_data)) {
    gpa$run_data$confound_input_file <- gpa$run_data$confound_input_file_present <- NA_character_
  } else {
    gpa$run_data$confound_input_file_present <- file.exists(get_mr_abspath(gpa$run_data, "confound_input_file"))
  }

  rhs_to_vars <- function(str) {
    if (is.null(str)) {
      NULL # return NULL if input is NULL
    } else if (checkmate::test_string(str)) {
      if (!grepl("^\\s*~", str)) str <- paste("~", str)
      all.vars(as.formula(paste("~", str)))
    } else if (checkmate::test_formula(str)) {
      all.vars(str)
    } else {
      stop("rhs_to_vars input is not a string or formula")
    }
  }

  # figure out all confound columns that will be used in the pipeline
  gpa$confound_settings$run_exclusion_columns <- rhs_to_vars(gpa$confound_settings$exclude_run)

  # figure out all confound columns that will be used in the pipeline
  gpa$confound_settings$run_truncation_columns <- rhs_to_vars(gpa$confound_settings$truncate_run)

  # TODO: Should this become 'id_exclusion_columns' and should we support session versus subject exclusion
  # (E.g., in longitudinal analysis)
  gpa$confound_settings$subject_exclusion_columns <- rhs_to_vars(gpa$confound_settings$exclude_subject)

  gpa$confound_settings$all_confound_columns <- unique(c(
    gpa$confound_settings$l1_confound_regressors,
    gpa$confound_settings$run_exclusion_columns,
    gpa$confound_settings$run_truncation_columns,
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

  # add onset + offset + isi data to run_data so that get_l1_confounds can handle run truncation appropriately
  gpa <- populate_last_events(gpa, lg)

  # for each run, calculate confounds, exclusions, and truncation
  run_list <- lapply(seq_len(nrow(gpa$run_data)), function(ii) {
    get_l1_confounds(run_df = gpa$run_data[ii, , drop = FALSE], gpa = gpa)
  })

  if (length(unique(sapply(run_list, length))) == 1L) {
    gpa$run_data <- data.table::rbindlist(run_list)
  } else {
    msg <- "Lengths of confound run_df elements have different lengths. Cannot recombine"
    lg$error(msg)
    stop(msg)
  }

  return(gpa)
}
