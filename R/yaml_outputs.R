#' Output GLM configuration to a '.yaml' or 'json' file
#' @param gpa the `glm_pipeline_arguments` object with relevant fields already populated, including
#'   `$l1_models`, `$l2_models`, and/or `$l3_models`.
#' @param file the filename for the exported configuration file
#' @return invisibly return the YAML/JSON syntax written to `file`
#' @details 
#'   This function exports the configuration of level 1, 2, and 3 GLM models in a `glm_pipeline_arguments`
#'   (aka 'gpa') object to a YAML or JSON file. This allows the user to setup their models from this configuration
#'   file, rather than having to go through the building process manually (with prompts).
#' 
#'   Typically, the user will first run `build_l1_models`, `build_l2_models` (for multi-run data), and `build_l3_models`,
#'   then run `export_glm_config` so that the decisions made during the build process are codified in the resulting file. 
#' 
#' @examples 
#' \dontrun{
#'   gpa <- setup_glm_pipeline(analysis_name="flanker_test",
#'     output_directory="/proj/mnhallqlab/no_backup/flanker_gnomes",
#'     n_expected_runs = 2, tr = 2.0
#'     run_data = run_df, subject_data = subject_df, trial_data = trial_df,
#'     l1_models=NULL, l2_models=NULL, l3_models=NULL
#'   )
#' 
#'  # setup level 1 models through menu system
#'  gpa <- build_l1_models(gpa)
#' 
#'  # for multi-run data, setup level 2 models through menu system
#'  gpa <- build_l2_models(gpa)
#' 
#'  # setup level 3 models through menu system
#'  gpa <- build_l3_models(gpa)
#' 
#'  # export model configuration to YAML file
#'  export_glm_config(gpa, file="my_glm_config.yaml")
#' 
#'  # in future, we can build the models from the file, rather than menus
#'  gpa <- build_l1_models(gpa, from_spec_file="my_glm_config.yaml")
#'  gpa <- build_l2_models(gpa, from_spec_file="my_glm_config.yaml")
#'  gpa <- build_l3_models(gpa, from_spec_file="my_glm_config.yaml")
#' }
#' @importFrom yaml as.yaml
#' @importFrom jsonlite toJSON
#' @export
export_glm_config <- function(gpa, file="glm_config.yaml") {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_string(file, null.ok = FALSE)

  have_l1 <- !is.null(gpa$l1_models)
  have_l2 <- isTRUE(gpa$multi_run) && !is.null(gpa$l2_models)
  have_l3 <- !is.null(gpa$l3_models)

  if (!any(c(have_l1, have_l2, have_l3))) {
    warning("In export_glm_config, nothing to export. Have you run build_l1_models, build_l2_models, and build_l3_models?")
    return(character(0))
  }
  
  l1_str <- get_l1_config(gpa)
  l2_str <- get_l2_config(gpa)
  l3_str <- get_l3_config(gpa)

  str_out <- c(l1_str, l2_str, l3_str)

  ext <- file_ext(file)
  if (is.na(ext)) { # missing file extension -- default to YAML
    message("Adding .yaml extension to output file: ", file)
    file <- paste0(file, ".yaml")
    ext <- ".yaml"
  }

  if (ext == ".yaml") {
    str_out <- as.yaml(str_out)
  } else if (ext == ".json") {
    str_out <- toJSON(str_out, pretty = TRUE)
  } else {
    stop("Unable to handle file extension: ", ext)
  }

  writeLines(str_out, file)

  return(invisible(str_out))
}


#' Function to output l1_models into '.yaml' object
#' @param gpa the \code{gpa} object
#' @return a yaml document of l1_models for an object
#' @details 
#'   This function exports the configuration of level 1 (run-level) GLM settings
#'   in a glm_pipeline_arguments (gpa) object to a YAML file. This allows the user to
#'   then setup their models from this configuration file, rather than having to go through
#'   the `build_l1_models` step manually.
#' 
#' @examples
#' \dontrun{
#'   gpa <- setup_glm_pipeline(analysis_name="flanker_test",
#'     output_directory="/proj/mnhallqlab/no_backup/flanker_gnomes",
#'     n_expected_runs = 2, tr = 2.0
#'     run_data = run_df, subject_data = subject_df, trial_data = trial_df,
#'     l1_models=NULL, l2_models=NULL, l3_models=NULL
#'   )
#' 
#'  # setup level 1 models through menu system
#'  gpa <- build_l1_models(gpa)
#' 
#'  # get level 1 model settings as a list
#'  l1_settings <- get_l1_config(gpa)
#' 
#'  # we should typically use export_glm_config to write to file, but we could convert manually
#'  l1_yaml <- as.yaml(l1_settings)
#'  writeLines(l1_yaml, "my_l1_config.yaml")
#' 
#'  # in future, we can populate the models using
#'  gpa <- build_l1_models(gpa, from_spec_file="my_l1_config.yaml")
#' }
#' @keywords internal
get_l1_config <- function(gpa) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  l1m <- gpa$l1_models

  if (is.null(l1m)) {
    warning("Cannot export l1 config because $l1_models is NULL. Run build_l1_models first?")
    return(character(0))
  }

  if (!checkmate::test_class(l1m, "l1_model_set")) {
    stop("$l1_models is not of class l1_model_set")
  }

  ## Variable mappings
  str_out <- list(
    onsets = l1m$onsets, # Onset columns
    durations = l1m$durations, # Duration columns
    isis = l1m$isis, # ISI/ITI columns
    wi_factors = l1m$wi_factors, # within-subject factors
    values = l1m$values # parametric modulators
  )

  # export events
  for (ee in l1m$events) {
    eobj <- list()
    checkmate::assert_string(ee$name)
    checkmate::assert_string(ee$onset)
    if (!checkmate::test_string(ee$duration)) {
      checkmate::assert_number(ee$duration, lower = 0)
    }
    checkmate::assert_string(ee$isi, null.ok = TRUE)
    eobj$onset <- ee$onset
    eobj$duration <- ee$duration
    eobj$isi <- ee$isi
    str_out$events[[ee$name]] <- eobj # populate event list
  }

  # export signals
  for (ss in l1m$signals) {
    sobj <- list()
    sobj$event <- ss$event
    sobj$trial_subset_expression <- ss$trial_subset_expression
    if (!is.null(ss$normalization)) {
      checkmate::assert_subset(ss$normalization, c("none", "evtmax_1", "durmax_1"))
      sobj$normalization <- ss$normalization
    }

    if (!is.null(ss$parametric_modulator)) {
      stopifnot(ss$parametric_modulator %in% l1m$values)
      sobj$parametric_modulator <- ss$parametric_modulator
      sobj$value_type <- "parametric"
      sobj$value_fixed <- NULL
    } else if (!is.null(ss$value_fixed)) {
      checkmate::assert_number(ss$value_fixed)
      sobj$value_fixed <- ss$value_fixed
      if (abs(ss$value_fixed - 1) < 1e-5) {
        sobj$value_type <- "unit"
      } else {
        sobj$value_type <- "number"
      }
    }

    # within-subject formula
    if (!is.null(ss$wi_formula)) {
      sobj$wi_formula <- ss$wi_formula
    }

    # within-subject factors
    if (!is.null(ss$wi_factors)) {
      sobj$wi_factors <- ss$wi_factors
    }

    # other signal flags
    if (!is.null(ss$demean_convolved)) {
      checkmate::assert_flag(ss$demean_convolved)
      sobj$demean_convolved <- ss$demean_convolved
    }

    if (!is.null(ss$add_deriv)) {
      checkmate::assert_flag(ss$add_deriv)
      sobj$add_deriv <- ss$add_deriv
    }

    if (!is.null(ss$beta_series)) {
      checkmate::assert_flag(ss$beta_series)
      sobj$beta_series <- ss$beta_series
    }

    if (!is.null(ss$convmax_1)) {
      checkmate::assert_flag(ss$convmax_1)
      sobj$convmax_1 <- ss$convmax_1
    }

    if (!is.null(ss$ts_multiplier) && !checkmate::test_logical(ss$ts_multiplier)) {
      checkmate::assert_string(ss$ts_multiplier)
      sobj$ts_multiplier <- ss$ts_multiplier
    }

    str_out$signals[[ss$name]] <- sobj
  }

  # L1 models
  for (mm in l1m$models) {
    mobj <- list()
    mobj$signals <- mm$signals
    mobj$contrasts <- list(
      diagonal = mm$contrast_spec$diagonal,
      cell_means = mm$contrast_spec$cell_means,
      cond_means = mm$contrast_spec$cond_means,
      pairwise_diffs = mm$contrast_spec$pairwise_diffs,
      overall_response = mm$contrast_spec$overall_response,
      weights = mm$contrast_spec$weights,
      delete = mm$contrast_spec$delete
    )
    if (!is.null(mm$contrast_spec$simple_slopes) && length(mm$contrast_spec$simple_slopes) > 0L) {
      mobj$contrasts$simple_slopes <- mm$contrast_spec$simple_slopes
    }

    if (!is.null(mm$wi_models)) {
      mobj$contrasts$wi_contrasts <- list()
      for (cc in names(mm$wi_models)) {
        mcopy <- mm$wi_models[[cc]]$contrast_spec
        mcopy$cat_vars <- NULL
        if (!is.null(mcopy$simple_slopes) && length(mm$contrast_spec$simple_slopes) == 0L) mcopy$simple_slopes <- NULL # remove empty
        # mcopy$contrasts <- NULL # do not export constructed matrix
        # mcopy$contrast_list <- NULL
        # mcopy$lmfit <- NULL # should be refit upon import
        # mcopy$wi_formula <- as.character.formula(mcopy$wi_formula)

        # mobj$wi_models[[cc]] <- list(
        #   wi_factors = mm$wi_models[[cc]]$wi_factors
        #   wi_formula = mm$wi_models[[cc]]$wi_factors
        # )
        mobj$contrasts$wi_contrasts[[cc]] <- mcopy
      }
    }

    str_out$l1_models[[mm$name]] <- mobj
  }

  return(invisible(str_out))
}

#' Helper function to output l2_models into '.yaml' object
#' @param gpa the \code{gpa} object
#' @importFrom yaml as.yaml
#' @return a yaml document of l2_models for an object
#' @keywords internal
get_l2_config <- function(gpa) {

  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  # capture whether we are building level 2 (subject) or level 3 (sample) models
  fname <- match.call()[[1]]
  if (fname == "get_l2_config") {
    model_set <- gpa$l2_models
    level <- 2L #model level inside model object
  } else {
    model_set <- gpa$l3_models
    level <- 3L
  }

  if (is.null(model_set)) {
    warning(glue("Cannot export level {level} config because $l{level}_models is NULL. Run build_l{level}_models first?"))
    return(character(0))
  }

  if (!checkmate::test_class(model_set, "hi_model_set")) {
    stop(glue::glue("$l{level}_models is not of class hi_model_set"))
  }

  str_out <- list()

  for (mm in model_set$models) {
    if (!inherits(mm, "hi_model_spec")) {
      warning(glue::glue("In export_l{level}_config, model is not of class 'hi_model_spec'.")) # would be nice if this were more informative, but we don't know the class
      next # skip
    }

    mobj <- list()
    mobj$level <- mm$level
    mobj$model_formula <- as.character(mm$model_formula)
    mobj$num2fac <- mm$num2fac # variables that should be converted from numbers to factors before contrasts are setup
    
    # convert named vectors for covariate transformation and reference levels to lists for YAML to capture the key:value pairs
    if (!is.null(mm$covariate_transform)) mobj$covariate_transform <- as.list(mm$covariate_transform) # settings for covariate transformation
    if (!is.null(mm$reference_level)) mobj$reference_level <- as.list(mm$reference_level) # settings for covariate transformation

    mobj$contrasts <- list(
      diagonal = mm$contrast_spec$diagonal,
      cell_means = mm$contrast_spec$cell_means,
      cond_means = mm$contrast_spec$cond_means,
      pairwise_diffs = mm$contrast_spec$pairwise_diffs,
      overall_response = mm$contrast_spec$overall_response,
      weights = mm$contrast_spec$weights,
      delete = mm$contrast_spec$delete
    )
    if (!is.null(mm$contrast_spec$simple_slopes) && length(mm$contrast_spec$simple_slopes) > 0L) {
      mobj$contrasts$simple_slopes <- mm$contrast_spec$simple_slopes
    }

    if (level == 3L) {
      mobj$fsl_outlier_deweighting <- mm$fsl_outlier_deweighting
    }

    str_out[[mm$name]] <- mobj

  }

  if (level == 2L) {
    str_out <- list(l2_models=str_out) # wrap inside $l2_models top-level element
  } else if (level == 3L) {
    str_out <- list(l3_models=str_out) # wrap inside $l3_models top-level element
  } else {
    stop("Unknown level in get_l2_config")
  }

  return(invisible(str_out))
}

#' Get a list containing configuration settings for l3 models
#' @rdname get_l3_config
#' @name get_l3_config
#' @keywords internal
get_l3_config <- get_l2_config
