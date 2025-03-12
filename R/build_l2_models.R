#' Interactive function to build an l2 model specification for setup_glm_pipeline
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing an analysis pipeline to which $l2_models
#'   should be added. If $l2_models is already present, these will be amended.
#' @param from_spec_file optional YAML or JSON file containing settings to populated into l2 models
#' @param regressor_cols an optional character vector of columns in \code{data} that should be considered
#'   as possible regressors
#'
#' @return a \code{hi_model_set} object containing a list of higher-level regression models
#' @author Michael Hallquist
#' @importFrom checkmate assert_data_frame assert_class assert_subset
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom dplyr select
#' @export
#'
build_l2_models <- function(gpa, from_spec_file = NULL, regressor_cols = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  # capture whether we are building level 2 (subject) or level 3 (sample) models
  fname <- match.call()[[1]]
  if (fname == "build_l2_models") {
    menu_desc <- "second-level (subject)"
    data <- gpa$run_data # L2
    model_set <- gpa$l2_models
    level <- 2L # model level inside model object
    out_field <- "l2_models" # where to put results in gpa
  } else {
    menu_desc <- "third-level (sample)"
    data <- gpa$subject_data # L3
    model_set <- gpa$l3_models
    level <- 3L
    out_field <- "l3_models" # where to put results in gpa
  }

  lg <- lgr::get_logger(paste0("glm_pipeline/l", level, "_setup"))
  lg$set_threshold(gpa$lgr_threshold)

  # allow deferred model specification upstream
  if (!is.null(model_set) && model_set[1L] == "prompt") model_set <- NULL

  checkmate::assert_class(model_set, "hi_model_set", null.ok = TRUE)
  checkmate::assert_subset(regressor_cols, names(data)) # make sure all parametric regressor columns are in the data frame

  possible_cols <- names(data)
  possible_cols <- possible_cols[!possible_cols %in% "mr_dir"]
  data <- data %>% dplyr::select(all_of(possible_cols))

  if (!is.null(regressor_cols)) { # select only columns requested by user
    data <- data %>% dplyr::select(regressor_cols)
  }

  if (!is.null(from_spec_file)) {
    checkmate::assert_file_exists(from_spec_file)
    lg$info("In build_l%d_models, populated fields from: %s", level, from_spec_file)
    ext <- file_ext(from_spec_file, withdot = FALSE)
    if (ext == "yaml") {
      spec_list <- yaml::read_yaml(from_spec_file)
    } else if (ext == "json") {
      spec_list <- jsonlite::read_json(from_spec_file, simplifyVector = TRUE)
    } else {
      msg <- sprintf("Cannot understand how to input spec file: %s", from_spec_file)
      lg$error(msg)
      stop(msg)
    }

    prompt_models <- FALSE

    # at present, l2 and l3 model specs don't necessarily have subsections -- each top-level section can be a model
    # given this, standardize the format so that if $l2_models or $l3_models are top-level sections, then
    # spec_list becomes only the relevant part for building l2/l3 models

    # first, for degenerate files, add $l2_models or $l3_models as top-level
    if (is.null(spec_list[[out_field]])) {
      spec_list <- list(spec_list) %>% setNames(out_field)
    }
    # if (level == 2L && is.null(spec_list$l2_models)) {
    #   spec_list <- list(l2_models = spec_list)
    # } else if (level == 3L && is.null(spec_list$l3_models)) {
    #   spec_list <- list(l3_models = spec_list)
    # }

    # next, populate model names
    spec_list <- propagate_spec_names(spec_list) # ensure that $name field is populated for each model

    # finally, subset spec_list down to just hte relevant part
    spec_list <- spec_list[[out_field]]

    nm <- names(spec_list)
    if (any(duplicated(nm))) {
      msg <- sprintf(
        "Duplicate %s model names in %s. Cannot proceed. Resolve duplicates: %s",
        out_field, from_spec_file, paste(unique(nm[duplicated[nm]]), collapse = ", ")
      )
      lg$error(msg)
      stop(msg)
    }

    # if this changes in future, amend propagate_spec_names to work for l2/l3
    # populate the $name field to make it easier for downstream functions to work
    # for (mm in seq_along(spec_list)) {
    #  spec_list[[mm]]$name <- names(spec_list)[mm]
    # }

    ## initialize overall higher-level model set object
    model_set <- list(models = NULL)
    class(model_set) <- c("list", "hi_model_set")

    for (mm in seq_along(spec_list)) {
      # convert covariate_transform and reference_level back to named vectors for object consistency
      # these are converted to lists during the export so that YAML/JSON preserves the key:value pairs
      spec_list[[mm]]$covariate_transform <- unlist(spec_list[[mm]]$covariate_transform)
      spec_list[[mm]]$reference_level <- unlist(spec_list[[mm]]$reference_level)
      model_set$models[[names(spec_list)[mm]]] <- create_new_hi_model(data, level = level, cur_model_names = NULL, spec_list = spec_list[[mm]], lg = lg)
    }
  } else {
    prompt_models <- TRUE
  }

  if (is.null(model_set)) {
    ## initialize overall higher-level model set object
    model_set <- list(models = NULL)
    class(model_set) <- c("list", "hi_model_set")
  }

  # regressor manager
  get_regressors <- function(data, regressor_cols = NULL) {
    done_regressors <- FALSE
    while (isFALSE(done_regressors)) {
      cat("Current regressors for this model:\n\n  ", paste(regressor_cols, collapse = ", "), "\n\n")
      action <- menu(c("Add/modify regressors", "Delete regressors", "Done with regressor selection"),
        title = "Would you like to modify the model regressors?"
      )

      if (action == 1L) { # Add/modify
        regressor_cols <- select.list(names(data),
          multiple = TRUE, preselect = regressor_cols,
          title = "Choose all model regressors\n(Command/Control-click to select multiple)"
        )
      } else if (action == 2L) { # Delete
        if (length(regressor_cols) == 0L) {
          message("No regressors yet. Please add at least one.")
        } else {
          which_del <- menu(regressor_cols, title = "Which regressor would you like to remove?")
          if (which_del > 0) {
            proceed <- menu(c("Proceed", "Cancel"),
              title = paste0("Are you sure you want to delete ", regressor_cols[which_del], "?")
            )
            if (proceed == 1) {
              cat("  Deleting ", regressor_cols[which_del], "\n")
              regressor_cols <- regressor_cols[-which_del]
            } else {
              cat("  Not deleting ", regressor_cols[which_del], "\n")
            }
          }
        }
      } else if (action == 3L) { # Done
        if (length(regressor_cols) == 0L) {
          message("No regressors yet. Please add at least one.")
        } else {
          done_regressors <- TRUE
          cat("The following columns were chosen as model regressors.\n\n")
          cat("  ", paste(regressor_cols, collapse = ", "), "\n\n")
        }
      }
    }
    return(regressor_cols)
  }

  ### -- BUILD MODELS ---

  summarize_l2_models <- function(ml) {
    if (length(ml) == 0L) {
      return(invisible(NULL))
    }
    lapply(seq_along(ml), function(ii) {
      this <- ml[[ii]]
      cat("--------\nModel ", ii, "\n\n")
      cat("  Name:", this$name, "\n")
      cat("  Regressors:", paste(this$model_variables, collapse = ", "), "\n")
      if (ncol(this$contrasts) < 30) {
        cat("  Contrasts:\n\n")
        print(round(this$contrasts, 3))
        cat("\n--------\n\n")
      } else {
        cat("More than 30 regressor columns. Not printing contrasts\n")
      }
    })
  }

  model_list <- model_set$models
  add_more <- ifelse(isTRUE(prompt_models), 1, 4)
  while (add_more != 4) {
    summarize_l2_models(model_list)

    add_more <- menu(c("Add model", "Modify model", "Delete model", paste("Done with", menu_desc, "model setup")),
      title = paste(sub("^(\\w{1})", "\\U\\1", menu_desc, perl = TRUE), "model setup menu")
    )

    if (add_more == 1L) { # add
      mobj <- create_new_hi_model(data, level = level, cur_model_names = names(model_list), lg = lg)
      if (mobj$name %in% names(model_set)) {
        warning("An model with the same name exists: ", mobj$name, ". Overwriting it.")
      }
      model_list[[mobj$name]] <- mobj # add to set
    } else if (add_more == 2L) { # modify
      if (is.null(model_list)) {
        message("No models available to modify. Add at least one model first")
      } else {
        res <- 0L
        while (res == 0L) {
          res <- menu(names(model_list), title = "Which model do you want to modify?")
        }
        model_list[[res]] <- create_new_hi_model(data, to_modify = model_list[[res]], cur_model_names = names(model_list), level = level, lg = lg)
      }
    } else if (add_more == 3L) { # delete
      which_del <- menu(names(model_list), title = "Which model would you like to delete?")
      if (which_del > 0) {
        proceed <- menu(c("Proceed", "Cancel"),
          title = paste0("Are you sure you want to delete ", names(model_list)[which_del], "?")
        )
        if (proceed == 1) {
          cat("  Deleting ", names(model_list)[which_del], "\n")
          model_list[[which_del]] <- NULL
        } else {
          cat("  Not deleting ", names(model_list)[which_del], "\n")
        }
      }
    }
  }

  model_set$models <- model_list
  model_set$n_contrasts <- sapply(model_list, function(mobj) { ncol(mobj$contrasts) })

  if (fname == "build_l2_models") {
    gpa$l2_models <- model_set
  } else if (fname == "build_l3_models") {
    gpa$l3_models <- model_set
  }

  return(gpa)
}

#' Interactive function to build an l3 model specification for setup_glm_pipeline
#' @rdname build_l3_models
#' @name build_l3_models
#' @export
build_l3_models <- build_l2_models


#' Internal function to walk through l2/l3 model setup, including populated from a specification (YAML) file
#' @param data a data.frame containing data that can be included in the model
#' @param to_modify an existing `hi_model_spec` object to be modified
#' @param level the model level: 2 or 3
#' @param cur_model_names a character vector of the names of all current models in the model list
#' @param spec_list a list containing model specifications read in from a YAML/JSON file
#' @param lg the current logger
#' @details if `spec_list` is passed in, the resulting model object will be populated from the fields
#'   of the specification list, rather than prompting the user for input
#' @return a `hi_model_spec` object containing the L2/L3 model to be added to gpa$l2_models or gpa$l3_models
#' @importFrom checkmate assert_integerish assert_list assert_class test_data_table
#' @keywords internal
create_new_hi_model <- function(data, to_modify = NULL, level = NULL, cur_model_names = NULL, spec_list = NULL, lg = NULL) {
    if (checkmate::test_data_table(data)) {
      data <- as.data.frame(data) # make subsetting syntax in this function consistent with standard data.frame conventions
    }

    checkmate::assert_class(to_modify, "hi_model_spec", null.ok = TRUE)
    if (is.null(to_modify)) {
      mobj <- list(level = level)
      class(mobj) <- c("list", "hi_model_spec")
      modify <- FALSE
    } else {
      mobj <- to_modify
      modify <- TRUE
    }

    checkmate::assert_integerish(level, lower = 2, upper = 3, len = 1L)
    checkmate::assert_character(cur_model_names, null.ok = TRUE)
    checkmate::assert_list(spec_list, null.ok = TRUE)

    if (level == 2) {
      id_cols <- c("id", "session", "run_number")
    } else {
      id_cols <- c("id", "session")
    }

    ### ------ model name ------
    if (isTRUE(modify)) {
      cat("Current model name:", mobj$name, "\n")
      res <- menu(c("No", "Yes"), title = "Change model name?")
      if (res == 2) mobj$name <- NULL # clear out so that it is respecified
    }

    if (!is.null(spec_list$name)) { # populate from specification, if requested
      mobj$name <- spec_list$name
    }

    while (is.null(mobj$name) || mobj$name == "") {
      res <- trimws(readline("Enter the model name: "))
      if (res != "") {
        res <- make.names(res)
        if (res %in% cur_model_names) {
          cat("\nModel name:", res, "already exists. Names must be unique.\n")
          cat("Current models:", paste(cur_model_names, collapse = ", "), "\n")
        } else {
          mobj$name <- res
        }
      }
    }

    if (isTRUE(modify)) {
      cat("Current model regressors:", paste(names(mobj$model_variables), collapse = ", "), "\n")
      res <- menu(c("No", "Yes"), title = "Respecify model regressors (and contrasts)?")
      if (res == 2) { # clear out so that it is respecified
        cat("Okay, resetting model regressors and contrasts\n")
        mobj$model_variables <- mobj$contrasts <- NULL
      }
    }

    # Setup lm formula: let user decide LM formula approach versus walkthrough
    if (!is.null(spec_list$model_formula)) {
      model_approach <- 1L
      model_formula <- as.formula(spec_list$model_formula) # maybe need better error trapping here
    } else {
      model_approach <- 0L
      model_formula <- NULL
    }
    
    while (model_approach == 0L) {
      model_approach <- menu(c("Specify formula", "Model builder"),
        title = "Do you want to specify the model formula or walk through the model builder?"
      )
    }

    if (model_approach == 1L) {
      # formula approach
      while (is.null(model_formula)) {
        cat(c(
          "\nSpecify the right-hand side of the model you wish to fit.",
          "Use variable names in the dataset provided to this function.",
          "Note that this syntax follows standard R model syntax. See ?lm for details.",
          "Example: ~ emotion * wmload + run_number\n",
          "For an intercept-only model, specify: ~1\n",
          "Available column names: \n"
        ), sep = "\n")
        cat(strwrap(paste(names(data), collapse = ", "), 70, exdent = 5), sep = "\n")

        res <- trimws(readline("Enter the model formula: "))
        # always trim any LHS specification
        res <- sub("^[^~]*", "", res, perl = TRUE)
        model_formula <- tryCatch(as.formula(res), error = function(e) {
          print(e)
          cat("Problem converting your syntax to formula. Try again\n")
          return(NULL)
        })
      }

      mobj$model_variables <- all.vars(model_formula)

      # if not an intercept-only model, ensure that all variables in formula are present in the data frame
      if (!identical(model_formula, ~1)) {
        checkmate::assert_subset(mobj$model_variables, names(data), empty.ok = FALSE)
      }

    } else if (model_approach == 2L) {
      # walkthrough approach (only support additive model for now)
      mobj$model_variables <- get_regressors(data, regressor_cols = mobj$model_variables)
      model_formula <- as.formula(paste("~", paste(mobj$model_variables, collapse = " + ")))
    }

    # convert integer-like variables to integer
    for (vv in mobj$model_variables) {
      # convert integer-like numbers to integers for possible conversion below
      if (is.numeric(data[[vv]]) && checkmate::test_integerish(data[[vv]])) {
        data[[vv]] <- as.integer(data[[vv]])
      }

      if (is.character(data[[vv]])) {
        cat("Converting character variable to factor for model setup: ", vv, "\n", sep = "")
        data[[vv]] <- as.factor(data[[vv]])
      }

      if (is.integer(data[[vv]]) && length(unique(data[[vv]]) < 20)) {
        if (vv %in% spec_list$num2fac) {
          conv <- 1L
        } else {
          cat(
            "\n",
            vv, " appears to be an integer with < 20 levels.\n",
            "If you convert it to a factor, R will generate dummy code variables in the model.\n",
            "If you leave it as a number, R will include the linear trend in the model.\n\n",
            sep = ""
          )
          conv <- menu(c("Convert to factor", "Leave as number"),
            title = paste("Do you want to convert", vv, "to a factor?")
          )
        }
        if (conv == 1L) {
          cat("Converting", vv, "into a factor.\n")
          data[[vv]] <- as.factor(data[[vv]])
          mobj$num2fac <- c(mobj$num2fac, vv)
        }
      }
    }

    # need to build a model LM-style
    if (!is.null(spec_list)) {
      mobj$covariate_transform <- spec_list$covariate_transform
      mobj$reference_level <- spec_list$reference_level
    } else if (length(mobj$model_variables) > 0L) {
      cat("Summary of variables included in model:\n\n")
      print(summary(data[, mobj$model_variables, drop = FALSE]))
      cat("\n")

      # handle mean centering and reference levels
      cont_vars <- sapply(data[, mobj$model_variables, drop=FALSE], class) %in% c("integer", "numeric")
      if (any(cont_vars)) {
        cat(
          "We will now ask you to indicate whether to transform any of the continuous covariates in the model.",
          "The most common transformation, which is often a good default, is to mean-center the covariate.",
          "Other options are to subtract the minimum (so that zero now represents the lowest covariate value),",
          "subtract the maximum (so that zero represents the highest covariate value),",
          "or to z-score/standardized the covariate to a mean of zero and standard deviation of one.",
          sep = "\n  "
        )
        cont_vars <- mobj$model_variables[cont_vars == TRUE]
        mobj$covariate_transform <- c() # reset any previous centering settings
        for (cc in cont_vars) {
          cat("\n")
          resp <- menu(
            c("Mean center", "Standardize/z-score", "Subtract minimum", "Subtract maximum", "No transformation"),
            title = paste0("Transform ", cc, "?")
          )
          if (resp == 1L) {
            cat(glue("\nMean centering {cc}. Mean centering is recomputed if runs or subjects are dropped in L2 or L3 analyses.\n", .trim = FALSE))
            mobj$covariate_transform[cc] <- "mean"
          } else if (resp == 2L) {
            cat(glue("\nStandardizing {cc}. Standardization is recomputed if runs or subjects are dropped in L2 or L3 analyses.\n", .trim = FALSE))
            mobj$covariate_transform[cc] <- "zscore"
          } else if (resp == 3L) {
            cat(glue("\nSubtracting minimum value of {cc}. This is recomputed if runs or subjects are dropped in L2 or L3 analyses.\n", .trim = FALSE))
            mobj$covariate_transform[cc] <- "min"
          } else if (resp == 4L) {
            cat(glue("\nSubtracting maximum value of {cc}. This is recomputed if runs or subjects are dropped in L2 or L3 analyses.\n", .trim = FALSE))
            mobj$covariate_transform[cc] <- "max"
          }
        }
      }

      # handle reference levels for factors
      cat_vars <- sapply(data[, mobj$model_variables, drop=FALSE], class) %in% c("factor")
      if (any(cat_vars)) {
        cat_vars <- mobj$model_variables[cat_vars == TRUE]
        cat(
          "We will now ask you to choose the reference level for dummy coding factors in the model.",
          "This affects how contrasts are coded by emmeans and it controls how you interpret the Intercept maps",
          "at levels 2 and 3, if you look at these. In general, this choice should not make a big difference,",
          "but we default to the factor level that is most common in the dataset.", sep="\n  "
        )
        for (cc in cat_vars) {
          # default to cell with largest N
          default_lev <- names(which.max(table(data[[cc]])))
          levs <- levels(data[[cc]])
          levs <- levs[c(which(levs == default_lev), which(levs != default_lev))]
          labs <- levs
          labs[1L] <- paste0("Recommended: ", labs[1])
          which_lev <- menu(labs, title = glue("Choose the reference level for the {cc} factor."))
          if (!exists("which_lev") || which_lev == 0L) {
            which_lev <- 1L # choose default on cancel or ctrl+c
          }
          mobj$reference_level[cc] <- levs[which_lev]
        }

      }
    }

    # fit linear model and populate model object
    mobj <- mobj_fit_lm(mobj, as.character(model_formula), data, id_cols, lg = lg)

    # walk through contrast generation for this model
    mobj <- specify_contrasts(mobj, spec_list = spec_list)

    # ask about outlier deweighting
    if (level == 3L) {
      if (!is.null(spec_list)) {
        mobj$fsl_outlier_deweighting <- spec_list$fsl_outlier_deweighting
      } else {
        cat(
          "",
          "Do you want to turn on outlier deweighting in FSL level 3 analysis? This is available for any mixed effects analysis,",
          "most commonly FLAME1 or FLAME1+2. In general, this is a good idea, but we have seen it break down for some legitimate inputs,",
          "generating many complaints about excessive outliers being detected. If you get those messages, you should probably disable this.",
          "",
          sep = "\n  "
        )

        res <- menu(c("No", "Yes"), title = "Turn on outlier deweighting?")
        if (res == 2L) {
          mobj$fsl_outlier_deweighting <- TRUE
        } else {
          mobj$fsl_outlier_deweighting <- FALSE
        }
      }
    }

    return(mobj)
  }