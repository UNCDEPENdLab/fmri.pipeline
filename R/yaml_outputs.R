#' Helper function to output l1_models into '.yaml' object
#' @param gpa the \code{gpa} object
#' @importFrom yaml as.yaml
#' @return a yaml document of l1_models for an object
#' @export
export_l1_config <- function(gpa, file = "l1_config.yaml") {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  l1m <- gpa$l1_models

  if (is.null(l1m)) {
    warning("Cannot export l1 config because $l1_models is NULL. Run build_l1_models first?")
  }

  if (!checkmate::test_class(l1m, "l1_model_set")) {
    stop("$l1_models is not of class l1_model_set")
  }

  ## Variable mappings
  end_yaml <- list(
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
    end_yaml$events[[ee$name]] <- eobj # populate event list
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
      sobj$wi_formula <- as.character.formula(ss$wi_formula)
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

    if (!is.null(ss$ts_multiplier)) {
      checkmate::assert_string(ss$ts_multiplier)
      sobj$ts_multiplier <- ss$ts_multiplier
    }

    end_yaml$signals[[ss$name]] <- sobj
  }

  # L1 models
  for (mm in l1m$models) {
    mobj <- list()
    mobj$signals <- mm$signals

    mobj$diagonal <- mm$contrast_spec$diagonal
    # diag(mm$contrasts) <- ifelse(mobj$diagonal==TRUE, 0, diag(mm$contrasts))
    mobj$cell_means <- mm$contrast_spec$cell_means
    mobj$overall_response <- mm$contrast_spec$overall_response
    if (!is.null(mm$contrast_spec$simple_slopes) && length(mm$contrast_spec$simple_slopes) > 0L) {
      mobj$simple_slopes <- mm$contrast_spec$simple_slopes
    }
    mobj$weights <- mm$contrast_spec$weights
    mobj$delete <- mm$contrast_spec$delete
    
    if (!is.null(mm$wi_models)) {
      mobj$wi_contrasts <- list()
      for (cc in names(mm$wi_models)) {
        mcopy <- mm$wi_models[[cc]]$contrast_spec
        mcopy$cat_vars <- NULL
        if (!is.null(mcopy$simple_slopes) && length(mm$contrast_spec$simple_slopes) == 0L) mcopy$simple_slopes <- NULL # remove empty
        #mcopy$contrasts <- NULL # do not export constructed matrix
        #mcopy$contrast_list <- NULL
        #mcopy$lmfit <- NULL # should be refit upon import
        #mcopy$wi_formula <- as.character.formula(mcopy$wi_formula)

        # mobj$wi_models[[cc]] <- list(
        #   wi_factors = mm$wi_models[[cc]]$wi_factors
        #   wi_formula = mm$wi_models[[cc]]$wi_factors
        # )
        mobj$wi_contrasts[[cc]] <- mcopy
      }
    }


    # cobj <- list()
    # for (cc in seq_along(colnames(mm$contrasts))) {
    #   name_c <- NULL
    #   name_c <- colnames(mm$contrasts)[cc]
    #   cobj$row <- rownames(mm$contrasts)[which(mm$contrasts[cc, ] != 0)]
    #   cobj$contrast <- mm$contrasts[cc, which(mm$contrasts[cc, ] != 0)]
    #   mobj$contrasts[[name_c]] <- cobj
    # }
    end_yaml$l1_models[[mm$name]] <- mobj
  }

  yaml_str <- as.yaml(end_yaml)
  
  if (!is.null(file)) {
    if (file_ext(file) != ".yaml") {
      message("Adding .yaml extension to output file: ", file)
      file <- paste0(file, ".yaml")
    }
    writeLines(yaml_str, file)
  }

  return(invisible(yaml_str))
}

l2_yaml <- function(gpa) {
  gpa2 <- gpa$l2_models
  ## Model Specifications
  end_yaml <- list()
  mobj <- list()
    for (mm in gpa2$models) {
    name_m <- mm$name
    mobj$diagonal <- mm$contrast_spec$diagonal
    #diag(mm$contrasts) <- ifelse(mobj$diagonal==TRUE, 0, diag(mm$contrasts))
    mobj$cell_means <- mm$contrast_spec$cell_means
    mobj$overall_response <- mm$contrast_spec$overall_response
    mobj$weights <- mm$contrast_spec$weights
    mobj$regressors <- mm$regressors
    mobj$signals <- mm$signals
    cobj <- list()
    for (cc in seq_along(colnames(mm$contrasts))) {
      name_c <- NULL
      name_c <- colnames(mm$contrasts)[cc]
      cobj$row <- rownames(mm$contrasts)[which(mm$contrasts[cc,]!=0)]
      cobj$contrast <- mm$contrasts[cc,which(mm$contrasts[cc,]!=0)]
      mobj$contrasts[[name_c]] <- cobj
    }
    end_yaml$l2_models[[name_m]] <- mobj
    mobj$contrasts <- NULL
  }
  yaml_str <- as.yaml(end_yaml)
  yaml_choice <- menu(c(
          "Console Output",
          "File Output",
          "Both",
          "Exit"), title= "How would you like to receive the YAML file?")
          if (yaml_choice == 1) {
            return(cat(yaml_str, "\nGoodbye.\n"))
          } else if (yaml_choice == 2) {
            var <- readline(prompt = "Enter Output File Name: ")
            var = as.character(var)
            var2 <- paste(c(var, "yaml"), collapse=".")
            writeLines(yaml_str, var2)
            return(cat("\nFile should be seen as", var2, "\nGoodbye.\n"))
          } else if (yaml_choice == 3) {
            var <- readline(prompt = "Enter Output File Name: ")
            var = as.character(var)
            var2 <- paste(c(var, "yaml"), collapse=".")
            writeLines(yaml_str, var2)
            return(cat(yaml_str, "\nFile should be seen as", var2, "\nGoodbye.\n"))
          } else (
            return(cat("\nGoodbye.\n"))
          )
}

