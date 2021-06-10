#' Interactive function to build an l1 model specification for setup_glm_pipeline
#'
#' @param trial_data a data.frame containing trial-level data for one or more subjects
#' @param l1_model_set optional existing l1_model_set to be modified
#' @param variable_mapping a vector of mappings between columns in \code{trial_data} and internal constructs
#' @param onset_cols an optional character vector of columns in \code{trial_data} that should be
#'   in the set of event onsets
#' @param onset_regex an optional PCRE-compatible regular expression for identifying potential
#'   event onset columns in \code{trial_data}
#' @param duration_regex an optional PCRE-compatible regular expression for identifying potential
#'   event duration columns in \code{trial_data}
#' @param value_cols an optional character vector of columns in \code{trial_data} that should be in the set of signal values
#'
#' @return a \code{l1_model_set} object containing events, signals, and models, compatible with build_design_matrix
#' @author Michael Hallquist
#' @importFrom checkmate assert_data_frame assert_class assert_subset
#' @export
#' 
build_l1_models <- function(gpa=NULL, trial_data=NULL, l1_model_set=NULL,
                           onset_cols=NULL, onset_regex=".*(onset|time).*", 
                           duration_regex=".*duration.*", value_cols=NULL) {

  # Maybe allow glm object to be passed in that would have trial_data and variable_mapping.
  # I guess that would be like "add_l1_model"
  lg <- lgr::get_logger("glm_pipeline/l1_setup")

  checkmate::assert_class(gpa, "glm_pipeline_arguments", null.ok = TRUE)
  if (!is.null(gpa)) {
    lg$info("In build_l1_models, using existing gpa object to build l1 models (ignoring trial_data argument etc.)")
    use_gpa <- TRUE
    trial_data <- gpa$trial_data
    l1_model_set <- gpa$l1_models
  } else {
    lg$info("In build_l1_models, using trial_data passed in, rather than gpa object")
    use_gpa <- FALSE
  }

  checkmate::assert_data_frame(trial_data, null.ok = FALSE)
  checkmate::assert_class(l1_model_set, "l1_model_set", null.ok=TRUE)
  checkmate::assert_subset(onset_cols, names(trial_data)) #make sure that all event columns are in the data frame
  checkmate::assert_string(onset_regex, null.ok=TRUE)
  checkmate::assert_string(duration_regex, null.ok=TRUE)
  checkmate::assert_subset(value_cols, names(trial_data)) #make sure all parametric regressor columns are in the data frame
  checkmate::assert_subset(c("id", "session", "run_number", "trial"), names(trial_data)) # required metadata in trial_data

  lg <- lgr::get_logger("glm_pipeline/build_l1_models")

  if (is.null(l1_model_set)) {
    ## initialize overall l1 design object (holds events, signals, and models)
    l1_model_set <- list(events=NULL, signals=NULL, models=NULL)
    class(l1_model_set) <- c("list", "l1_model_set")
  }

  #onset manager
  get_onsets <- function(trial_data, onset_cols=NULL, onset_regex=NULL) {
    if (!is.null(onset_regex)) {
      detected_onsets <- grep(onset_regex, names(trial_data), value=TRUE, perl=TRUE)
      if (length(detected_onsets) > 0L) {
        cat("Detected the following possible event onset columns:\n\n  ", paste(detected_onsets, collapse=", "), "\n\n")
        res <- menu(c("Yes", "No"), title="Add these columns to possible onsets?")
        if (res == 1) {
          #add these to any that were manually specified/current
          onset_cols <- unique(c(onset_cols, detected_onsets))
        }
      }
    }

    done_onsets <- FALSE
    while (isFALSE(done_onsets)) {
      cat("Current onset columns:\n\n  ", paste(onset_cols, collapse=", "), "\n\n")
      action <- menu(c("Add/modify onset columns", "Delete onset columns", "Done with onset selection"),
        title="Would you like to modify the event onset columns?")

      if (action == 1L) { #Add/modify
        onset_cols <- select.list(names(trial_data), multiple=TRUE, preselect=onset_cols,
          title="Choose all columns denoting event onset times\n(Command/Control-click to select multiple)")
      } else if (action == 2L) { #Delete
        if (length(onset_cols) == 0L) {
          message("No onsets yet. Please add at least one.")
        } else {
          which_del <- menu(onset_cols, title="Which onset column would you like to remove?")
          if (which_del > 0) {
            proceed <- menu(c("Proceed", "Cancel"),
              title=paste0("Are you sure you want to delete ", onset_cols[which_del], "?"))
            if (proceed==1) {
              cat("  Deleting ", onset_cols[which_del], "\n")
              onset_cols <- onset_cols[-which_del]
            } else {
              cat("  Not deleting ", onset_cols[which_del], "\n")
            }
          }
        }
      } else if (action == 3L) { #Done
        done_onsets <- TRUE
        cat("The following columns were chosen as event onset times.",
          "These will be used as possible onset times for each regressor.\n\n", sep="\n")
        cat("  ", paste(onset_cols, collapse=", "), "\n\n")
      }
    }
    return(onset_cols)
  }

  #ask user to setup onsets, build from existing onsets, if available
  if (!is.null(l1_model_set$events)) { onset_cols <- unique(names(l1_model_set$events), onset_cols) }
  onset_cols <- get_onsets(trial_data, onset_cols, onset_regex)

  #basal data frame for each event
  metadata_df <- trial_data %>%
    dplyr::select(id, session, run_number, trial)

  #build a list of data frames, one per event (to be rbind'ed later)
  event_list <- lapply(onset_cols, function(xx) {
    metadata_df %>% bind_cols(trial_data %>%
      select(all_of(xx)) %>%
      setNames("onset") %>%
      mutate(event = xx))
  }) %>% setNames(onset_cols)


  #populate existing events into structure
  if (!is.null(l1_model_set$events)) {
    extant_events <- names(l1_model_set$events)[names(l1_model_set$events) %in% names(event_list)]
    event_list[extant_events] <- l1_model_set$events[extant_events]
  } else {
    extant_events <- c()
  }

  #handle durations
  for (oo in onset_cols) {
    cat("\n----\nSpecify a fixed duration value or column for the event: ", oo, "\n")

    if (oo %in% extant_events) {
      cat("\nCurrent duration value summary, mean [min -- max]:",
        round(mean(event_list[[oo]]$duration, na.rm=TRUE), 2), "[",
        round(min(event_list[[oo]]$duration, na.rm=TRUE), 2), "--",
        round(max(event_list[[oo]]$duration, na.rm=TRUE), 2), "]\n")
      cat("First 6 values:", paste(head(event_list[[oo]]$duration), collapse=", "), "\n\n")
      reselect <- menu(c("Yes", "No (keep durations)"), title="Do you want to respecify the duration?")

      if (reselect == 2) { next } #skip over this column
    }

    choices <- c("Specify fixed duration", names(trial_data))
    oval <- menu(choices=choices)
    if (oval==0) {
      message("Skipping out of duration?") #TODO
    } else if (oval==1) {
      duration <- NULL
      while (!checkmate::test_number(duration, lower=0, upper=5000)) {
        duration <- as.numeric(readline(paste0("Enter the duration value (in seconds) for ", oo, ": ")))
      }
      if (duration > 50) {
        lg$warn("Duration more than 50s specified. Make sure that your durations are in seconds, not milliseconds!")
      }
      event_list[[oo]] <- event_list[[oo]] %>%
        mutate(duration = duration)
    } else {
      event_list[[oo]] <- event_list[[oo]] %>%
        mutate(duration = trial_data[[choices[oval]]])
    }
  }

  #populate into set object
  l1_model_set$events <- event_list

  ####### setup signals
  cat("\nNow, we will build up a set of signals that can be included as regressors in the level 1 model.\n")
  cat("First, select all columns that contain parametric signals to be used for regressors\n")
  cat("If you wish to use a fixed value (e.g., 1.0, a unit-height regressor), this will be entered in the next step\n\n")

  reselect <- 2
  while (reselect == 2L) {
    value_cols <- select.list(c("None (no parametric modulators)", names(trial_data)), multiple=TRUE, preselect=value_cols,
      title="Choose all columns denoting regressor values\n(Command/Control-click to select multiple)")

    if (any(value_cols %in% c(0, 1))) {
      value_cols <- c() #no parametric modulators
      reselect <- 0 #drop loop
    } else {
      cat("The following columns were chosen as regressor values.\nThese will be used as possible values for each signal.\n\n")
      cat("  ", paste(value_cols, collapse=", "), "\n\n")

      reselect <- menu(c("Yes", "No (reselect values)"), title="Are you done selecting all possible regressor values?")
    }
  }

  #helper function to print signal setup
  summarize_signals <- function(sl) {
    if (length(sl) == 0L) { return(invisible(NULL)) }
    lapply(seq_along(sl), function(ii) {
      this <- sl[[ii]]
      cat("--------\nSignal ", ii, "\n\n")
      cat("  Name:", this$name, "\n")
      cat("  Event alignment:", this$event, "\n")
      if (length(this$value) == 1L && is.numeric(this$value[1L])) {
        cat("  Regressor value (constant): ", this$value[1L], "\n")
      } else {
        cat("  Parametric value, mean [min -- max]:",
          round(mean(this$value$value, na.rm=TRUE), 2), "[",
          round(min(this$value$value, na.rm=TRUE), 2), "--",
          round(max(this$value$value, na.rm=TRUE), 2), "]\n")
      }
      cat("  HRF Normalization:", this$normalization, "\n")
      cat("  Add signal derivative:", as.character(this$add_deriv), "\n")
      cat("  Demean convolved signal:", as.character(this$demean_convolved), "\n")
      cat("  Generate beta series:", as.character(this$beta_series), "\n")
      cat("  Multiply convolved regressor against time series:",
        ifelse(this$ts_multipliers == FALSE || is.null(this$ts_multipliers),
          "FALSE", this$ts_multipliers), "\n")
      cat("\n--------\n")
    })
  }

  signal_list <- l1_model_set$signals
  add_more <- 1
  while (add_more != 4) {
    summarize_signals(signal_list)

    add_more <- menu(c("Add signal", "Modify signal", "Delete signal", "Done with signal setup"),
      title="Signal setup menu")

    if (add_more==4) {
      break
    } else if (add_more == 3) {
      if (length(signal_list) == 0L) {
        message("No signals added yet. Please add one first")
      } else {
        which_del <- menu(names(signal_list), title="Which signal would you like to delete?")
        if (which_del > 0) {
          proceed <- menu(c("Proceed", "Cancel"),
            title=paste0("Are you sure you want to delete ", names(signal_list)[which_del], "?"))
          if (proceed==1) {
            cat("  Deleting ", names(signal_list)[which_del], "\n")
            signal_list[[which_del]] <- NULL
          } else {
            cat("  Not deleting ", names(signal_list)[which_del], "\n")
          }
        }
      }
    } else if (add_more %in% c(1, 2)) {
      modify <- FALSE

      #prompt for signal details
      if (add_more == 2) {
        if (length(signal_list) == 0L) {
          message("No signals available to modify")
          ss <- list()
        } else {
          res <- 0L
          while (res == 0L) { res <- menu(names(signal_list), title="Which signal do you want to modify?") }
          ss <- signal_list[[res]]
          signal_list[[res]] <- NULL #clear out old settings
          modify <- TRUE
        }
      } else {
        ss <- list()
      }

      ### ------ name ------
      complete <- FALSE
      while (isFALSE(complete)) {
        prompt <- "Enter the signal name: "
        if (isTRUE(modify)) {
          cat("Current signal name:", ss$name, "\n")
          prompt <- "Enter the signal name (press return to leave unchanged): "
        }

        nm <- readline(prompt)
        if (nm != "") {
          ss$name <- nm
          complete <- TRUE
        } else if (nm == "" && isTRUE(modify)) {
          complete <- TRUE
        }
      }

      ### ------ event alignment ------
      if (isTRUE(modify)) {
        cat("Current signal alignment:", ss$event, "\n")
        res <- menu(c("No", "Yes"), title="Change signal alignment?")
        if (res == 2) { ss$event <- NULL } #clear out event so that it is respecified
      }

      while (is.null(ss$event)) {
        res <- menu(names(event_list), title="With which event is this signal aligned?")
        if (res > 0) { ss$event <- names(event_list)[res] }
      }

      ### ------ value of regressor ------
      if (isTRUE(modify)) {
        cat("Current signal value:",
          ifelse(length(ss$value) == 1L && is.numeric(ss$value[1L]),
            ss$value[1L],
            paste0(round(mean(ss$value$value, na.rm=TRUE), 2), "[",
              round(min(ss$value$value, na.rm=TRUE), 2), "--",
              round(max(ss$value$value, na.rm=TRUE), 2), "]\n")
          ), "\n")
        res <- menu(c("No", "Yes"), title="Change signal value?")
        if (res == 2) { ss$value <- NULL } #clear out event so that it is respecified
      }

      while (is.null(ss$value) || ss$value == 0) {
        regtype <- menu(c("Unit height (1.0)", "Other fixed value (will prompt for details)",
          "Parametric modulator (will prompt for details)"),
          title="What should be the value of regressor (pre-convolution)?")

        if (regtype == 1) {
          ss$value <- 1.0
        } else if (regtype == 2) {
          ss$value <- NULL
          while (!test_number(ss$value)) {
            ss$value <- as.numeric(readline("Enter the regressor value/height (pre-convolution): "))
          }
        } else if (regtype == 3) {
          val <- 0L
          while (val == 0L) {
            val <- menu(value_cols, c("Which value should be used for this signal?"))
            if (val > 0) {
              #TODO: have build_design matrix support a simple value vector, which requires
              #same number of rows as metadata (avoid redundancy)
              ss$value <- metadata_df %>% dplyr::bind_cols(value = trial_data[[ value_cols[val] ]])
            }
          }
        }
      }

      ### ------ hrf normalization ------
      if (isTRUE(modify)) {
        cat("Current HRF normalization:", ss$normalization, "\n")
        res <- menu(c("No", "Yes"), title="Change HRF normalization?")
        if (res == 2) { ss$normalization <- NULL } #clear out so that it is respecified
      }

      while (is.null(ss$normalization)) {
        opt <- c("none", "evtmax_1", "durmax_1")
        res <- menu(c("none",
          "evtmax_1 (aka dmUBLOCK(1); HRF max of 1.0 for each event, regardless of duration)",
          "durmax_1 (aka dmUBLOCK; HRF maxing at 1.0 as events become longer (1.0 around 15 sec)"),
          title="How should the HRF be normalized in convolution?")
        if (res > 0) { ss$normalization <- opt[res] }
      }

      accept_defaults <- 1 # "No"
      prompt_advanced <- FALSE
      if (isTRUE(modify)) {
        cat("Current advanced option settings:\n",
          "  - Temporal derivative:", as.character(ss$add_deriv), "\n",
          "  - Demean convolved signal:", as.character(ss$demean_convolved), "\n",
          "  - Beta series:", as.character(ss$beta_series), "\n",
          "  - Time series multiplier:", as.character(ss$ts_multiplier), "\n\n")

        res <- menu(c("No", "Yes"), title="Change advanced settings?")
        if (res == 2) { #clear out so that it is respecified
          ss$add_deriv <- ss$demean_convolved <- ss$beta_series <- ss$ts_multiplier <- NULL
          prompt_advanced <- TRUE
        }
      } else {
        cat("Advanced options defaults:",
          "  - No temporal derivative",
          "  - Demean convolved signal",
          "  - No beta series",
          "  - No time series multiplier (PPI)\n", sep="\n")
        accept_defaults <- menu(c("Yes", "No"), title="Accept default advanced options for this signal?")
        if (accept_defaults == 2L) { prompt_advanced <- TRUE }
      }

      if (accept_defaults == 1L) {
        cat("Using defaults for signal: ", ss$name, "\n")
        ss$add_deriv <- FALSE
        ss$demean_convolved <- TRUE
        ss$beta_series <- FALSE
      } else if (isTRUE(prompt_advanced)) {
        #derivative
        while (is.null(ss$add_deriv)) {
          res <- menu(c("No", "Yes"), title="Add temporal derivative?")
          if (res == 1L) {
            ss$add_deriv <- FALSE
          } else if (res == 2L) {
            ss$add_deriv <- TRUE
          }
        }

        #demean
        while (is.null(ss$demean_convolved)) {
          res <- menu(c("No", "Yes"), title="Demean signal post-convolution?")
          if (res == 1L) {
            ss$demean_convolved <- FALSE
          } else if (res == 2L) {
            ss$demean_convolved <- TRUE
          }
        }

        #beta series
        while (is.null(ss$beta_series)) {
          res <- menu(c("No", "Yes"),
            title="Generate beta series for this signal (one regressor per trial)?")
          if (res == 1L) {
            ss$beta_series <- FALSE
          } else if (res == 2L) {
            ss$beta_series <- TRUE
          }
        }

        #ts multiplier [not quite there]
        ## while (is.null(ss$beta_series)) {
        ##   res <- menu(c("No", "Yes"),
        ##     title="Generate beta series for this signal (one regressor per trial)?")
        ##   if (res == 1L) { ss$beta_series <- FALSE
        ##   } else if (res == 2L) { ss$beta_series <- TRUE }
        ## }

      }

      signal_list[[ss$name]] <- ss
    }
  }

  #populate back into model set
  l1_model_set$signals <- signal_list

  ############### BUILD MODELS FROM SIGNALS AND EVENTS

  summarize_models <- function(ml) {
    if (length(ml) == 0L) { return(invisible(NULL)) }
    lapply(seq_along(ml), function(ii) {
      this <- ml[[ii]]
      cat("--------\nModel ", ii, "\n\n")
      cat("  Name:", this$name, "\n")
      cat("  Signals:", paste(this$signals, collapse=", "), "\n")
      if (ncol(this$contrasts) < 30) {
        cat("  Contrasts:\n\n")
        print(round(this$contrasts, 3))
        cat("\n--------\n\n")
      } else { cat("More than 30 regressor columns. Not printing contrasts\n") }
    })
  }

  create_new_model <- function(signal_list, to_modify=NULL) {
    checkmate::assert_class(to_modify, "l1_model_spec", null.ok=TRUE)
    if (is.null(to_modify)) {
      mm <- list(level = 1L) #first-level model
      class(mm) <- c("list", "l1_model_spec")
      modify <- FALSE
    } else {
      mm <- to_modify
      modify <- TRUE
    }

    ### ------ model name ------
    if (isTRUE(modify)) {
      cat("Current model name:", mm$name, "\n")
      res <- menu(c("No", "Yes"), title="Change model name?")
      if (res == 2) { mm$name <- NULL } #clear out so that it is respecified
    }

    while (is.null(mm$name) || mm$name == "") {
      res <- trimws(readline("Enter the model name: "))
      if (res != "") {
        res <- make.names(res)
        if (res %in% names(model_list)) {
          cat("\nModel name:", res, "already exists. Names must be unique.\n")
          cat("Current models:", paste(names(model_list), collapse=", "), "\n")
        } else {
          mm$name <- res
        }
      }
    }

    if (isTRUE(modify)) {
      cat("Current model signals:", paste(names(mm$signals), collapse=", "), "\n")
      res <- menu(c("No", "Yes"), title="Change model signals (and contrasts)?")
      if (res == 2) { #clear out so that it is respecified
        mm$signals <- mm$regressors <- mm$contrasts <- NULL
      }
    }

    #signals
    summarize_signals(signal_list) #print summary

    while (is.null(mm$signals)) {
      signals <- select.list(names(signal_list), multiple=TRUE, preselect=mm$signals,
        title="Choose all signals to include in this model\n(Command/Control-click to select multiple)")

      if (length(signals) == 0L) {
        proceed <- menu(c("Yes", "No"), title="Nothing entered. Do you want to cancel model setup?")
        if (proceed == 1L) {
          return(invisible(NULL)) #return nothing from function
        }
      } else {
        mm$signals <- signals
      }
    }

    #look up what the regressors will be for this.
    if (is.null(mm$regressors)) {
      mm$regressors <- unlist(lapply(mm$signals, function(nn) {
        if (isTRUE(signal_list[[nn]]$add_deriv)) {
          return(c(nn, paste0("d_", nn))) # regressor and temporal derivative
        } else if (isTRUE(signal_list[[nn]]$beta_series)) {
          if (is.data.frame(signal_list[[nn]]$value)) {
            #TODO: this approach is imperfect if there are jumps in trials for a subject
            #this assumes that all subjects have all trials
            trials <- sort(unique(signal_list[[nn]]$value$trial)) #vector of trials for parametric signal
          } else {
            #trials will be in corresponding event in case value is a scalar
            trials <- sort(unique(event_list[[ signal_list[[nn]]$event ]]$trial))
          }

          return(paste(nn, sprintf("%03d", trials), sep="_t")) # signal_t001 etc.
        } else {
          return(nn) #just the signal name itself (no deriv, no bs)
        }
      }))
    }

    #contrast editor
    mm <- specify_contrasts(mm)

    return(mm)
  }

  model_list <- l1_model_set$models
  add_more <- 1
  while (add_more != 4) {
    summarize_models(model_list)

    add_more <- menu(c("Add model", "Modify model", "Delete model", "Done with l1 model setup"),
      title="Level 1 model setup menu")

    if (add_more == 1L) { #add
      mm <- create_new_model(signal_list)
      if (mm$name %in% names(l1_model_set)) { warning("A model with the same name exists: ", mm$name, ". Overwriting it.") }
      model_list[[mm$name]] <- mm #add to set
    } else if (add_more == 2L) { #modify
      if (is.null(model_list)) {
        message("No models available to modify. Add at least one model first")
      } else {
        res <- 0L
        while (res == 0L) { res <- menu(names(model_list), title="Which model do you want to modify?") }
        model_list[[res]] <- create_new_model(signal_list, to_modify=model_list[[res]])
      }
    } else if (add_more == 3L) { #delete
      which_del <- menu(names(model_list), title="Which model would you like to delete?")
      if (which_del > 0) {
        proceed <- menu(c("Proceed", "Cancel"),
          title=paste0("Are you sure you want to delete ", names(model_list)[which_del], "?"))
        if (proceed==1) {
          cat("  Deleting ", names(model_list)[which_del], "\n")
          model_list[[which_del]] <- NULL
        } else {
          cat("  Not deleting ", names(model_list)[which_del], "\n")
        }
      }
    }

  }

  l1_model_set$models <- model_list

  if (isTRUE(use_gpa)) {
    gpa$l1_models <- l1_model_set
    return(gpa)
  } else {
    return(l1_model_set)
  }
  
}
