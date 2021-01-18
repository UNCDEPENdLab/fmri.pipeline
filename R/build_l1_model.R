#' helper function to build an l1 model specification interactively
#'
#' @param trial_data A data.frame of trial-level events
build_l1_model <- function(trial_data, variable_mapping=c(id="id", run="run", trial="trial", run_trial="trial", mr_dir="mr_dir"),
                           onset_cols=NULL, onset_regex=".*(onset|time).*", duration_regex=".*duration.*", value_cols=NULL) {


  
  #maybe allow glm object to be passed in that would have trial_data and variable_mapping. I guess that would be like "add_l1_model"
  
  checkmate::assert_data_frame(trial_data) #yeah, move toward allowing the broader model specification object here
  checkmate::assert_subset(onset_cols, names(trial_data)) #make sure that all event columns are in the data frame
  checkmate::assert_string(onset_regex, null.ok=TRUE)
  checkmate::assert_string(duration_regex, null.ok=TRUE)
  checkmate::assert_subset(value_cols, names(trial_data)) #make sure all parametric regressor columns are in the data frame
  
  
  possible_cols <- names(trial_data)
  possible_cols <- possible_cols[!names(possible_cols) %in% variable_mapping]
  
  if (is.null(onset_cols)) {
    if (!is.null(onset_regex)) {
      onset_cols <- grep(onset_regex, names(trial_data), value=TRUE, perl=TRUE)
      if (length(onset_cols) > 0L) {
        cat("Detected the following possible event onset columns:\n\n  ", paste(onset_cols, collapse=", "), "\n\n")
        done <- menu(c("Yes (add columns)", "No (done with onsets)"), title="Would you like to select additional event onset columns?")
        reselect <- ifelse(done==1, 2, 1) #a value of 2 will let the user tack on columns
      }
    } else {
      reselect <- 2
      onset_cols <- c()
    }

    while(reselect==2) {
      onset_cols <- select.list(names(trial_data), multiple=TRUE, preselect=onset_cols,
        title="Choose all columns denoting event onset times\n(Command/Control-click to select multiple)")
      cat("The following columns were chosen as event onset times.\n  These will be used as possible onset times for each regressor.\n\n")
      cat("  ", paste(onset_cols, collapse=", "), "\n\n")
      
      reselect <- menu(c("Yes", "No (reselect events)"), title="Are you done selecting all event onset times?")
    }
  }

  #basal data frame for each event
  metadata_df <- trial_data %>% dplyr::select(!!variable_mapping[c("id", "run", "run_trial")]) %>%
    setNames(c("id", "run", "trial"))
  
  #build a list of data frames, one per event (to be rbind'ed later)
  event_list <- lapply(onset_cols, function(xx) {
    metadata_df %>% bind_cols(trial_data %>% select(all_of(xx)) %>% setNames("onset") %>% mutate(event=xx))
  }) %>% setNames(onset_cols)
      
  #handle durations
  for (oo in onset_cols) {
    cat("Specify a duration value or column for the event onset: ", oo, "\n")
    choices <- c("Specify fixed duration", names(trial_data))
    oval <- menu(choices=choices)
    if (oval==1) {
      duration <- NULL
      while(!test_number(duration, lower=0, upper=5000)) {
        duration <- as.numeric(readline(paste0("Enter the duration value (in seconds) for ", oo, ": ")))
      }
      if (duration > 50) { warning("Duration more than 50s specified. Make sure that your durations are in seconds, not milliseconds!") }
      event_list[[oo]] <- event_list[[oo]] %>% mutate(duration=duration)
    } else {
      event_list[[oo]] <- event_list[[oo]] %>% mutate(duration=trial_data[[choices[oval] ]])
    }
  }
  
  ## setup signals
  cat("Now, we will build up a set of signals to be included as regressor in the level 1 model\n")
  cat("First, select all columns that contain parametric signals to be used for regressors\n")
  cat("If you wish to use a fixed value (e.g., 1.0, a unit-height regressor), this will be entered in the next step\n")


  reselect <- 2
  while(reselect==2) {
    value_cols <- select.list(c("None (no parametric modulators)", names(trial_data)), multiple=TRUE, preselect=value_cols,
      title="Choose all columns denoting regressor values\n(Command/Control-click to select multiple)")

    if (any(value_cols %in% c(0, 1))) {
      value_cols <- c() #no parametric modulators
    } else {     
      cat("The following columns were chosen as regressor values.\n  These will be used as possible values for each signal.\n\n")
      cat("  ", paste(value_cols, collapse=", "), "\n\n")
      
      reselect <- menu(c("Yes", "No (reselect values)"), title="Are you done selecting all possible regressor values?")
    }
  }

  summarize_signals <- function(sl) {
    if (length(sl) == 0L) { return(invisible(NULL)) }
    lapply(1:length(sl), function(ii) {
      this <- sl[[ii]]
      cat("--------\nSignal ", ii, "\n\n")
      cat("  Name: ", this$name, "\n")
      cat("  Event alignment: ", this$event, "\n")
      if (length(this$value) == 1L && is.numeric(this$value[1L])) {
        cat("  Regressor value (constant): ", this$value[1L], "\n")
      } else {
        cat("  Parametric value, mean [min -- max]: ",
          round(mean(this$value$value, na.rm=TRUE), 2), " [",
          round(min(this$value$value, na.rm=TRUE), 2), " -- ",
          round(max(this$value$value, na.rm=TRUE), 2), "]\n")
      }
      cat("  HRF Normalization: ", this$normalization, "\n")
      cat("  Add signal derivative: ", as.character(this$add_deriv), "\n")
      cat("  Demean convolved signal: ", as.character(this$demean_convolved), "\n")
      cat("  Generate beta series: ", as.character(this$beta_series), "\n")
      cat("  Multiply convolved regressor against time series: ",
        ifelse(this$ts_multipliers == FALSE || is.null(this$ts_multipliers),
          "FALSE", this$ts_multipliers), "\n")
      cat("\n--------\n")      
    })
  }
  
  signal_list <- list()
  add_more <- 1
  while (add_more != 3) {
    summarize_signals(signal_list)
    
    add_more <- menu(c("Add signal", "Delete signal", "Done with signal setup"),
      title="Signal setup menu")

    if (add_more==3) {      
      break
    } else if (add_more == 2) {
      which_del <- menu(names(signal_list), title="Which signal would you like to delete?")
      if (which_del > 0) {
        proceed <- menu(c("Proceed", "Cancel"),
          title="Are you sure you want to delete ", names(signal_list)[which_del])
        if (proceed==1) {
          cat("  Deleting ", names(signal_list)[which_del], "\n")
          signal_list[[which_del]] <- NULL
        } else {
          cat("  Not deleting ", names(signal_list)[which_del], "\n")
        }
      }
    } else if (add_more == 1) {
      #prompt for signal details
      ss <- list()

      #name
      while(is.null(ss$name) || ss$name == "") {
        ss$name <- readline("Enter the signal name: ")
      }

      #alignment
      while (is.null(ss$event) || ss$event == 0) {
        ss$event <- menu(names(event_list), title="With which event is this signal aligned?")
      }

      #value
      while (is.null(ss$value) || ss$value == 0) {
        regtype <- menu(c("Unit height (1.0)", "Other fixed value (will prompt for details)",
          "Parametric modulator (will prompt for details)"),
          title="What should be the value of regressor (pre-convolution)?")

        if (regtype == 1) {
          ss$value <- 1.0
        } else if (regtype == 2) {
          ss$value <- -1
          while(!test_number(ss$value)) {
            ss$value <- as.numeric(readline("Enter the regressor value/height (pre-convolution): "))
          }          
        } else if (regtype == 3) {
          val <- 0L
          while(valfound == 0L) {
            val <- menu(value_cols, c("Which value should be used for this signal?"))
            if (val > 0) {
              #TODO: have build_design matrix support a simple value vector, which requires
              #same number of rows as metadata (avoid redundancy)
              ss$value <- metadata_df %>% bind_cols(trial_data[[ value_cols[val] ]])
            }
          }
        }
      }
      
      #hrf normalization
      while (is.null(ss$normalization) || ss$normalization == 0) {
        ss[["normalization"]] <- menu(c("none",
          "evtmax_1 (aka dmUBLOCK(1); HRF max of 1.0 for each event, regardless of duration)",
          "durmax_1 (aka dmUBLOCK; HRF maxing at 1.0 as events become longer (1.0 around 15 sec)"),
          title="How should the HRF be normalized in convolution?")
      }

      #derivative
      while (is.null(ss$add_deriv) || ss$add_deriv == 0L) {
        ss$add_deriv <- menu(c("No", "Yes"), title="Add temporal derivative?")
        if (ss$add_deriv == 1L) { ss$add_deriv <- FALSE
        } else if (ss$add_deriv == 2L) { ss$add_deriv <- TRUE }
      }

      #demean
      while (is.null(ss$demean_convolved) || ss$demean_convolved == 0L) {
        ss$demean_convolved <- menu(c("No", "Yes"), title="Demean signal post-convolution?")
        if (ss$demean_convolved == 1L) { ss$demean_convolved <- FALSE
        } else if (ss$demean_convolved == 2L) { ss$demean_convolved <- TRUE }
      }      

      #beta series
      while (is.null(ss$beta_series) || ss$beta_series == 0L) {
        ss$beta_series <- menu(c("No", "Yes"),
          title="Generate beta series for this signal (one regressor per trial)?")
        if (ss$beta_series == 1L) { ss$beta_series <- FALSE
        } else if (ss$beta_series == 2L) { ss$beta_series <- TRUE }
      }

      #ts multiplier [not quite there]
      ## while (is.null(ss$beta_series) || ss$beta_series == 0L) {
      ##   ss$beta_series <- menu(c("No", "Yes"),
      ##     title="Generate beta series for this signal (one regressor per trial)?")
      ##   if (ss$beta_series == 1L) { ss$beta_series <- FALSE
      ##   } else if (ss$beta_series == 2L) { ss$beta_series <- TRUE }
      ## }      
      
    }
  }
  
  return(event_list)
}
