#' internal function to populate l1_model_set onsets from spec specification
fields_from_spec <- function(l1_model_set, slist, trial_data, field_names=NULL) {
  checkmate::assert_class(l1_model_set, "l1_model_set", null.ok = FALSE)
  checkmate::assert_list(slist)
  checkmate::assert_data_frame(trial_data)
  checkmate::assert_subset(field_names, c("onsets", "durations", "isis", "values", "wi_factors"))

  field_mapping <- c(
    "onsets" = "onset",
    "durations" = "duration",
    "isis" = "ISI/ITI",
    "values" = "parametric value",
    "wi_factors" = "within-subject factor"
  )

  for (fn in field_names) {
    field_vals <- slist[[fn]]
    if (is.null(field_vals)) {
      next # nothing in yaml
    } else {
      checkmate::assert_subset(field_vals, names(trial_data))
    }

    l1_model_set <- bl1_get_cols(l1_model_set, trial_data,
      field_name = fn, field_desc = field_mapping[fn],
      select_cols = field_vals, prompt_input = FALSE
    )
  }

  return(l1_model_set)
}

events_from_spec <- function(l1_model_set, slist, trial_data) {
  checkmate::assert_class(l1_model_set, "l1_model_set", null.ok = FALSE)
  checkmate::assert_list(slist)
  checkmate::assert_data_frame(trial_data)

  if (is.null(slist$events)) return(l1_model_set)

  for (ee in slist$events) {
    eobj <- list()
    checkmate::assert_string(ee$name)
    checkmate::assert_string(ee$onset)
    if (!checkmate::test_string(ee$duration)) {
      checkmate::assert_number(ee$duration, lower=0)
    }
    checkmate::assert_string(ee$isi, null.ok = TRUE)
    eobj$name <- ee$name
    eobj$onset <- ee$onset
    eobj$duration <- ee$duration
    eobj$isi <- ee$isi
    eobj <- populate_event_data(eobj, trial_data)
    class(eobj) <- c("list", "l1_model_set_events") # set 'l1_model_set_events' class
    l1_model_set$events[[eobj$name]] <- eobj # this will overwrite existing specification
  }

  return(l1_model_set)
}

signals_from_spec <- function(l1_model_set, slist, trial_data, lg=NULL) {
  checkmate::assert_class(l1_model_set, "l1_model_set", null.ok = FALSE)
  checkmate::assert_list(slist)
  checkmate::assert_data_frame(trial_data)

  if (is.null(slist$signals)) return(l1_model_set)

  signal_list <- list() # holding tank for populated signals
  for (ss in slist$signals) {
    # initialize defaults
    sobj <- list(
      normalization = "none", trial_subset = FALSE, add_deriv = FALSE, ts_multiplier = FALSE,
      demean_convolved = TRUE, beta_series = FALSE, value_type = "unit", value_fixed = 1
    )
    checkmate::assert_string(ss$name)
    checkmate::assert_string(ss$event)
    if (!ss$event %in% names(l1_model_set$events)) {
      stop(sprintf("Cannot locate event %s in l1_model_set$events", ss$event))
    }

    sobj$name <- ss$name
    sobj$event <- ss$event

    trial_set <- NULL
    if (!is.null(ss$trial_subset_expression)) {
      trial_set <- tryCatch(with(trial_data, eval(parse(text = ss$trial_subset_expression))),
        error = function(e) {
          lg$error("Problem evaluating trial subsetting expression: %s. Error: %s", ss$trial_subset_expression, as.character(e))
          return(NULL)
        }
      )
      if (!is.null(trial_set)) {
        sobj$trial_subset_expression <- ss$trial_subset_expression
        sobj$trial_subset <- TRUE
      }
    }

    # this gets populated in signals even with empty trial subset expression
    sobj <- get_trial_subset_stats(sobj, trial_data, trial_set)

    if (!is.null(ss$parametric_modulator)) {
      stopifnot(ss$parametric_modulator %in% slist$values)
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

    if (!is.null(ss$normalization)) {
      checkmate::assert_subset(ss$normalization, c("none", "evtmax_1", "durmax_1"))
      sobj$normalization <- ss$normalization
    }

    if (!is.null(ss$demean_convolved)) {
      checkmate::assert_logical(ss$demean_convolved, len = 1L)
      sobj$demean_convolved <- ss$demean_convolved
    }

    if (!is.null(ss$add_deriv)) {
      checkmate::assert_logical(ss$add_deriv, len = 1L)
      sobj$add_deriv <- ss$add_deriv
    }

    # TODO: this is a bit redundant with parts of bl1_specify_wi_factors... would be nice to unify
    if (!is.null(ss$wi_factors)) {
      checkmate::assert_subset(ss$wi_factors, l1_model_set$wi_factors)
      sobj$wi_factors <- ss$wi_factors

      ff <- NULL
      if (!is.null(ss$wi_formula)) {
        ff <- tryCatch(as.formula(ss$wi_formula), error = function(e) {
          wi_formula <- tryCatch(as.formula(wi_formula), error = function(e) {
            print(e)
            cat("Problem converting your syntax to formula. Defaulting to additive model\n")
            return(NULL)
          })
        })
      }

      if (is.null(ff)) {
        ff <- as.formula(paste("~", paste(ss$wi_factors, collapse = " + ")))
      }

      sobj$wi_formula <- as.character(ff)

      model_vars <- all.vars(ff)
      if (!all(model_vars %in% sobj$wi_factors)) {
        lg$warn("At least one variable in the spec wi_formula expression was not in wi_factors. Adding it!\n")
        sobj$wi_factors <- union(sobj$wi_factors, model_vars)
      }

      sobj$value <- get_value_df(sobj, trial_data, wi_factors = sobj$wi_factors)

      # fit dummy model to populate a set of dummy coefficients, then save those to the object
      sobj <- fit_wi_model(sobj)
    } else {
      sobj$value <- get_value_df(sobj, trial_data)
    }

    signal_list[[sobj$name]] <- sobj # this will overwrite existing specification
  }

  class(signal_list) <- c("list", "l1_model_set_signals") # set 'l1_model_set_signals' class
  l1_model_set$signals <- signal_list

  return(l1_model_set)
}

populate_event_data <- function(eobj, trial_data) {
  # basal columns for each event
  meta_cols <- c("id", "session", "run_number", "trial")

  if (is.numeric(eobj$duration)) {
    edata <- trial_data %>%
      dplyr::select(!!meta_cols, onset = !!eobj$onset) %>%
      #setNames(c("onset")) %>%
      dplyr::mutate(duration = !!eobj$duration, event = !!eobj$name)
  } else {
    edata <- trial_data %>%
      dplyr::select(!!meta_cols, onset = !!eobj$onset, duration = !!eobj$duration) %>%
      #setNames(c("onset", "duration")) %>%
      dplyr::mutate(event = !!eobj$name)
  }

  if (!is.null(eobj$isi)) {
    if (is.numeric(eobj$isi)) {
      edata$isi <- eobj$isi
    } else {
      idata <- trial_data %>%
        dplyr::select(!!eobj$isi) %>%
        setNames("isi")
      edata <- edata %>% bind_cols(idata) # add isi column
    }
  } else {
    edata$isi <- NA_real_ # populate NA
  }

  eobj$data <- edata #metadata_df %>% dplyr::bind_cols(edata)
  return(eobj)
}

# Helper function to propagate names of each element into the $name field for that element
# For example, the events should be named in the YAML specification. The user can also supply a name field for the event, which is almost
# always redundant. Thus, if the name field is not specified for an event, fall back to the name of the event in the event list
# Example:
#   events:
#     clock:
#       name: clock # this is redundant!
propagate_spec_names <- function(spec_list) {
  fields <- c("signals", "events", "l1_models")
  for (ff in fields) {
    this_l <- spec_list[[ff]]
    for (ii in seq_along(this_l)) {
      if (is.null(this_l[[ii]]$name)) {
        this_l[[ii]]$name <- names(this_l)[ii]
      }
    }
    spec_list[[ff]] <- this_l
  }
  return(spec_list)
}