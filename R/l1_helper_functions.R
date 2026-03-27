# helper function to return a logical vector of trials included in subset
get_trial_set_from_signal <- function(ss, trial_data=NULL) {
  checkmate::assert_data_frame(trial_data)
  # if user has a cached trial expression, need to evaluate the trial set before getting value df
  if (isTRUE(ss$trial_subset)) {
    # calculate trial set from cached expression
    trial_set <- with(trial_data, eval(parse(text = ss$trial_subset_expression)))
  } else {
    # when FALSE: keep all trials. This is a local variable that is calculated every time this function is called (e.g., modification)
    trial_set <- rep(TRUE, nrow(trial_data))
  }
  return(trial_set)
}

# helper function to compute statistics for how much data is excluded by trial subset
get_trial_subset_stats <- function(signal, trial_data, trial_set=NULL) {
  if (is.null(trial_set)) {
    trial_set <- get_trial_set_from_signal(signal, trial_data)
  }

  overall <- sum(trial_set == TRUE) / length(trial_set)
  tmp <- trial_data %>% bind_cols(trial_set = trial_set)
  by_id <- tmp %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(pct_true = sum(trial_set == TRUE) / n(), .groups = "drop") %>%
    dplyr::summarize(
      mean = mean(pct_true, na.rm = T), sd = sd(pct_true, na.rm = T),
      min = min(pct_true, na.rm = T), max = max(pct_true, na.rm = T)
    ) %>%
    mutate(overall = overall) %>%
      unlist()

  signal$trial_subset_statistics <- by_id
  return(signal)
}

# helper function to figure out expected regressor columns for a given signal based on whether
# the signal has within-subject factors, is a beta series signal, and/or includes a temporal derivative
get_regressors_from_signal <- function(sig) {
  # in terms of design, always add derivative columns en bloc after the corresponding non-derivative columns
  # handle expansion of PPI signals slightly differently -- expand w/i factors within psych and PPI parts if needed
  if (checkmate::test_string(sig$ts_multiplier)) {
    if (!is.null(sig$wi_factors)) {
      cols <- paste(sig$name, names(coef(sig$wi_model)), sep = ".")
      cols <- c(cols, sig$ts_multiplier, sub(paste0("^", sig$name), paste(sig$name, "ppi", sep = "."), cols))
    } else {
      cols <- c(sig$name, sig$ts_multiplier, paste(sig$name, "ppi", sep = "."))
    }
  } else if (!is.null(sig$wi_factors)) {
    # use the lmfit object in the signal to determine the columns that will be included
    # for within-subject factor modulation, always include the signal as the prefix on the column names
    cols <- paste(sig$name, names(coef(sig$wi_model)), sep = ".")
  } else if (isTRUE(sig$beta_series)) {
    if (is.data.frame(sig$value)) {
      # TODO: this approach is imperfect if there are jumps in trials for a subject
      # this assumes that all subjects have all trials
      trials <- sort(unique(sig$value$trial)) # vector of trials for parametric signal
    } else {
      # trials will be in corresponding event in case value is a scalar
      trials <- sort(unique(event_list[[sig$event]]$trial))
    }

    cols <- paste(sig$name, sprintf("%03d", trials), sep = "_t") # signal_t001 etc.
  } else {
    cols <- sig$name # nothing special
  }

  if (isTRUE(sig$add_deriv)) {
    cols <- c(cols, paste0("d_", cols)) # add temporal derivative for each column
  }
  return(cols)
}

# cf. get_regressors_from_signal -- these have to align!
# vv <- expand_signal(gpa$l1_models$signals$pe)
expand_signal <- function(sig) {
  if (checkmate::test_string(sig$ts_multiplier) && is.null(sig$ts_expanded)) {
    # need to have two copies of sig and expand each of those -- one with the ts_multiplier and one without
    sig$ts_expanded <- TRUE # used to indicate that it shouldn't be recursively expanded below
    sig_psych <- sig
    sig_psych$ts_multiplier <- NULL # remove this to prevent multiplication by ppi data
    sig_psych <- expand_signal(sig_psych) # handle expansion due to wi_factors etc.

    sig_physio <- sig
    sig_physio$name <- sig_physio$ts_multiplier # name signal after physio column
    sig_physio$physio_only <- TRUE # physio only marker
    sig_physio <- list(sig_physio)
    names(sig_physio) <- sig_physio[[1]]$name

    # interaction signal should be expanded, too, but with physio multiplication
    sig_ppi <- sig
    sig_ppi$name <- paste(sig_ppi$name, "ppi", sep = ".")
    sig_ppi <- expand_signal(sig_ppi)

    s_list <- c(sig_psych, sig_physio, sig_ppi)
  } else if (is.null(sig$wi_factors) && isFALSE(sig$beta_series)) {
    # nothing to expand
    s_list <- list(sig) # wrap as list to ensure we always get a list return
    names(s_list) <- sig$name
  } else if (!is.null(sig$wi_factors) && isTRUE(sig$beta_series)) {
    stop("In expand_signal, cannot sort out what to do when beta series is enabled and there is a within-subject factor!")
  } else if (!is.null(sig$wi_factors)) {
    # extract within-subject regression model matrix (has dummy coding)
    model_df <- as.data.frame(model.matrix(sig$wi_model))
    if (nrow(model_df) != nrow(sig$value)) {
      warning("Number of rows in wi_model data.frame does not match signal value data.frame")
    }

    # for each column of the model matrix, pull the signal value data.frame where
    # the dummy code is 1 (i.e., rows to which the dummy pertains). This divides up the value
    # data.frame into smaller data.frames, one per dummy combination
    dummy_value <- lapply(model_df, function(x) {
      which_1 <- which(abs(x - 1) < 1e-5)
      sig$value %>%
        dplyr::slice(!!which_1) %>%
        dplyr::select(-!!sig$wi_factors)
    })

    # coefficients to loop over
    coefs <- names(coef(sig$wi_model))

    # expand regressor names to match l1 model specification
    exp_names <- paste(sig$name, names(coef(sig$wi_model)), sep = ".")

    # loop over each coefficient and build the signal copy
    s_list <- lapply(seq_along(exp_names), function(ee) {
      sig_copy <- sig

      # update signal name to reflect level of within-subject factor(s)
      sig_copy$name <- exp_names[ee]

      # copy across relevant data.frame for this level
      sig_copy$value <- dummy_value[[ coefs[ee] ]]

      # this gets tricky because the $value has already been trimmed based on the parent subset
      # if (isTRUE(sig$trial_subset)) {
      #   ss <- get_trial_subset_stats(sig_copy, trial_data, trial_set)
      # }

      # remove residual within-subject modulation settings
      sig_copy$wi_formula <- sig_copy$wi_model <- sig_copy$wi_factors <- NULL
      return(sig_copy)
    })
  } else if (isTRUE(sig$beta_series)) {
    # for now, only support data.frame inputs for signal value
    # in general, these are expanded in setup_l1_models
    checkmate::assert_data_frame(sig$value)
    trials <- sort(unique(sig$value$trial)) # vector of trials for parametric signal
    exp_names <- paste(sig$name, sprintf("%03d", trials), sep = "_t") # signal_t001 etc.
    s_list <- lapply(seq_along(exp_names), function(ee) {
      sig_copy <- sig

      # update signal name to reflect level of within-subject factor(s)
      sig_copy$name <- exp_names[ee]

      # copy across relevant data.frame for this level
      sig_copy$value <- sig$value %>% filter(trial == !!trials[ee])
      sig_copy$beta_series <- FALSE # each copy is a single-trial event
      return(sig_copy)
    })
  }

  # create derivative signals en bloc
  if (isTRUE(sig$add_deriv)) {
    s_list_deriv <- lapply(s_list, function(ss) {
      ss$name <- paste0("d_", ss$name)
      ss$add_deriv <- FALSE
      ss$return_deriv <- TRUE # internal flag to tell convolution to compute and return deriv
      return(ss)
    })
    s_list <- c(s_list, s_list_deriv)
  }

  return(s_list)
}

fit_wi_model <- function(sobj) {
  if (!all(c("wi_factors", "wi_formula") %in% names(sobj))) {
    return(sobj) # cannot proceed with fitting
  }

  wi_df <- sobj$value %>%
    mutate(dummy = 1:n()) %>%
    mutate(across(!!sobj$wi_factors, factor)) # always force wi_factors to be stored as factor to make contrasts straightforward
  
  ffit <- update.formula(as.formula(sobj$wi_formula), "dummy ~ .")
  sobj$wi_model <- lm(ffit, wi_df)
  return(sobj)
}
  
#' S3 method to convert formula to character
#' @param x a formula to convert to character
#' @param ... not used
#' @details copied from formula.tools to avoid NAMESPACE import
#' @export
as.character.formula <- function(x, ...) { 
  form <- paste( deparse(x), collapse=" " )
  form <- gsub( "\\s+", " ", form, perl=FALSE ) # remove multiple spaces
  return(form)
}

#obsolete
# # helper function to expand beta series
# expand_signals <- function(orig) {
#   # use levels of wi_factors here to build out multiple signals, one per level... maybe just copy-paste these inside object, subsetting value?
#   signals_expanded <- lapply(orig, function(ss) {
#     if (is.null(ss$wi_factors)) {
#       return(ss)
#     } else {
#       # extract within-subject regression model matrix (has dummy coding)
#       model_df <- as.data.frame(model.matrix(ss$wi_model))

#       # for each column of the model matrix, pull the signal value data.frame where
#       # the dummy code is 1 (i.e., rows to which the dummy pertains). This divides up the value
#       # data.frame into smaller data.frames, one per dummy combination
#       dummy_value <- lapply(model_df, function(x) {
#         which_1 <- which(x == 1)
#         ss$value %>%
#           dplyr::slice(!!which_1) %>%
#           dplyr::select(-!!ss$wi_factors)
#       })

#       # create a list of signals, one per dummy level, by copying the master signal
#       # settings, then amending the value field to include only relevant trials/rows
#       signal_list <- lapply(seq_along(dummy_value), function(x) {
#         sig_copy <- ss
#         # update signal name to reflect level of within-subject factor(s)
#         sig_copy$name <- paste(sig_copy$name, make.names(names(dummy_value)[x]), sep = "_")

#         # copy across relevant data.frame for this level
#         sig_copy$value <- dummy_value[[x]]

#         # remove residual within-subject modulation settings
#         sig_copy$wi_formula <- sig_copy$wi_model <- sig_copy$wi_factors <- NULL
#         sig_copy
#       })
#     }
#   })

#   # now handle beta series
# }

