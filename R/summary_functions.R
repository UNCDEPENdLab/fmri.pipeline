#' Helper function to better summarize GPA object
#' @param gpa the \code{gpa} object
#' @return a readable summary of the given gpa object
#' @export
summary.glm_pipeline_arguments <- function(gpa) {
      cat("\nSummary of GLM Pipeline Analysis: \n\n")
      if (!is.null(gpa)) {
          cat("Analysis Name: ", gpa$analysis_name[1])
          cat("\n  Node Name: ", gpa$nodename[1])
          cat("\n  Finalize Complete: ", gpa$finalize_complete[1])
        # GLM and Confound Settings
          cat("\n--------\n Settings: \n")
          cat(" -  Scheduler Name: ", gpa$scheduler, "\n")
          cat(" -  TR: ", ifelse(is.null(gpa$run_data$tr), gpa$tr, ol_unique(gpa$run_data$tr)), "\n")
          cat(" -  Multi-Level Run: ", gpa$multi_run[1], "\n")
          #cat(" -  GLM Settings: ", gpa$glm_settings, "\n")
          if (!is.null(gpa$confound_settings$na_strings)) {
            cat(" -  NA Strings: ", paste(gpa$confound_settings$na_strings, collapse = ", "), "\n")
          }
          if (!is.null(gpa$confound_settings$exclude_run)) {
            cat(" -  Exclude Runs: ", gpa$confound_settings$exclude_run, "\n")
          }
          if (!is.null(gpa$confound_settings$exclude_subject)) {
            cat(" -  Exclude Subjects: ", gpa$confound_settings$exclude_subject, "\n")
          } 
          if (!is.null(gpa$confound_settings$truncate_run)) {
            cat(" -  Truncate Run: ", gpa$confound_settings$truncate_run, "\n")
          }
          if (!is.null(gpa$confound_settings$spike_volumes)) {
            cat(" -  Spike Volumes: ", gpa$confound_settings$spike_volumes, "\n")
          }
          # VM/Key Columns
          cat("--------\n Key Columns: \n")
          cat("  ", paste(gpa$vm, collapse = ", "), "\n")
          cat("  n Expected Runs: ", gpa$n_expected_runs, "\n")
          # Models: L1/2/3
          cat("--------\n Models: \n")
          if (!is.null(gpa$l1_models)) {
            cat(" -  L1 Models: ", "Present \n")
          } else {
            cat(" -  L1 Models: NA\n")
          }
          if (!is.null(gpa$l1_models)) {
            cat(" -  L2 Models: ", "Present \n")
          } else {
            cat(" -  L2 Models: NA\n")
          }
          if (!is.null(gpa$l1_models)) {
            cat(" -  L3 Models: ", "Present \n")
          } else {
            cat(" -  L3 Models: NA\n")
          }
          # Computation Info
          cat("--------\n Computation Info: \n")
          cat(" FSL Level 1:\n")
          cat("  -  Feat Time: ", gpa$parallel$fsl$l1_feat_time, "\n")
          cat("  -  Feat Memory (GB): ", gpa$parallel$fsl$l1_feat_memgb, "\n")
          cat(" FSL Level 2:\n")
          cat("  -  Feat Time: ", gpa$parallel$fsl$l2_feat_time, "\n")
          cat("  -  Feat Memory (GB): ", gpa$parallel$fsl$l2_feat_memgb, "\n")
          cat(" FSL Level 3:\n")
          cat("  -  Feat Time: ", gpa$parallel$fsl$l3_feat_time, "\n")
          cat("  -  Feat Memory (GB): ", gpa$parallel$fsl$l3_feat_memgb, "\n")
          cat("  -  Feat CPUs per Job: ", gpa$parallel$fsl$l3_feat_cpusperjob, "\n")
          # Compute Environment Info
          cat("--------\n Computing Environment: \n")
          cat(" Global:\n")
          cat("  -  ", gpa$parallel$compute_environment$global, "\n")
          cat(" R:\n")
          cat("  -  ", gpa$parallel$compute_environment$r[2], "\n")
          cat(" FSL:\n")
          cat("  -  ", gpa$parallel$compute_environment$fsl[2], "\n")
          cat(" AFNI:\n")
          cat("  -  ", gpa$parallel$compute_environment$afni[2], "\n")
          # Output Directory and System Info
          cat("--------\n System Info: \n")
          cat(" -  Output Directory: ", gpa$output_directory[1], "\n")
          cat(" -  System Name: ", gpa$sys_info[1], "\n")
          cat(" -  Release: ", gpa$sys_info[2], "\n")
          cat(" -  Version: ", gpa$sys_info[3], "\n")
          cat(" -  Node Name: ", gpa$sys_info[4], "\n")
          cat(" -  Machine: ", gpa$sys_info[5], "\n")
          cat(" -  Login: ", gpa$sys_info[6], "\n")
          cat(" -  User: ", gpa$sys_info[7], "\n")
          cat(" -  Effective User: ", gpa$sys_info[8], "\n")
      } else {
      cat("No info in gpa object\n")
    }
}

#' Helper function to better summarize GPA object's L1 models
#' @param gpa the \code{gpa} object
#' @return a readable summary of the given gpa object's L1 models
#' @export
summary.l1_model_set <- function(l1_ms) {
  complete <- FALSE
  if (!is.null(l1_ms)) {
      while (isFALSE(complete)) {
        model_object <- menu(c(
          "Events",
          "Signals",
          "Models",
          "All of the Above",
          "Exit Summary Build"), title= "Choose L1 Model Object")
          if (model_object == 1) {
            cat("\n \n")
            summary(l1_ms$events)
          } else if (model_object == 2) {
            summarize_l1_signals(l1_ms$signals)
          } else if (model_object == 3) {
            cat("\n \n")
            summary(l1_ms$models)
          } else if (model_object == 4) {
            cat("\n \n")
            summary(l1_ms$events)
            summarize_l1_signals(l1_ms$signals)
            summary(l1_ms$models)
          } else {
            cat("\n\nGoodbye.\n")
          }
          complete <- TRUE
    }
  }
}

#' Helper function to better summarize GPA object's l1_models' events
#' @param gpa the \code{gpa} object
#' @return a readable summary of the given gpa object's l1_models' events
#' @export
summary.l1_model_set_events <- function(l1_model_set_events) {
    cat("Summary of events available in l1 models:\n--------\n")
    if (!is.null(l1_model_set_events)) {
      lapply(l1_model_set_events, function(x) {
        cat("Event name:     ", x$name, "\n")
        cat("  Event onset:    ", x$onset, "\n")
        cat(
          "    mean [min -- max]: ",
          round(mean(x$data$onset, na.rm = TRUE), 2), "[",
          round(min(x$data$onset, na.rm = TRUE), 2), "--",
          round(max(x$data$onset, na.rm = TRUE), 2), "]\n"
        )
        cat("    First 6 values:", paste(head(x$data$onset), collapse = ", "), "\n\n")
        cat("  Event duration: ", x$duration, "\n")
        cat(
          "    mean [min -- max]: ",
          round(mean(x$data$duration, na.rm = TRUE), 2), "[",
          round(min(x$data$duration, na.rm = TRUE), 2), "--",
          round(max(x$data$duration, na.rm = TRUE), 2), "]\n"
        )
        cat("    First 6 values:", paste(round(head(x$data$duration), 2), collapse = ", "), "\n\n")
        if (!is.null(x$isi)) {
          cat("  Event ISI/ITI:    ", x$isi, "\n")
          cat(
            "    mean [min -- max]: ",
            round(mean(x$data$isi, na.rm = TRUE), 2), "[",
            round(min(x$data$isi, na.rm = TRUE), 2), "--",
            round(max(x$data$isi, na.rm = TRUE), 2), "]\n"
          )
        }
      })
    } else {
      cat("No events in l1 models\n")
    }
}

#' Helper function to better summarize GPA object's l1_models' models
#' @param gpa the \code{gpa} object
#' @return a readable summary of the given gpa object's l1_models' models
#' @export
summary.l1_model_set_models <- function(ml) {
  if (length(ml) == 0L) { return(invisible(NULL)) }
  cat("Summary of l1 models:\n--------\n")
  lapply(seq_along(ml), function(ii) {
    this <- ml[[ii]]
    cat("--------\nModel ", ii, "\n\n")
    cat("  Name:", this$name, "\n")
    cat("  Signals:", paste(this$signals, collapse=", "), "\n")
    if (ncol(this$contrasts) < 20) {
      cat("  Contrasts:\n\n")
      print(round(this$contrasts, 3))
      cat("\n--------\n\n")
    } else { cat("More than 20 regressor columns. Not printing contrasts\n") }
  })
}

#' Helper function to better summarize GPA object's l1_models' signals
#' @param gpa the \code{gpa} object
#' @return a readable summary of the given gpa object's l1_models' signals
#' @export
summary.l1_model_set_signals <- function(sl) {
  if (length(sl) == 0L) {
    return(invisible(NULL))
  }
  cat("Summary of signals available in l1 models:\n--------\n")
  lapply(seq_along(sl), function(ii) {
    this <- sl[[ii]]
    cat("--------\nSignal ", ii, "\n\n")
    cat("  Name:", this$name, "\n")
    cat("  Event alignment:", this$event, "\n")
    if (isFALSE(this$trial_subset)) {
      cat("  Trial subset: FALSE", "\n")
    } else if (isTRUE(this$trial_subset)) {
      cat(glue("  Trial subset: {this$trial_subset_expression}"), "\n")
      cat(
        glue("  Proportion of trials included: overall = {round(this$trial_subset_statistics['overall'], 2)}"),
        glue(
          "    By id: mean = {round(this$trial_subset_statistics['mean'], 2)}, ", "sd = {round(this$trial_subset_statistics['sd'], 2)}, ",
          "min = {round(this$trial_subset_statistics['min'], 2)}, ", "max = {round(this$trial_subset_statistics['max'], 2)}"
        ), sep="\n"
      )
    }
    if (this$value_type %in% c("unit", "number")) {
      cat("  Regressor value (constant): ", this$value_fixed[1L], "\n")
    } else {
      cat(
        "  Parametric value: ", this$parametric_modulator, ", mean [min -- max]: ",
        round(mean(this$value$value, na.rm = TRUE), 2), " [ ",
        round(min(this$value$value, na.rm = TRUE), 2), " -- ",
        round(max(this$value$value, na.rm = TRUE), 2), " ] \n", sep=""
      )
    }
    cat("  HRF Normalization:", this$normalization, "\n")
    cat("  Add signal derivative:", as.character(this$add_deriv), "\n")
    cat("  Demean convolved signal:", as.character(this$demean_convolved), "\n")
    cat("  Within-subject factor:", c_string(this$wi_factors), "\n")
    cat("  Generate beta series:", as.character(this$beta_series), "\n")
    cat(
      "  Multiply convolved regressor against time series:",
      ifelse(is.null(this$ts_multiplier || isFALSE(this$ts_multiplier)),
        "FALSE", this$ts_multiplier
      ), "\n"
    )
    cat("\n--------\n")
  })
}

