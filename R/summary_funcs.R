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
          cat("  -  ", gpa$parallel$compute_environment[2], "\n")
          cat(" R:\n")
          cat("  -  ", gpa$parallel$compute_environment[1], "\n")
          cat(" FSL:\n")
          cat("  -  ", gpa$parallel$fsl$compute_environment, "\n")
          cat(" AFNI:\n")
          cat("  -  ", gpa$parallel$afni$compute_environment, "\n")
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


model_summary <- function(gpa_object) {
  complete <- FALSE
      while (isFALSE(complete)) {
        model_object <- menu(c(
          "Events",
          "Signals",
          "Models",
          "All of the Above",
          "Exit Summary Build"), title= "Choose L1 Model Object")
          if (model_object == 1) {
            cat("\n \n")
            summarize_events(gpa_object$l1_models)
          } else if (model_object == 2) {
            cat("\n \n")
            summ_l1_sig(gpa_object$l1_models$signals)
            cat("\n")
            y_n <- menu(c("Yes", "No"), title="More Intensive Signals Summary?")
            if (y_n == 1) {
              summarize_l1_signals(gpa_object$l1_models$signals)
            } else if (y_n == 2) {
                cat("\nGoodbye.\n\n")
            }
          } else if (model_object == 3) {
            cat("\n \n")
            summarize_l1_models(gpa_object$l1_models$models)
          } else if (model_object == 4) {
            cat("\n \n")
            summarize_events(gpa_object$l1_models)
            summarize_l1_signals(gpa_object$l1_models$signals)
            summarize_l1_models(gpa_object$l1_models$models)
          } else {
            cat("\n\nGoodbye.\n")
          }
          complete <- TRUE
        }
      }

summ_l1_sig <- function(sl) {
  if (length(sl) == 0L) {
    return(invisible(NULL))
  }
  cat("Summary of signals available in l1 models:")
  cat("\n--------\n")
  lapply(seq_along(sl), function(ii) {
    this <- sl[[ii]]
    cat("Signal", ii, "of", length(sl), ":\n")
    cat("  - Name:", this$name, "\n")
    cat("  - Event alignment:", this$event, "\n")
    if (isFALSE(this$trial_subset)) {
      cat("  - Trial subset: FALSE", "\n")
    } else if (isTRUE(this$trial_subset)) {
      cat(glue("  - Trial subset: {this$trial_subset_expression}"), "\n")
    }
    if (this$value_type %in% c("unit", "number")) {
      cat("  - Regressor value (constant): ", this$value_fixed[1L], "\n")
    } else {
      cat(
        "  - Parametric value: ", this$parametric_modulator, ", mean [min -- max]: ",
        round(mean(this$value$value, na.rm = TRUE), 2), " [ ",
        round(min(this$value$value, na.rm = TRUE), 2), " -- ",
        round(max(this$value$value, na.rm = TRUE), 2), " ] \n", sep=""
      )
    }
    }
  )
}

install.packages("yaml")
library(yaml)
l1_yaml <- function(gpa){
  gpa2 <- gpa$l1_models
  ## Model Specifications
  end_yaml <- list(
  onset = gpa2$onset,
  durations = gpa2$durations,
  isis = gpa2$isis,
  wi_factors = gpa2$wi_factors,
  values = gpa2$values
  )
  eobj <- list()
  for (ee in gpa2$events) {
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
    end_yaml$events[[eobj$name]] <- eobj # this will overwrite existing specification
  }
  sobj <- list()
  for (ss in gpa2$signals) {
    sobj$name <- ss$name
    sobj$event <- ss$event
    sobj$trial_subset_expression <- ss$trial_subset_expression
    if (!is.null(ss$normalization)) {
      checkmate::assert_subset(ss$normalization, c("none", "evtmax_1", "durmax_1"))
      sobj$normalization <- ss$normalization
    }
    if (!is.null(ss$parametric_modulator)) {
      stopifnot(ss$parametric_modulator %in% gpa2$values)
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
    end_yaml$signals[[sobj$name]] <- sobj
  }
  mobj <- list()
  for (mm in gpa2$models) {
    mobj$name <- mm$name
    mobj$signals <- mm$signals
    end_yaml$models[[mobj$name]] <- mobj
  }
  endyaml3 <- as.yaml(end_yaml)
  yaml_choice <- menu(c(
          "Console Output",
          "File Output",
          "Both",
          "Exit"), title= "How would you like to receive the YAML file?")
          if (yaml_choice == 1) {
            return(cat(endyaml3, "\nGoodbye.\n"))
          } else if (yaml_choice == 2) {
            writeLines(endyaml3, "output.yaml")
            return(cat("\nFile should be seen as \"output.yaml\"."))
          } else if (yaml_choice == 3) {
            writeLines(endyaml3, "output.yaml")
            return(endyaml3, "\nFile should be seen as \"output.yaml\".")
          } else (
            return(cat("\nGoodbye.\n"))
          )
}