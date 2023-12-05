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
          cat(" -  GLM Settings: ", gpa$glm_settings, "\n")
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
          cat("  -  Feat Time: ", "\n")
          cat("  -  Feat Memory (GB): ", "\n")
          cat(" FSL Level 2:\n")
          cat("  -  Feat Time: ", "\n")
          cat("  -  Feat Memory (GB): ", "\n")
          cat(" FSL Level 3:\n")
          cat("  -  Feat Time: ", "\n")
          cat("  -  Feat Memory (GB): ", "\n")
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