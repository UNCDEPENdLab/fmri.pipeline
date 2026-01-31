#' Functions for extracting and summarizing collinearity diagnostics across the GLM pipeline
#'
#' These functions provide structured summaries of regressor correlations and variance
#' inflation factors (VIF) for all subjects and runs following model setup.
#'
#' @name collinearity_diagnostics
NULL

#' Extract collinearity diagnostics from all L1 models
#'
#' This function reads the saved build_design_matrix objects (d_obj) from all subjects

#' and models, extracting collinearity diagnostics (correlations and VIFs) into a
#' structured format suitable for review and quality control.
#'
#' @param gpa A glm_pipeline_arguments object that has been processed through setup_l1_models
#' @param l1_model_names Character vector of L1 model names to extract. If NULL (default),
#'   extracts from all available models.
#' @param vif_threshold Numeric threshold for flagging high VIF values. Default is 5.
#' @param cor_threshold Numeric threshold for flagging high correlations. Default is 0.8.
#'
#' @return A list with class "l1_collinearity_summary" containing:
#' \itemize{
#'   \item \code{vif_summary}: A data.frame with VIF values for each regressor, subject, session, run, and model.
#'     Includes a flag for VIFs exceeding the threshold.
#'   \item \code{correlation_summary}: A data.frame with pairwise correlations between regressors
#'     for each subject, session, run, and model. Includes a flag for correlations exceeding the threshold.
#'   \item \code{high_vif}: A data.frame subset of vif_summary where VIF exceeds the threshold.
#'   \item \code{high_correlation}: A data.frame subset of correlation_summary where |r| exceeds the threshold.
#'   \item \code{summary_stats}: Aggregate statistics across all subjects/runs.
#'   \item \code{thresholds}: The VIF and correlation thresholds used.
#' }
#'
#' @details
#' This function iterates through all subject/model combinations in the gpa object,
#' loading the cached design matrix objects and extracting collinearity information.
#' The resulting summaries can be used to identify problematic regressors or runs
#' that may have estimation issues due to multicollinearity.
#'
#' Common causes of high collinearity include:
#' \itemize{
#'   \item Parametric regressors that are highly correlated with task indicator regressors
#'   \item Events that occur too close together in time
#'   \item Insufficient variation in parametric modulators within a run
#' }
#'
#' @examples
#' \dontrun{
#'   # After running setup_l1_models
#'   collin_summary <- extract_l1_collinearity(gpa)
#'   
#'   # View high VIF regressors
#'   print(collin_summary$high_vif)
#'   
#'   # View highly correlated regressor pairs
#'   print(collin_summary$high_correlation)
#'   
#'   # Get overall summary
#'   print(collin_summary)
#' }
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr filter mutate arrange
#' @export
extract_l1_collinearity <- function(gpa, l1_model_names = NULL, vif_threshold = 5, cor_threshold = 0.8) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_number(vif_threshold, lower = 1)
  checkmate::assert_number(cor_threshold, lower = 0, upper = 1)
  
  # Default to all models if not specified
  if (is.null(l1_model_names)) {
    l1_model_names <- names(gpa$l1_models$models)
  }
  checkmate::assert_subset(l1_model_names, names(gpa$l1_models$models))
  
  # Initialize collectors
  vif_records <- list()
  cor_records <- list()
  
  # Get unique subject/session combinations
  subjects <- unique(gpa$subject_data[, c("id", "session")])
  
  for (i in seq_len(nrow(subjects))) {
    subj_id <- subjects$id[i]
    subj_session <- subjects$session[i]
    
    for (model_name in l1_model_names) {
      # Construct path to the BDM cache file
      l1_output_dir <- get_output_directory(
        id = subj_id, session = subj_session,
        l1_model = model_name, gpa = gpa, what = "l1"
      )
      
      bdm_file <- file.path(l1_output_dir, paste0(model_name, "_bdm_setup.RData"))
      
      if (!file.exists(bdm_file)) {
        warning(sprintf("BDM file not found for id=%s, session=%s, model=%s: %s",
                        subj_id, subj_session, model_name, bdm_file))
        next
      }
      
      # Load the d_obj
      d_obj <- NULL
      tryCatch({
        load(bdm_file)
      }, error = function(e) {
        warning(sprintf("Failed to load BDM file: %s. Error: %s", bdm_file, e$message))
      })
      
      if (is.null(d_obj) || is.null(d_obj$collin_convolved)) {
        warning(sprintf("No collinearity data in BDM file for id=%s, session=%s, model=%s",
                        subj_id, subj_session, model_name))
        next
      }
      
      # Extract VIF and correlation data for each run
      for (run_name in names(d_obj$collin_convolved)) {
        run_collin <- d_obj$collin_convolved[[run_name]]
        run_number <- as.integer(gsub("run_number", "", run_name))
        
        # Extract VIF data
        if (!is.null(run_collin$vif) && !all(is.na(run_collin$vif))) {
          vif_df <- data.frame(
            id = subj_id,
            session = subj_session,
            run_number = run_number,
            l1_model = model_name,
            regressor = names(run_collin$vif),
            vif = as.numeric(run_collin$vif),
            stringsAsFactors = FALSE
          )
          vif_records[[length(vif_records) + 1]] <- vif_df
        }
        
        # Extract correlation data (upper triangle only, excluding diagonal)
        if (!is.null(run_collin$r) && !all(is.na(run_collin$r))) {
          cor_mat <- run_collin$r
          reg_names <- colnames(cor_mat)
          
          # Get upper triangle indices
          for (j in seq_len(ncol(cor_mat) - 1)) {
            for (k in (j + 1):ncol(cor_mat)) {
              cor_df <- data.frame(
                id = subj_id,
                session = subj_session,
                run_number = run_number,
                l1_model = model_name,
                regressor1 = reg_names[j],
                regressor2 = reg_names[k],
                correlation = cor_mat[j, k],
                stringsAsFactors = FALSE
              )
              cor_records[[length(cor_records) + 1]] <- cor_df
            }
          }
        }
      }
    }
  }
  
  # Combine all records
  vif_summary <- if (length(vif_records) > 0) {
    data.table::rbindlist(vif_records) |> as.data.frame()
  } else {
    data.frame()
  }
  
  correlation_summary <- if (length(cor_records) > 0) {
    data.table::rbindlist(cor_records) |> as.data.frame()
  } else {
    data.frame()
  }
  
  # Add threshold flags
  if (nrow(vif_summary) > 0) {
    vif_summary$high_vif <- vif_summary$vif > vif_threshold
  }
  
  if (nrow(correlation_summary) > 0) {
    correlation_summary$abs_correlation <- abs(correlation_summary$correlation)
    correlation_summary$high_correlation <- correlation_summary$abs_correlation > cor_threshold
  }
  
  # Extract flagged cases
  high_vif <- if (nrow(vif_summary) > 0) {
    vif_summary[vif_summary$high_vif, , drop = FALSE]
  } else {
    data.frame()
  }
  
  high_correlation <- if (nrow(correlation_summary) > 0) {
    correlation_summary[correlation_summary$high_correlation, , drop = FALSE]
  } else {
    data.frame()
  }
  
  # Compute summary statistics
  summary_stats <- list(
    n_subjects = length(unique(paste(vif_summary$id, vif_summary$session))),
    n_runs = nrow(unique(vif_summary[, c("id", "session", "run_number", "l1_model")])),
    n_models = length(unique(vif_summary$l1_model)),
    n_regressors = length(unique(vif_summary$regressor)),
    vif_mean = if (nrow(vif_summary) > 0) mean(vif_summary$vif, na.rm = TRUE) else NA,
    vif_max = if (nrow(vif_summary) > 0) max(vif_summary$vif, na.rm = TRUE) else NA,
    vif_n_flagged = sum(vif_summary$high_vif, na.rm = TRUE),
    cor_mean_abs = if (nrow(correlation_summary) > 0) mean(correlation_summary$abs_correlation, na.rm = TRUE) else NA,
    cor_max_abs = if (nrow(correlation_summary) > 0) max(correlation_summary$abs_correlation, na.rm = TRUE) else NA,
    cor_n_flagged = sum(correlation_summary$high_correlation, na.rm = TRUE)
  )
  
  result <- list(
    vif_summary = vif_summary,
    correlation_summary = correlation_summary,
    high_vif = high_vif,
    high_correlation = high_correlation,
    summary_stats = summary_stats,
    thresholds = list(vif = vif_threshold, correlation = cor_threshold)
  )
  
  class(result) <- c("l1_collinearity_summary", "list")
  return(result)
}

#' Print method for l1_collinearity_summary objects
#'
#' @param x An l1_collinearity_summary object
#' @param ... Additional arguments (ignored)
#' @export
print.l1_collinearity_summary <- function(x, ...) {
  cat("L1 Collinearity Diagnostics Summary\n")
  cat("====================================\n\n")
  
  cat(sprintf("Subjects analyzed: %d\n", x$summary_stats$n_subjects))
  cat(sprintf("Total runs: %d\n", x$summary_stats$n_runs))
  cat(sprintf("L1 models: %d\n", x$summary_stats$n_models))
  cat(sprintf("Regressors: %d\n\n", x$summary_stats$n_regressors))
  
  cat("Variance Inflation Factors (VIF):\n")
  cat(sprintf("  Threshold: %.1f\n", x$thresholds$vif))
  cat(sprintf("  Mean VIF: %.2f\n", x$summary_stats$vif_mean))
  cat(sprintf("  Max VIF: %.2f\n", x$summary_stats$vif_max))
  cat(sprintf("  Flagged (VIF > %.1f): %d\n\n", x$thresholds$vif, x$summary_stats$vif_n_flagged))
  
  cat("Regressor Correlations:\n")
  cat(sprintf("  Threshold: |r| > %.2f\n", x$thresholds$correlation))
  cat(sprintf("  Mean |r|: %.3f\n", x$summary_stats$cor_mean_abs))
  cat(sprintf("  Max |r|: %.3f\n", x$summary_stats$cor_max_abs))
  cat(sprintf("  Flagged pairs: %d\n\n", x$summary_stats$cor_n_flagged))
  
  if (nrow(x$high_vif) > 0) {
    cat("High VIF regressors (showing first 10):\n")
    print(head(x$high_vif[order(-x$high_vif$vif), 
                          c("id", "session", "run_number", "l1_model", "regressor", "vif")], 10),
          row.names = FALSE)
    cat("\n")
  }
  
  if (nrow(x$high_correlation) > 0) {
    cat("Highly correlated regressor pairs (showing first 10):\n")
    print(head(x$high_correlation[order(-x$high_correlation$abs_correlation),
                                   c("id", "session", "run_number", "l1_model", 
                                     "regressor1", "regressor2", "correlation")], 10),
          row.names = FALSE)
    cat("\n")
  }
  
  invisible(x)
}

#' Summarize VIF by regressor across all subjects
#'
#' Creates a summary table showing VIF statistics for each regressor across all
#' subjects and runs.
#'
#' @param collin_summary An l1_collinearity_summary object from extract_l1_collinearity
#' @param by_model Logical. If TRUE, summarize separately by L1 model. Default FALSE.
#'
#' @return A data.frame with columns: regressor, mean_vif, sd_vif, min_vif, max_vif,
#'   n_runs, n_flagged, pct_flagged
#'
#' @export
summarize_vif_by_regressor <- function(collin_summary, by_model = FALSE) {
  checkmate::assert_class(collin_summary, "l1_collinearity_summary")
  
  if (nrow(collin_summary$vif_summary) == 0) {
    return(data.frame())
  }
  
  group_vars <- if (by_model) c("l1_model", "regressor") else "regressor"
  
  result <- collin_summary$vif_summary |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarize(
      mean_vif = mean(vif, na.rm = TRUE),
      sd_vif = sd(vif, na.rm = TRUE),
      min_vif = min(vif, na.rm = TRUE),
      max_vif = max(vif, na.rm = TRUE),
      n_runs = dplyr::n(),
      n_flagged = sum(high_vif, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(pct_flagged = 100 * n_flagged / n_runs) |>
    dplyr::arrange(dplyr::desc(mean_vif)) |>
    as.data.frame()
  
  return(result)
}

#' Summarize correlations by regressor pair across all subjects
#'
#' Creates a summary table showing correlation statistics for each regressor pair
#' across all subjects and runs.
#'
#' @param collin_summary An l1_collinearity_summary object from extract_l1_collinearity
#' @param by_model Logical. If TRUE, summarize separately by L1 model. Default FALSE.
#'
#' @return A data.frame with columns: regressor1, regressor2, mean_cor, sd_cor,
#'   min_cor, max_cor, mean_abs_cor, n_runs, n_flagged, pct_flagged
#'
#' @export
summarize_correlations_by_pair <- function(collin_summary, by_model = FALSE) {
  checkmate::assert_class(collin_summary, "l1_collinearity_summary")
  
  if (nrow(collin_summary$correlation_summary) == 0) {
    return(data.frame())
  }
  
  group_vars <- if (by_model) {
    c("l1_model", "regressor1", "regressor2")
  } else {
    c("regressor1", "regressor2")
  }
  
  result <- collin_summary$correlation_summary |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarize(
      mean_cor = mean(correlation, na.rm = TRUE),
      sd_cor = sd(correlation, na.rm = TRUE),
      min_cor = min(correlation, na.rm = TRUE),
      max_cor = max(correlation, na.rm = TRUE),
      mean_abs_cor = mean(abs_correlation, na.rm = TRUE),
      n_runs = dplyr::n(),
      n_flagged = sum(high_correlation, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(pct_flagged = 100 * n_flagged / n_runs) |>
    dplyr::arrange(dplyr::desc(mean_abs_cor)) |>
    as.data.frame()
  
  return(result)
}

#' Export collinearity diagnostics to CSV files
#'
#' Writes the collinearity summary data to CSV files for external review or
#' documentation purposes.
#'
#' @param collin_summary An l1_collinearity_summary object from extract_l1_collinearity
#' @param output_dir Directory where CSV files should be written. If NULL,
#'   uses the gpa output directory.
#' @param prefix Prefix for output file names. Default is "collinearity_".
#'
#' @return Invisibly returns a character vector of file paths that were written.
#'
#' @export
export_collinearity_to_csv <- function(collin_summary, output_dir, prefix = "collinearity_") {
  checkmate::assert_class(collin_summary, "l1_collinearity_summary")
  checkmate::assert_directory_exists(output_dir)
  checkmate::assert_string(prefix)
  
  files_written <- character()
  
  # Write VIF summary
  if (nrow(collin_summary$vif_summary) > 0) {
    vif_file <- file.path(output_dir, paste0(prefix, "vif_all.csv"))
    write.csv(collin_summary$vif_summary, file = vif_file, row.names = FALSE)
    files_written <- c(files_written, vif_file)
    message("Wrote VIF data to: ", vif_file)
  }
  
  # Write correlation summary
  if (nrow(collin_summary$correlation_summary) > 0) {
    cor_file <- file.path(output_dir, paste0(prefix, "correlations_all.csv"))
    write.csv(collin_summary$correlation_summary, file = cor_file, row.names = FALSE)
    files_written <- c(files_written, cor_file)
    message("Wrote correlation data to: ", cor_file)
  }
  
  # Write flagged cases
  if (nrow(collin_summary$high_vif) > 0) {
    high_vif_file <- file.path(output_dir, paste0(prefix, "high_vif.csv"))
    write.csv(collin_summary$high_vif, file = high_vif_file, row.names = FALSE)
    files_written <- c(files_written, high_vif_file)
    message("Wrote high VIF cases to: ", high_vif_file)
  }
  
  if (nrow(collin_summary$high_correlation) > 0) {
    high_cor_file <- file.path(output_dir, paste0(prefix, "high_correlations.csv"))
    write.csv(collin_summary$high_correlation, file = high_cor_file, row.names = FALSE)
    files_written <- c(files_written, high_cor_file)
    message("Wrote high correlation cases to: ", high_cor_file)
  }
  
  # Write regressor summaries
  vif_by_reg <- summarize_vif_by_regressor(collin_summary)
  if (nrow(vif_by_reg) > 0) {
    vif_reg_file <- file.path(output_dir, paste0(prefix, "vif_by_regressor.csv"))
    write.csv(vif_by_reg, file = vif_reg_file, row.names = FALSE)
    files_written <- c(files_written, vif_reg_file)
    message("Wrote VIF by regressor summary to: ", vif_reg_file)
  }
  
  cor_by_pair <- summarize_correlations_by_pair(collin_summary)
  if (nrow(cor_by_pair) > 0) {
    cor_pair_file <- file.path(output_dir, paste0(prefix, "correlations_by_pair.csv"))
    write.csv(cor_by_pair, file = cor_pair_file, row.names = FALSE)
    files_written <- c(files_written, cor_pair_file)
    message("Wrote correlations by pair summary to: ", cor_pair_file)
  }
  
  invisible(files_written)
}

#' Create a visualization of collinearity diagnostics
#'
#' Generates plots summarizing VIF values and correlations across subjects.
#'
#' @param collin_summary An l1_collinearity_summary object from extract_l1_collinearity
#' @param plot_type Character string specifying the type of plot:
#'   "vif" for VIF boxplots, "correlation" for correlation heatmap,
#'   "both" for both plots. Default is "both".
#'
#' @return A ggplot object or list of ggplot objects
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline geom_tile scale_fill_gradient2
#'   theme_minimal labs coord_flip theme element_text
#' @export
plot_collinearity <- function(collin_summary, plot_type = "both") {
  checkmate::assert_class(collin_summary, "l1_collinearity_summary")
  checkmate::assert_choice(plot_type, c("vif", "correlation", "both"))
  
  plots <- list()
  
  # VIF boxplot
  if (plot_type %in% c("vif", "both") && nrow(collin_summary$vif_summary) > 0) {
    p_vif <- ggplot2::ggplot(collin_summary$vif_summary, 
                              ggplot2::aes(x = reorder(regressor, vif, FUN = median), y = vif)) +
      ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.7) +
      ggplot2::geom_hline(yintercept = collin_summary$thresholds$vif, 
                          linetype = "dashed", color = "red", linewidth = 1) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Variance Inflation Factors by Regressor",
        subtitle = sprintf("Red line indicates threshold (VIF = %.1f)", collin_summary$thresholds$vif),
        x = "Regressor",
        y = "VIF"
      ) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
    
    plots$vif <- p_vif
  }
  
  # Correlation heatmap (average across subjects)
  if (plot_type %in% c("correlation", "both") && nrow(collin_summary$correlation_summary) > 0) {
    cor_avg <- summarize_correlations_by_pair(collin_summary)
    
    # Create a symmetric matrix for the heatmap
    all_regs <- unique(c(cor_avg$regressor1, cor_avg$regressor2))
    
    # Build data for heatmap including both directions
    heatmap_data <- rbind(
      cor_avg[, c("regressor1", "regressor2", "mean_cor")],
      data.frame(
        regressor1 = cor_avg$regressor2,
        regressor2 = cor_avg$regressor1,
        mean_cor = cor_avg$mean_cor
      )
    )
    
    p_cor <- ggplot2::ggplot(heatmap_data, 
                              ggplot2::aes(x = regressor1, y = regressor2, fill = mean_cor)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, limits = c(-1, 1),
        name = "Mean r"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Average Regressor Correlations Across Subjects",
        x = "", y = ""
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = ggplot2::element_text(size = 9)
      )
    
    plots$correlation <- p_cor
  }
  
  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }
}
