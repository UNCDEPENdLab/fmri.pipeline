#' Helper function to  summarize l1 models outside of build_l1_models()
#' @param l1_model_set an l1_model_set object that may have extant events in it
#' @param gpa the \code{gpa} object
#' @return a summary of l1_model_set with events added/updated
#' @importFrom dplyr select mutate
#' @importFrom checkmate test_number
#' @importFrom glue glue
#' @export
summarize_l1_models <- function(gpa) {

    checkmate::assert_list(gpa, null.ok = FALSE)

    l1_model_set <- gpa$l1_models

   summarize_events(l1_model_set)
    print("hi")
}

#' Copied over from build_l1_models.R, since it is an inner function
#' inside build_l1_models() and did not want to refract to that
summarize_events <- function(l1_model_set) {
      cat("Summary of events available in l1 models:\n--------\n")
      if (!is.null(l1_model_set$events)) {
      lapply(l1_model_set$events, function(x) {
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