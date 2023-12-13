# helper functions for UI

#' Round and format the mean, min, and max as a single character string
#' @param x the numeric variable to summarize
#' @param digits the number of digits to use in rounding
#' @param prespace An integer indicating the number of spaces to include at the beginning of the string. Default: 0
#' @param showNA whether to print information about NAs. As with \code{table}, the
#'   options are "no", "ifany", and "always".
#' @param newline if TRUE, include a newline character \code{\n} at the end of the string
#' @return a formatted character string with mean, min, and max
#' @keywords internal
ol_mean_min_max <- function(x, digits=2, prespace=0, showNA="ifany", newline=TRUE) {
  checkmate::assert_integerish(digits, lower = 0, upper = 100, len=1L)
  checkmate::assert_integerish(prespace, lower = 0, upper = 100, len=1L)
  checkmate::assert_subset(showNA, c("no", "ifany", "always"))
  checkmate::assert_flag(newline)

  ret <- paste(
    paste(rep(" ", prespace), collapse = ""),
    "mean [min -- max]: ",
    round(mean(x, na.rm = TRUE), digits), "[",
    round(min(x, na.rm = TRUE), digits), "--",
    round(max(x, na.rm = TRUE), digits), "]"
  )

  nas <- sum(is.na(x))
  if (showNA == "always" || (showNA == "ifany" && nas > 0)) {
    ret <- paste0(ret, "  (", nas, " NAs)")
  }

  if (isTRUE(newline)) ret <- paste0(ret, "\n")

  return(ret)
}

#' Print the top n unique values and frequencies of an atomic variable
#' @param x the variable to tabulate
#' @param max_show the maximum number of unique values to show
#' @param max_allowed the maximum number of allowed unique values before the function summarizes as mean, min, max
#' @param prespace An integer indicating the number of spaces to include at the beginning of the string. Default: 0
#' @param showNA whether to print information about NAs. As with \code{table}, the
#'   options are "no", "ifany", and "always".
#' @param newline if TRUE, include a newline character \code{\n} at the end of the string
#' @return a formatted character string with frequencies for the top \code{max_show} categories
#' @importFrom dplyr slice_max
#' @keywords internal
ol_unique <- function(x, max_show = 5, max_allowed = 25, prespace = 0, showNA = "ifany", newline = TRUE) {
  checkmate::assert_atomic(x)
  checkmate::assert_integerish(max_show, lower = 1, len = 1L)
  checkmate::assert_integerish(max_allowed, lower = 1, len = 1L)
  checkmate::assert_integerish(prespace, lower = 0, upper = 100, len = 1L)
  checkmate::assert_subset(showNA, c("no", "ifany", "always"))
  checkmate::assert_flag(newline)

  n_uniq <- length(unique(x))
  # revert to mean, min, max if too many unique values (suggestive of continuous variable)
  if (n_uniq > max_allowed) {
    return(ol_mean_min_max(x, prespace=prespace, showNA=showNA, newline=newline))
  } else if (n_uniq == 1L) {
    return(as.character(x[1L])) # just print the single value
  }

  top_x <- as.data.frame(table(x)) %>%
    slice_max(Freq, n = max_show) %>%
    mutate(Freq = paste0("(n=", Freq, ")")) %>%
    tidyr::unite("val_n", c(x, Freq), sep = " ") %>%
    pull(val_n) %>%
    paste(collapse=", ")

  ret <- paste0(
    paste(rep(" ", prespace), collapse = ""),
    paste0(
      "Top ", min(max_show, n_uniq), " values",
      ifelse(n_uniq > max_show, paste0(" (of ", n_uniq, "): "), ": ")
    ),
    top_x
  )

  nas <- sum(is.na(x))
  if (showNA == "always" || (showNA == "ifany" && nas > 0)) {
    ret <- paste0(ret, "  (", nas, " NAs)")
  }

  if (isTRUE(newline)) ret <- paste0(ret, "\n")

  return(ret)
}