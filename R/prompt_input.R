#' Obtain user input from the console
#' @param prompt The character string to place on the line preceding the user input prompt. For example, "Enter location"
#' @param prompt_eol The character string to place at the end of the prompt line. For example, ">"
#' @param instruct The instructions to display above the prompt.
#' @param lower For numeric inputs, the lowest valid value
#' @param upper For numeric inputs, the highest valid value
#' @param len The number of expected values to be returned. If NULL, the user can enter any number of values.
#' @param min.len The minimum number of values to be returned. If NULL, the user can enter any number of values.
#' @param max.len The maximum number of values to be returned. If NULL, the user can enter any number of values.
#' @param split The character(s) to split the input string into multiple values. Only relevant if len > 1.
#' @param among A vector of valid values for the input. If NULL, any value is accepted.
#' @param required If TRUE, the user must provide a value. If FALSE, the user can skip the input by pressing Enter.
#' @param uniq If TRUE, all entries must be unique.
#' @return The user input, converted to the appropriate type (numeric, integer, or character).
#' @details The function will keep prompting the user until valid input is provided. It will also display
#'   instructions and feedback about the expected input format.
#' @note This function is intended for interactive use and may not work as expected in non-interactive
#'   environments (e.g., R scripts run in batch mode).
#' @importFrom glue glue
#' @importFrom checkmate test_string assert_string assert_subset assert_number
#' @keywords internal
prompt_input <- function(prompt = "", prompt_eol=">", instruct = NULL, type = "character", lower = -Inf, upper = Inf, 
                         len = NULL, min.len=NULL, max.len=NULL, split = NULL, among = NULL, required = TRUE, uniq=FALSE, default = NULL) {
  
  if (!interactive()) stop("prompt_input() requires an interactive session.")
  
  if (is.null(prompt)) prompt <- ""
  if (is.null(prompt_eol)) prompt_eol <- ""
  
  checkmate::assert_string(prompt)
  checkmate::assert_string(prompt_eol)
  checkmate::assert_string(instruct, null.ok = TRUE)
  checkmate::assert_subset(type, c("numeric", "integer", "character", "file", "flag"))
  checkmate::assert_number(lower)
  checkmate::assert_number(upper)
  checkmate::assert_number(len, lower = 1L, null.ok = TRUE)
  checkmate::assert_number(min.len, lower = 1L, null.ok = TRUE)
  checkmate::assert_number(max.len, lower = 1L, null.ok = TRUE)
  if (!is.null(min.len) && !is.null(max.len) && max.len < min.len) {
    stop("max.len must be greater than or equal to min.len")
  }
  
  if (type == "flag" && !is.null(len) && len > 1L) {
    warning("Ignoring len > 1 for type 'flag' -- only one return supported")
    len <- 1L
  }
  
  if (!is.null(len) && len > 1 && !checkmate::test_string(split)) {
    stop("When multiple return values are required, you must specify character(s) to split the input.")
  }
  
  # setup feedback about the number and type of inputs expected
  inp_expect <- ""
  plural <- "s"
  if (!is.null(len)) {
    n_expect <- ifelse(len == 1L, "a", glue("{len}"))
    if (len==1L) plural <- ""
  } else if ((is.null(min.len) || min.len == 1L) && (is.null(max.len) || is.infinite(max.len))) {
    n_expect <- ""
  } else if (is.null(min.len) || min.len == 1L) {
    # max only
    n_expect <- glue("no more than {max.len}")
  } else if (is.null(max.len) || is.infinite(max.len)) {
    # min only
    n_expect <- glue("at least {min.len}")
  } else {
    # min and max
    n_expect <- glue("{min.len}-{max.len}")
  }
  
  if (nchar(n_expect) > 0L) n_expect <- paste0(n_expect, " ") # add spacing for formatting
  
  if (type=="integer") {
    if (n_expect=="a ") n_expect <- "an "
    inp_expect <- glue("Input must be {n_expect}integer{plural} between {lower} and {upper}\n")
  } else if (type == "numeric") {
    inp_expect <- glue("Input must be {n_expect}number{plural} between {lower} and {upper}\n")
  } else if (type == "character") {
    inp_expect <- glue("Input must be {n_expect}string{plural} separated by '{split}'\n")
  }
  
  # add trailing space
  if (stringr::str_length(prompt) != 0 && !grepl("\\s$", prompt)) {
    prompt <- paste0(prompt, " ")
  }
  
  # add options for flag prompt
  if (type == "flag") {
    # always ask user for yes/no input
    prompt <- paste0(prompt, ifelse(required, "(yes/no)", "(yes/no; press Enter to skip)"))
  } else if (!is.null(default)) {
    prompt <- glue::glue("{prompt}(Press enter to accept default: {default})") # let user know how to skip optional input
  } else if (!required) {
    prompt <- paste0(prompt, "(Press enter to skip)") # let user know how to skip optional input
  }
  
  # always add trailing space to make prompt clear
  if (stringr::str_length(prompt) != 0 && !grepl("\\s$", prompt)) prompt <- paste0(prompt, " ")
  if (!grepl("\\s$", prompt_eol)) prompt_eol <- paste0(prompt_eol, " ") # also ensure that prompt_eol has trailing space
  prompt <- paste0(prompt, prompt_eol)
  
  # Validate default value
  if (!is.null(default)) {
    valid_default <- switch(type,
                            "integer" = checkmate::test_integerish(default, lower = lower, upper = upper, len = len, min.len = min.len, max.len = max.len),
                            "numeric" = checkmate::test_numeric(default, lower = lower, upper = upper, len = len, min.len = min.len, max.len = max.len),
                            "character" = checkmate::test_character(default, len = len, min.len = min.len, max.len = max.len),
                            "flag" = is.logical(default) && length(default) == 1,
                            "file" = all(sapply(default, checkmate::test_file_exists)),
                            FALSE
    )
    if (!valid_default) stop("Default value does not meet the input requirements.")
  }
  
  # print instructions
  if (checkmate::test_string(instruct)) cat(instruct, "\n")
  
  # obtain user input
  res <- ""
  while (is.na(res[1L]) || res[1L] == "") {
    r <- readline(prompt)
    if (!is.null(split)) r <- strsplit(r, split, perl = TRUE)[[1]]
    
    if (!is.null(default) && r[1L] == "") {
      return(default)
    } else if (isFALSE(required) && r[1L] == "") {
      empty <- switch(type,
                      "integer" = NA_integer_,
                      "numeric" = NA_real_,
                      "character" = NA_character_,
                      "flag" = NA,
                      "file" = NA_character_
      )
      return(empty) # empty input and not required
    } else if (isTRUE(uniq) && length(unique(r)) != length(r)) {
      cat("All entries must be unique.\n")
    } else if (type == "flag") {
      r <- tolower(r)
      if (!r[1L] %in% c("yes", "no", "y", "n")) {
        cat("Please respond yes or no.\n")
      } else {
        res <- substr(r[1L], 1, 1) == "y" # TRUE if yes, FALSE if no
      }
    } else if (type == "integer") {
      r <- type.convert(r, as.is = TRUE) # convert to apparent atomic type for validation
      if (!checkmate::test_integerish(r, lower = lower, upper = upper, len = len, min.len = min.len, max.len = max.len)) {
        cat(inp_expect)
      } else {
        if (!is.null(among) && !all(r %in% among)) {
          cat(glue("Input must be integers in the set: {paste(among, collapse=', ')}\n"))
        } else {
          res <- r
        }
      }
    } else if (type == "numeric") {
      r <- type.convert(r, as.is = TRUE) # convert to apparent atomic type for validation
      if (!checkmate::test_numeric(r, lower = lower, upper = upper, len = len, min.len = min.len, max.len = max.len)) {
        cat(inp_expect)
      } else {
        res <- r
      }
    } else if (type == "character") {
      if (!checkmate::test_character(r, len = len, min.len = min.len, max.len = max.len)) {
        cat(inp_expect)
      } else {
        if (!is.null(among) && !all(r %in% among)) {
          cat(glue("Input must be {len} strings in the set: {paste(among, collapse=', ')}\n"))
        } else {
          res <- r
        }
      }
    } else if (type == "file") {
      # should probably think harder about unquoted filenames containing spaces throwing off len
      exist <- sapply(r, checkmate::test_file_exists)
      if (!all(exist)) {
        cat(glue("The following files could not be found: {paste(r[!exist], collapse=', ')}\n"))
      } else {
        res <- r
      }
    }
  }
  
  return(res)
}


# Pretty print a list with indentation and line wrapping
pretty_print_list <- function(x, indent = 0, width = 80) {
  indent_str <- strrep("  ", indent)
  
  for (name in names(x)) {
    value <- x[[name]]
    
    if (is.list(value)) {
      cat(sprintf("%s%s:\n", indent_str, name))
      pretty_print_list(value, indent + 1, width)
    } else {
      # Format value as character string
      value_str <- paste(value, collapse = ", ")
      
      # Wrap long lines
      wrapped <- strwrap(value_str,
                         width = width - nchar(indent_str) - nchar(name) - 2,
                         exdent = 2, simplify = FALSE
      )[[1]]
      
      # Print wrapped lines
      cat(sprintf("%s%s: %s\n", indent_str, name, wrapped[1]))
      if (length(wrapped) > 1) {
        for (line in wrapped[-1]) {
          cat(sprintf("%s%s  %s\n", indent_str, strrep(" ", nchar(name)), line))
        }
      }
    }
  }
}
