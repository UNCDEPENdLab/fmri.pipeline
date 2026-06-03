format_3dlmer_datatable_keys <- function(df, id_col = "id", session_col = "session", max_items = 8L) {
  if (nrow(df) == 0L) return("<none>")

  keys <- paste0(
    "id=", as.character(df[[id_col]]),
    ", session=", as.character(df[[session_col]])
  )
  keys <- unique(keys)
  if (length(keys) > max_items) {
    keys <- c(keys[seq_len(max_items)], sprintf("... (%d more)", length(keys) - max_items))
  }

  paste(keys, collapse = "; ")
}

format_3dlmer_datatable_paths <- function(paths, max_items = 8L) {
  paths <- unique(as.character(paths))
  if (length(paths) == 0L) return("<none>")
  if (length(paths) > max_items) {
    paths <- c(paths[seq_len(max_items)], sprintf("... (%d more)", length(paths) - max_items))
  }

  paste(paths, collapse = "; ")
}

is_missing_3dlmer_field <- function(x) {
  is.na(x) | !nzchar(trimws(as.character(x)))
}

find_3dlmer_rows_with_missing_fields <- function(df, fields) {
  if (length(fields) == 0L) return(logical(nrow(df)))

  missing_matrix <- vapply(
    fields,
    function(field) is_missing_3dlmer_field(df[[field]]),
    logical(nrow(df))
  )
  if (is.null(dim(missing_matrix))) return(missing_matrix)

  rowSums(missing_matrix) > 0L
}

find_3dlmer_duplicate_keys <- function(df, id_col = "id", session_col = "session") {
  keys <- paste(as.character(df[[id_col]]), as.character(df[[session_col]]), sep = "\r")
  duplicated(keys) | duplicated(keys, fromLast = TRUE)
}

#' Internal function to build the data table for AFNI 3dLMEr
#'
#' @param subject_data a data.frame containing subject/session-level covariates with unique id/session rows
#' @param input_files a data.frame with columns 'id', 'session', and 'InputFile'
#' @param model_variables a character vector of variables that must be included in the table
#'
#' @return a data.frame formatted for 3dLMEr -dataTable
#' @keywords internal
build_3dlmer_datatable <- function(subject_data, input_files, model_variables) {
  checkmate::assert_data_frame(subject_data)
  checkmate::assert_data_frame(input_files)
  checkmate::assert_character(model_variables)
  checkmate::assert_subset(c("id", "session", model_variables), names(subject_data))
  checkmate::assert_subset(c("id", "session", "InputFile"), names(input_files))

  subject_data <- subject_data %>%
    dplyr::select(dplyr::all_of(c("id", "session", model_variables)))

  input_files <- input_files %>%
    dplyr::select(dplyr::all_of(c("id", "session", "InputFile")))

  if (nrow(input_files) == 0L) {
    stop("AFNI 3dLMEr dataTable has no input rows.", call. = FALSE)
  }
  if (nrow(subject_data) == 0L) {
    stop("AFNI 3dLMEr subject_data has no rows.", call. = FALSE)
  }

  missing_input_rows <- find_3dlmer_rows_with_missing_fields(input_files, c("id", "session", "InputFile"))
  if (any(missing_input_rows)) {
    stop(
      sprintf(
        "AFNI 3dLMEr input_files has missing id, session, or InputFile values for rows: %s.",
        paste(which(missing_input_rows), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_subject_keys <- find_3dlmer_rows_with_missing_fields(subject_data, c("id", "session"))
  if (any(missing_subject_keys)) {
    stop(
      sprintf(
        "AFNI 3dLMEr subject_data has missing id or session values for rows: %s.",
        paste(which(missing_subject_keys), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_covariates <- find_3dlmer_rows_with_missing_fields(subject_data, model_variables)
  if (any(missing_covariates)) {
    stop(
      sprintf(
        "AFNI 3dLMEr subject_data has missing model variable values for id/session keys: %s.",
        format_3dlmer_datatable_keys(subject_data[missing_covariates, , drop = FALSE])
      ),
      call. = FALSE
    )
  }

  duplicate_input_keys <- find_3dlmer_duplicate_keys(input_files)
  if (any(duplicate_input_keys)) {
    stop(
      sprintf(
        "Duplicate AFNI 3dLMEr input rows for id/session keys: %s. Each contrast dataTable must have exactly one InputFile per id/session.",
        format_3dlmer_datatable_keys(input_files[duplicate_input_keys, , drop = FALSE])
      ),
      call. = FALSE
    )
  }

  duplicate_subject_keys <- find_3dlmer_duplicate_keys(subject_data)
  if (any(duplicate_subject_keys)) {
    stop(
      sprintf(
        "Duplicate AFNI 3dLMEr subject_data rows for id/session keys: %s. Covariates must be one-to-one with AFNI input rows.",
        format_3dlmer_datatable_keys(subject_data[duplicate_subject_keys, , drop = FALSE])
      ),
      call. = FALSE
    )
  }

  input_paths <- as.character(input_files$InputFile)
  missing_files <- input_paths[!file.exists(input_paths)]
  if (length(missing_files) > 0L) {
    stop(
      sprintf(
        "AFNI 3dLMEr InputFile entries do not exist: %s.",
        format_3dlmer_datatable_paths(missing_files)
      ),
      call. = FALSE
    )
  }

  missing_subject_rows <- input_files %>%
    dplyr::anti_join(subject_data %>% dplyr::select(dplyr::all_of(c("id", "session"))), by = c("id", "session"))
  if (nrow(missing_subject_rows) > 0L) {
    stop(
      sprintf(
        "AFNI 3dLMEr dataTable cannot be built because subject_data is missing id/session rows for input keys: %s.",
        format_3dlmer_datatable_keys(missing_subject_rows)
      ),
      call. = FALSE
    )
  }

  # Merge input files with subject data to get covariates
  # subject_data should contain 'id' and 'session' for longitudinal
  dt <- input_files %>%
    dplyr::left_join(subject_data, by = c("id", "session"))

  if (nrow(dt) != nrow(input_files)) {
    stop(
      sprintf(
        "AFNI 3dLMEr dataTable row count changed during merge: %d input row(s), %d merged row(s).",
        nrow(input_files), nrow(dt)
      ),
      call. = FALSE
    )
  }

  # Keep only necessary columns: Subj, model_variables, InputFile
  # 3dLMEr requires the column name 'Subj' for the subject ID
  dt <- dt %>%
    dplyr::rename(Subj = id) %>%
    dplyr::select(dplyr::all_of(c("Subj", "session", model_variables, "InputFile")))

  duplicate_final_keys <- find_3dlmer_duplicate_keys(dt, id_col = "Subj")
  if (any(duplicate_final_keys)) {
    stop(
      sprintf(
        "AFNI 3dLMEr dataTable contains duplicate Subj/session rows after merge: %s.",
        format_3dlmer_datatable_keys(dt[duplicate_final_keys, , drop = FALSE], id_col = "Subj")
      ),
      call. = FALSE
    )
  }

  return(dt)
}

build_3dlmer_intersection_mask <- function(mask_files, output_file, lg = NULL) {
  checkmate::assert_character(mask_files, min.len = 1L)
  checkmate::assert_file_exists(mask_files)
  checkmate::assert_string(output_file)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  mask_files <- unique(mask_files)
  if (is.null(lg)) lg <- lgr::get_logger("fmri.pipeline/afni_3dlmer/mask")
  if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

  mask_imgs <- lapply(mask_files, RNifti::readNifti)
  mask_dims <- vapply(mask_imgs, function(img) paste(dim(img), collapse = "x"), character(1))
  if (length(unique(mask_dims)) != 1L) {
    stop(
      sprintf(
        "Cannot build AFNI 3dLMEr intersection mask from mismatched FSL masks: %s",
        paste(unique(mask_dims), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  mask_4d <- do.call(abind::abind, c(mask_imgs, list(along = 4)))
  mask_sum <- rowSums(mask_4d, dims = 3)
  intersect_img <- 1L * (mask_sum == dim(mask_4d)[4L])

  header <- RNifti::niftiHeader(mask_files[1L])
  intersect_nifti <- RNifti::asNifti(intersect_img, header)
  RNifti::writeNifti(intersect_nifti, file = output_file)

  lg$info(
    "Computed AFNI 3dLMEr intersection mask from %d FSL L2 masks: %s",
    length(mask_files), output_file
  )

  output_file
}

resolve_3dlmer_mask <- function(target_dir, harvested_inputs, explicit_mask = NULL, requires_l2 = FALSE, lg = NULL) {
  checkmate::assert_string(target_dir)
  checkmate::assert_data_frame(harvested_inputs)
  checkmate::assert_string(explicit_mask, null.ok = TRUE)
  checkmate::assert_flag(requires_l2)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (!is.null(explicit_mask)) return(explicit_mask)
  if (!isTRUE(requires_l2) || !"feat_dir" %in% names(harvested_inputs)) return(NULL)

  mask_files <- unique(file.path(as.character(harvested_inputs$feat_dir), "mask.nii.gz"))
  if (length(mask_files) == 0L || anyNA(mask_files) || any(!file.exists(mask_files))) {
    if (!is.null(lg)) {
      lg$warn("Could not derive AFNI 3dLMEr intersection mask because one or more FSL L2 masks were missing.")
    }
    return(NULL)
  }

  build_3dlmer_intersection_mask(
    mask_files = mask_files,
    output_file = file.path(target_dir, "intersection_mask.nii.gz"),
    lg = lg
  )
}

#' Internal function to construct the 3dLMEr command string
#'
#' @param prefix output prefix for 3dLMEr
#' @param model_formula fixed effects formula string
#' @param qVars character vector of quantitative variables
#' @param glt_codes list of named GLT code strings
#' @param data_table_file path to the data table text file
#' @param mask path to the brain mask file
#' @param njobs number of parallel jobs for 3dLMEr
#' @param ss_type sum of squares type (default 3)
#'
#' @return a character string containing the 3dLMEr command
#' @keywords internal
# Collapse whitespace because AFNI counts spaces inside -model formulas.
sanitize_3dlmer_model_formula <- function(model_formula) {
  checkmate::assert_string(model_formula)
  gsub("[[:space:]]+", "", model_formula)
}

validate_3dlmer_formula_datatable <- function(model_formula, datatable, context = NULL) {
  checkmate::assert_string(model_formula)
  checkmate::assert_data_frame(datatable)
  checkmate::assert_string(context, null.ok = TRUE)

  formula_text <- trimws(model_formula)
  if (!grepl("~", formula_text, fixed = TRUE)) {
    formula_text <- paste("~", formula_text)
  }

  parsed_formula <- tryCatch(
    stats::as.formula(formula_text, env = baseenv()),
    error = function(e) {
      stop(
        sprintf(
          "Could not parse AFNI 3dLMEr lmer_formula%s: %s",
          if (!is.null(context)) paste0(" for ", context) else "",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  formula_vars <- unique(all.vars(parsed_formula))
  missing_vars <- setdiff(formula_vars, names(datatable))
  if (length(missing_vars) > 0L) {
    msg <- sprintf(
      "AFNI 3dLMEr lmer_formula%s references variable(s) not present in dataTable: %s. Available dataTable columns: %s.",
      if (!is.null(context)) paste0(" for ", context) else "",
      paste(missing_vars, collapse = ", "),
      paste(names(datatable), collapse = ", ")
    )
    if ("id" %in% missing_vars && "Subj" %in% names(datatable)) {
      msg <- paste(
        msg,
        "Use 'Subj' instead of 'id' in 3dLMEr random-effects terms because the AFNI dataTable renames subject IDs to Subj."
      )
    }
    stop(msg, call. = FALSE)
  }

  invisible(formula_vars)
}

sanitize_3dlmer_glt_names <- function(glt_names) {
  checkmate::assert_character(glt_names, any.missing = TRUE)

  glt_names[is.na(glt_names)] <- ""
  glt_names <- trimws(glt_names)
  glt_names <- gsub("[^A-Za-z0-9_]+", "_", glt_names)
  glt_names <- gsub("_+", "_", glt_names)
  glt_names <- gsub("^_+|_+$", "", glt_names)

  empty_names <- !nzchar(glt_names)
  if (any(empty_names)) {
    glt_names[empty_names] <- paste0("glt", which(empty_names))
  }

  starts_badly <- !grepl("^[A-Za-z]", glt_names)
  glt_names[starts_badly] <- paste0("glt_", glt_names[starts_badly])
  make.unique(glt_names, sep = "_")
}

empty_3dlmer_glt_table <- function() {
  data.frame(
    label_raw = character(0),
    label_afni = character(0),
    code = character(0),
    source = character(0),
    contrast_type = character(0),
    valid = logical(0),
    message = character(0),
    stringsAsFactors = FALSE
  )
}

as_3dlmer_glt_table <- function(glt_codes, source = "user", contrast_type = "raw") {
  if (is.null(glt_codes) || length(glt_codes) == 0L) return(empty_3dlmer_glt_table())

  if (is.data.frame(glt_codes)) {
    required <- c("label_raw", "code")
    missing_cols <- setdiff(required, names(glt_codes))
    if (length(missing_cols) > 0L) {
      stop("3dLMEr GLT table is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    out <- as.data.frame(glt_codes, stringsAsFactors = FALSE)
    if (!"label_afni" %in% names(out)) out$label_afni <- sanitize_3dlmer_glt_names(out$label_raw)
    if (!"source" %in% names(out)) out$source <- source
    if (!"contrast_type" %in% names(out)) out$contrast_type <- contrast_type
    if (!"valid" %in% names(out)) out$valid <- TRUE
    if (!"message" %in% names(out)) out$message <- ""
    return(out[, names(empty_3dlmer_glt_table()), drop = FALSE])
  }

  if (is.character(glt_codes)) {
    glt_codes <- as.list(glt_codes)
  }
  checkmate::assert_list(glt_codes)
  labels <- names(glt_codes)
  if (is.null(labels)) labels <- paste0("glt", seq_along(glt_codes))
  labels[is.na(labels) | !nzchar(labels)] <- paste0("glt", which(is.na(labels) | !nzchar(labels)))

  data.frame(
    label_raw = labels,
    label_afni = sanitize_3dlmer_glt_names(labels),
    code = as.character(unlist(glt_codes, use.names = FALSE)),
    source = source,
    contrast_type = contrast_type,
    valid = TRUE,
    message = "",
    stringsAsFactors = FALSE
  )
}

parse_3dlmer_glt_blocks <- function(code, datatable_names) {
  checkmate::assert_string(code)
  checkmate::assert_character(datatable_names, min.len = 1L)

  var_pattern <- paste(gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", datatable_names), collapse = "|")
  matches <- gregexpr(paste0("(?<![A-Za-z0-9_.])(", var_pattern, ")\\s*:"), code, perl = TRUE)[[1L]]
  if (identical(matches[1L], -1L)) {
    return(list(blocks = list(), error = "No 'variable :' block found."))
  }

  pre <- substr(code, 1L, matches[1L] - 1L)
  if (nzchar(trimws(pre))) {
    return(list(blocks = list(), error = sprintf("Unexpected text before first GLT variable block: %s", trimws(pre))))
  }

  lengths <- attr(matches, "match.length")
  headers <- regmatches(code, list(matches))[[1L]]
  vars <- trimws(sub(":\\s*$", "", headers))
  blocks <- vector("list", length(matches))
  for (ii in seq_along(matches)) {
    start <- matches[ii] + lengths[ii]
    end <- if (ii < length(matches)) matches[ii + 1L] - 1L else nchar(code)
    blocks[[ii]] <- list(var = vars[ii], expr = trimws(substr(code, start, end)))
  }

  list(blocks = blocks, error = "")
}

validate_3dlmer_glt_table <- function(glt_table, datatable, qVars = NULL, context = NULL) {
  checkmate::assert_data_frame(datatable)
  checkmate::assert_character(qVars, null.ok = TRUE)
  checkmate::assert_string(context, null.ok = TRUE)

  glt_table <- as_3dlmer_glt_table(glt_table)
  if (nrow(glt_table) == 0L) return(glt_table)

  qVars <- qVars %||% character(0)
  datatable_names <- names(datatable)
  numeric_pattern <- "[+-]?(?:(?:[0-9]+(?:\\.[0-9]*)?)|(?:\\.[0-9]+))(?:[eE][+-]?[0-9]+)?"
  errors <- character(nrow(glt_table))

  for (ii in seq_len(nrow(glt_table))) {
    code <- glt_table$code[ii]
    row_errors <- character(0)

    if (is.na(code) || !nzchar(trimws(code))) {
      row_errors <- c(row_errors, "GLT code is empty.")
    } else {
      parsed <- parse_3dlmer_glt_blocks(code, datatable_names)
      if (nzchar(parsed$error)) {
        row_errors <- c(row_errors, parsed$error)
      } else {
        for (block in parsed$blocks) {
          var <- block$var
          expr <- block$expr
          if (!var %in% datatable_names) {
            row_errors <- c(row_errors, sprintf("Variable '%s' is not in dataTable.", var))
            next
          }

          is_qvar <- var %in% qVars || is.numeric(datatable[[var]])
          tokens <- strsplit(expr, "\\s+")[[1L]]
          tokens <- tokens[nzchar(tokens)]
          if (length(tokens) == 0L) {
            row_errors <- c(row_errors, sprintf("Variable '%s' has no GLT weight expression.", var))
            next
          }

          if (is_qvar) {
            bad <- !grepl(paste0("^", numeric_pattern, "$"), tokens, perl = TRUE)
            if (any(bad) || length(tokens) != 1L) {
              row_errors <- c(row_errors, sprintf(
                "Quantitative variable '%s' must use one numeric weight, e.g. '%s : 1'.",
                var, var
              ))
            }
          } else {
            bad_syntax <- !grepl(paste0("^", numeric_pattern, "\\*.+$"), tokens, perl = TRUE)
            if (any(bad_syntax)) {
              row_errors <- c(row_errors, sprintf(
                "Categorical variable '%s' has invalid level weight token(s): %s.",
                var, paste(tokens[bad_syntax], collapse = ", ")
              ))
              next
            }
            levels_seen <- sub(paste0("^", numeric_pattern, "\\*"), "", tokens, perl = TRUE)
            known_levels <- unique(as.character(datatable[[var]]))
            missing_levels <- setdiff(levels_seen, known_levels)
            if (length(missing_levels) > 0L) {
              row_errors <- c(row_errors, sprintf(
                "Categorical variable '%s' references level(s) not present in dataTable: %s.",
                var, paste(missing_levels, collapse = ", ")
              ))
            }
            if (identical(glt_table$contrast_type[ii], "pairwise_diff")) {
              weights <- as.numeric(sub("\\*.*$", "", tokens, perl = TRUE))
              if (length(weights) < 2L || !isTRUE(all.equal(sum(weights), 0, tolerance = 1e-8))) {
                row_errors <- c(row_errors, sprintf(
                  "Pairwise AFNI 3dLMEr GLT '%s' must use explicit cell-level weights that sum to zero for '%s' (for example, '%s : 1*A -1*B'), not treatment-coded coefficient weights.",
                  glt_table$label_raw[ii], var, var
                ))
              }
            }
          }
        }
      }
    }

    errors[ii] <- paste(row_errors, collapse = " ")
  }

  glt_table$valid <- !nzchar(errors)
  glt_table$message <- errors
  if (any(!glt_table$valid)) {
    bad <- glt_table[!glt_table$valid, , drop = FALSE]
    stop(
      sprintf(
        "Invalid AFNI 3dLMEr GLT code%s: %s",
        if (!is.null(context)) paste0(" for ", context) else "",
        paste(sprintf("%s: %s", bad$label_raw, bad$message), collapse = " | ")
      ),
      call. = FALSE
    )
  }

  glt_table
}

build_3dlmer_command <- function(prefix, model_formula, qVars = NULL, glt_codes = NULL,
                                 data_table_file, mask = NULL, njobs = 1, ss_type = 3) {
  checkmate::assert_string(prefix)
  checkmate::assert_string(model_formula)
  checkmate::assert_character(qVars, null.ok = TRUE)
  if (!is.null(glt_codes) &&
      !checkmate::test_list(glt_codes) &&
      !checkmate::test_data_frame(glt_codes) &&
      !checkmate::test_character(glt_codes)) {
    stop("glt_codes must be a list, character vector, data.frame, or NULL.", call. = FALSE)
  }
  checkmate::assert_string(data_table_file)
  checkmate::assert_string(mask, null.ok = TRUE)
  checkmate::assert_integerish(njobs, lower = 1)
  checkmate::assert_integerish(ss_type, lower = 1, upper = 3)

  model_formula <- sanitize_3dlmer_model_formula(model_formula)

  cmd <- glue::glue("3dLMEr -prefix {shQuote(prefix)} -model {shQuote(model_formula)}")

  if (!is.null(qVars) && length(qVars) > 0) {
    cmd <- paste(cmd, glue::glue("-qVars {shQuote(paste(qVars, collapse=\",\"))}"))
  }

  if (!is.null(mask)) {
    cmd <- paste(cmd, glue::glue("-mask {shQuote(mask)}"))
  }

  if (njobs > 1) {
    cmd <- paste(cmd, glue::glue("-jobs {njobs}"))
  }

  cmd <- paste(cmd, glue::glue("-SS_type {ss_type}"))

  if (!is.null(glt_codes) && length(glt_codes) > 0) {
    glt_table <- as_3dlmer_glt_table(glt_codes)
    for (i in seq_len(nrow(glt_table))) {
      glt_name <- glt_table$label_afni[i]
      glt_code <- glt_table$code[i]
      # glt_code should already be in AFNI format: "Factor : 1*Level1 -1*Level2"
      cmd <- paste(cmd, glue::glue("-gltCode {glt_name} {shQuote(glt_code)}"))
    }
  }

  cmd <- paste(cmd, glue::glue("-dataTable {shQuote(paste0('@', data_table_file))}"))

  return(cmd)
}

#' Helper to write 3dLMEr script and data table
#'
#' @param output_dir directory where files should be written
#' @param dt the data table data.frame
#' @param cmd the 3dLMEr command string
#' @param script_name name of the shell script to write
#'
#' @return list with paths to script and data table
#' @keywords internal
write_3dlmer_files <- function(output_dir, dt, cmd, script_name = "run_3dlmer.sh") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  dt_file <- file.path(output_dir, "dataTable.txt")
  # write.table with proper formatting for AFNI (space separated, header)
  write.table(dt, file = dt_file, row.names = FALSE, quote = FALSE, sep = "\t")

  script_file <- file.path(output_dir, script_name)
  # Basic script structure
  script_content <- c(
    "#!/bin/bash",
    "set -uo pipefail",
    "",
    "script_dir=\"$(cd \"$(dirname \"$0\")\" && pwd)\"",
    "cd \"$script_dir\"",
    "",
    "complete_file=\".afni_complete\"",
    "fail_file=\".afni_fail\"",
    "start_file=\".afni_start\"",
    "log_file=\"3dLMEr.log\"",
    "rm -f \"$complete_file\" \"$fail_file\"",
    "start_time=\"$(date -Is)\"",
    "printf '%s\\n' \"$start_time\" > \"$start_file\"",
    "",
    "{",
    "  printf '%s\\n' \"[$start_time] Starting AFNI 3dLMEr\"",
    paste0("  printf '%s\\n' ", shQuote(paste("Command:", cmd))),
    paste0("  ", cmd),
    "} >> \"$log_file\" 2>&1",
    "exit_code=$?",
    "end_time=\"$(date -Is)\"",
    "",
    "if [ $exit_code -eq 0 ]; then",
    "  printf '%s\\n%s\\n' \"$start_time\" \"$end_time\" > \"$complete_file\"",
    "  rm -f \"$fail_file\"",
    "else",
    "  printf '%s\\n%s\\n' \"$start_time\" \"$end_time\" > \"$fail_file\"",
    "fi",
    "",
    "exit $exit_code"
  )
  writeLines(script_content, script_file)
  Sys.chmod(script_file, "0755")

  return(list(script = script_file, data_table = dt_file))
}

#' Status checker for 3dLMEr
#'
#' @param output_file expected output file from 3dLMEr
#' @param lg optional logger
#' @param prefix optional prefix for returned column names
#' @return data.frame indicating output existence and completion status
#' @keywords internal
get_3dlmer_status <- function(output_file, lg = NULL, prefix = NULL) {
  checkmate::assert_string(output_file)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)
  checkmate::assert_string(prefix, null.ok = TRUE)
  if (is.null(lg)) lg <- lgr::get_logger()

  resolved_output <- output_file
  candidates <- unique(c(
    output_file,
    paste0(output_file, ".HEAD"),
    paste0(output_file, ".BRIK"),
    paste0(output_file, "+tlrc.HEAD"),
    paste0(output_file, "+orig.HEAD")
  ))

  if (grepl("\\+tlrc$|\\+orig$", output_file)) {
    candidates <- unique(c(
      output_file,
      paste0(output_file, ".HEAD"),
      paste0(output_file, ".BRIK")
    ))
  } else if (grepl("\\.(nii|nii\\.gz|HEAD|BRIK)$", output_file, perl = TRUE)) {
    sans_ext <- file_sans_ext(output_file)
    if (!is.na(sans_ext)) {
      candidates <- unique(c(
        output_file,
        paste0(sans_ext, "+tlrc.HEAD"),
        paste0(sans_ext, "+orig.HEAD")
      ))
    }
  }

  existing <- candidates[file.exists(candidates)]
  output_exists <- length(existing) > 0L
  if (output_exists) {
    resolved_output <- existing[1L]
  } else {
    lg$debug("No AFNI 3dLMEr output found for '%s'. Checked: %s", output_file, paste(candidates, collapse = ", "))
  }

  output_dir <- dirname(output_file)
  start_file <- file.path(output_dir, ".afni_start")
  complete_file <- file.path(output_dir, ".afni_complete")
  fail_file <- file.path(output_dir, ".afni_fail")
  log_file <- file.path(output_dir, "3dLMEr.log")
  start_file_exists <- file.exists(start_file)
  complete_file_exists <- file.exists(complete_file)
  fail_file_exists <- file.exists(fail_file)

  execution_start <- as.POSIXct(NA)
  execution_end <- as.POSIXct(NA)
  execution_min <- NA_real_
  timing_file <- NULL

  if (complete_file_exists) {
    timing_file <- complete_file
    if (fail_file_exists) {
      lg$warn("Both .afni_complete and .afni_fail files exist in %s", output_dir)
      lg$warn(
        "Treating this AFNI 3dLMEr run as %s based on output files.",
        ifelse(isTRUE(output_exists), "complete", "failed")
      )
    }
    if (!isTRUE(output_exists)) {
      lg$warn("Found .afni_complete in %s, but expected AFNI 3dLMEr output is missing.", output_dir)
    }
    afni_complete <- isTRUE(output_exists)
    afni_failed <- !isTRUE(output_exists)
  } else if (fail_file_exists) {
    lg$debug("Detected AFNI 3dLMEr failure marker in: %s", output_dir)
    timing_file <- fail_file
    afni_complete <- FALSE
    afni_failed <- TRUE
  } else {
    afni_complete <- output_exists
    afni_failed <- if (output_exists) FALSE else NA
    if (start_file_exists) timing_file <- start_file
  }

  if (!is.null(timing_file)) {
    timing <- readLines(timing_file, warn = FALSE)
    if (length(timing) > 0L) {
      timing <- anytime::anytime(timing)
      execution_start <- timing[1L]
      if (length(timing) >= 2L) {
        execution_end <- timing[2L]
        execution_min <- as.numeric(difftime(timing[2L], timing[1L], units = "mins"))
      } else if (!identical(timing_file, start_file)) {
        lg$warn("Did not find two timing entries in %s.", timing_file)
        lg$warn("File contents: %s", paste(timing, collapse = "; "))
      }
    }
  }

  out <- data.frame(
    output_file = output_file,
    afni_output_file_resolved = resolved_output,
    afni_output_file_exists = output_exists,
    afni_start_file = start_file,
    afni_start_file_exists = start_file_exists,
    afni_complete_file = complete_file,
    afni_complete_file_exists = complete_file_exists,
    afni_fail_file = fail_file,
    afni_fail_file_exists = fail_file_exists,
    afni_log_file = log_file,
    afni_log_file_exists = file.exists(log_file),
    afni_execution_start = execution_start,
    afni_execution_end = execution_end,
    afni_execution_min = execution_min,
    afni_complete = afni_complete,
    afni_failed = afni_failed,
    stringsAsFactors = FALSE
  )
  if (!is.null(prefix)) names(out) <- paste0(prefix, names(out))
  out
}

#' Translator from emmeans-style contrast spec to 3dLMEr -gltCode
#'
#' @param mobj high-level model specification object (hi_model_spec)
#' @param data the data.frame used for the model (to check factor levels)
#'
#' @return a list of 3dLMEr gltCode strings
#' @keywords internal
emmeans_to_3dlmer_glt <- function(mobj, data, qVars = NULL, raw_glt_codes = NULL, context = NULL) {
  checkmate::assert_class(mobj, "hi_model_spec")
  checkmate::assert_data_frame(data)
  checkmate::assert_character(qVars, null.ok = TRUE)
  checkmate::assert_string(context, null.ok = TRUE)
  lg <- lgr::get_logger("fmri.pipeline/afni_3dlmer/emmeans_to_3dlmer_glt")

  if (is.null(mobj$lmfit)) {
    lg$warn("No lmfit in mobj. Cannot generate emmeans-to-3dlmer GLTs.")
    generated <- empty_3dlmer_glt_table()
  } else {
    generated <- empty_3dlmer_glt_table()
  }

  spec <- mobj$contrast_spec
  if (!is.null(spec) && !is.null(mobj$lmfit)) {
    add_generated <- function(label, code, contrast_type) {
      generated <<- rbind(
        generated,
        data.frame(
          label_raw = label,
          label_afni = label,
          code = code,
          source = "generated",
          contrast_type = contrast_type,
          valid = TRUE,
          message = "",
          stringsAsFactors = FALSE
        )
      )
    }

    reject_interaction <- function(vv, type) {
      if (grepl(":", vv, fixed = TRUE)) {
        stop(
          sprintf(
            "Automatic AFNI 3dLMEr GLT generation supports only one-factor %s. Contrast '%s'%s is an interaction; provide it explicitly with lmer_glt_codes.",
            type,
            vv,
            if (!is.null(context)) paste0(" for ", context) else ""
          ),
          call. = FALSE
        )
      }
    }

    if (length(spec$cond_means) > 0L) {
      for (vv in spec$cond_means) {
        reject_interaction(vv, "condition means")
        if (!vv %in% names(data)) {
          stop(sprintf("Cannot generate AFNI 3dLMEr GLTs for '%s' because it is not in the dataTable.", vv), call. = FALSE)
        }
        levs <- unique(as.character(data[[vv]]))
        levs <- levs[!is.na(levs)]
        for (lev in levs) {
          add_generated(
            label = paste(vv, lev, sep = "."),
            code = paste0(vv, " : 1*", lev),
            contrast_type = "cond_mean"
          )
        }
      }
    }

    if (length(spec$pairwise_diffs) > 0L) {
      for (vv in spec$pairwise_diffs) {
        reject_interaction(vv, "pairwise differences")
        if (!vv %in% names(data)) {
          stop(sprintf("Cannot generate AFNI 3dLMEr GLTs for '%s' because it is not in the dataTable.", vv), call. = FALSE)
        }
        levs <- unique(as.character(data[[vv]]))
        levs <- levs[!is.na(levs)]
        if (length(levs) >= 2L) {
          pairs <- utils::combn(levs, 2L, simplify = FALSE)
          for (pp in pairs) {
            add_generated(
              label = paste0(vv, ".pw.", pp[1L], ".vs.", pp[2L]),
              code = paste0(vv, " : 1*", pp[1L], " -1*", pp[2L]),
              contrast_type = "pairwise_diff"
            )
          }
        }
      }
    }

    if (!is.null(mobj$contrast_list$custom)) {
      cmat <- mobj$contrast_list$custom
      for (i in seq_len(nrow(cmat))) {
        con_name <- rownames(cmat)[i]
        weights <- cmat[i, ]
        non_zero <- which(abs(weights) > 1e-8)
        vars <- names(weights)[non_zero]
        if (!all(vars %in% names(data))) {
          stop(
            sprintf(
              "Cannot convert custom contrast '%s' to AFNI 3dLMEr GLT because regressor(s) are not direct dataTable variables: %s. Provide lmer_glt_codes manually.",
              con_name, paste(setdiff(vars, names(data)), collapse = ", ")
            ),
            call. = FALSE
          )
        }
        code_parts <- paste0(vars, " : ", as.numeric(weights[non_zero]))
        add_generated(con_name, paste(code_parts, collapse = " "), "custom")
      }
    }
  }

  raw_glt_codes <- raw_glt_codes %||% mobj$lmer_glt_codes
  raw <- as_3dlmer_glt_table(raw_glt_codes, source = "user", contrast_type = "raw")
  out <- rbind(generated, raw)
  if (nrow(out) == 0L) {
    return(out)
  }
  out$label_afni <- sanitize_3dlmer_glt_names(out$label_raw)

  validate_3dlmer_glt_table(out, datatable = data, qVars = qVars, context = context)
}
