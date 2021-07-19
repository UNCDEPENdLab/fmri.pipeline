#' This is a small helper function that returns the inputs provided in the feat_files field for
#' a set of .fsf files.
#' 
#' @details 
#' One can also pass in .feat or .gfeat directories and the function will
#' use the design.fsf files within each of these to find the inputs
#' 
#' @param feat_list a vector of design.fsf filenames, and/or .feat/.gfeat directory names
#' @param recursive a boolean indicating whether to drill down and find lower-level inputs for .gfeat/.feat inputs
#' @author Michael Hallquist
#' @importFrom checkmate assert_file_exists assert_directory_exists
#' @export
read_feat_inputs <- function(feat_list, recursive=FALSE) {
  is_feat_dir <- grepl(".g?feat$", feat_list, perl=TRUE)

  if (any(!is_feat_dir)) {
    check_fsfs <- feat_list[!is_feat_dir]
    sapply(check_fsfs, checkmate::assert_file_exists)
    if (!all(goodfiles <- grepl(".fsf$", check_fsfs, perl=TRUE))) {
      error("The following inputs do not end in .fsf, .gfeat, or .feat: ", paste(check_fsfs[!goodfiles], collapse=", "))
    }
  }

  #add the /design.fsf suffix to directory inputs
  if (any(is_feat_dir)) { 
    sapply(feat_list[is_feat_dir], checkmate::assert_directory_exists)
    feat_list[is_feat_dir] <- file.path(feat_list[is_feat_dir], "design.fsf") 
  }

  inputs <- lapply(feat_list, function(ff) {
    fsf <- readLines(ff)
    infile <- grep("^set feat_files\\(\\d+\\)", fsf, perl=TRUE, value=TRUE)
    #input <- paste0(sub("set feat_files\\(d+\\) \"([^\"]+)\"", "\\1", nifti, perl=TRUE), ".nii.gz")
    fnumbers <- as.numeric(sub("set feat_files\\((\\d+)\\) \"[^\"]+\"", "\\1", infile, perl=TRUE)) #always return in sorted order
    infile <- sub("set feat_files\\(\\d+\\) \"([^\"]+)\"", "\\1", infile, perl=TRUE)
    return(infile[fnumbers])
  })

  inputs_are_featdirs <- grepl(".g?feat$", inputs, perl=TRUE)
  if (isTRUE(recursive) && any(inputs_are_featdirs)) {
    inputs <- lapply(1:length(inputs), function(ii) {
      if (isTRUE(inputs_are_featdirs[ii])) {
        return(read_feat_inputs(inputs[ii], recursive=recursive))
      } else {
        return(list(inputs[ii])) #single-element list
      }
    })
  }

  #need some sort of wrap-up function here in the recursive case
  return(inputs)
}


#' helper function to look at whether feat ingredients exist for a .feat/.gfeat directory
#' 
#' @param feat_dir an expected .feat/.gfeat directory for an analysis
#' @param feat_fsf an expected .fsf file corresponding to the feat analysis
#' @param lg an optional lgr logger object to be used for logging. If not passed, the root
#'   logger will be used
#'
#' @return a data.frame containing information about whether the .feat analysis is complete
#'   and whether various ingredients are present
#'
#' @importFrom anytime anytime
#' @keywords internal
get_feat_status <- function(feat_dir, feat_fsf=NULL, lg=NULL, prefix=NULL) {
  checkmate::assert_string(feat_dir)
  checkmate::assert_string(feat_fsf, null.ok = TRUE)
  if (is.null(lg)) lg <- lgr::get_logger() #just use root logger

  # handle feat_dir with no extension
  if (!grepl(".feat$", feat_dir)) {
    feat_dir <- paste0(feat_dir, ".feat")
    if (!dir.exists(feat_dir) && dir.exists(sub(".feat$", ".gfeat", feat_dir))) {
      feat_dir <- sub(".feat$", ".gfeat", feat_dir)
    }
  }

  if (is.null(feat_fsf)) {
    # assume that fsf of same name as feat_dir exists at same level of filesystem hierarchy
    feat_fsf <- sub("\\.g?feat$", ".fsf", feat_dir)
  }

  if (dir.exists(feat_dir) && !file.exists(feat_fsf) && file.exists(file.path(feat_dir, "design.fsf"))) {
    # fall back to design.fsf inside the .feat folder if .fsf is not present in parent folder
    feat_fsf <- file.path(feat_dir, "design.fsf")
  }

  feat_checks <- list()
  feat_checks$feat_fsf <- feat_fsf
  feat_checks$feat_fsf_modified_date <- if (file.exists(feat_fsf)) file.info(feat_fsf)$mtime else as.POSIXct(NA)
  feat_checks$feat_fsf_exists <- file.exists(feat_fsf)
  feat_checks$feat_dir <- feat_dir
  feat_checks$feat_dir_exists <- dir.exists(feat_dir)
  feat_checks$feat_execution_start <- as.POSIXct(NA)
  feat_checks$feat_execution_end <- as.POSIXct(NA)
  feat_checks$feat_execution_min <- as.POSIXct(NA)
  feat_checks$feat_complete <- FALSE # default: FALSE for anything but a .feat_complete outcome
  feat_checks$feat_failed <- NA # default: NA if feat hasn't even been run


  if (dir.exists(feat_dir)) {
    if (file.exists(file.path(feat_dir, ".feat_complete"))) {
      lg$debug("Feat directory is complete: %s", feat_dir)
      timing_file <- file.path(feat_dir, ".feat_complete")
      if (file.exists(file.path(feat_dir, ".feat_fail"))) {
        lg$warn("Both .feat_complete and .feat_fail objects exist in %s", feat_dir)
        lg$warn("Assuming that .feat_complete reflects a successful completion of feat")
      }
      feat_checks$feat_complete <- TRUE
      feat_checks$feat_failed <- FALSE
    } else if (file.exists(file.path(feat_dir, ".feat_fail"))) {
      lg$debug("Detected feat failure in: %s", feat_dir)
      timing_file <- file.path(feat_dir, ".feat_fail")
      feat_checks$feat_failed <- TRUE
    } else {
      timing_file <- NULL
    }

    if (!is.null(timing_file)) {
      timing <- readLines(timing_file)
       if (length(timing) > 0L) {
         # convert to POSIXct object to allow for any date calculations
         timing <- anytime::anytime(timing)
         feat_checks$feat_execution_start <- timing[1L]
         if (length(timing) == 2L) {
           feat_checks$feat_execution_end <- timing[2L]
           feat_checks$feat_execution_min <- as.numeric(difftime(timing[2L], timing[1L], units = "mins"))
         } else {
           lg$warn("Did not find two timing entries in %s.", timing_file)
           lg$warn("File contents: %s", timing)
         }
       }
    }
  }

  df <- as.data.frame(feat_checks)
  if (!is.null(prefix)) { #support naming prefix like "l1_"
    names(df) <- paste0(prefix, names(df))
  }
  return(df)
}

#' Helper function to find/replace or insert additional arguments to include
#'   in FSF syntax for a FEAT model.
#'
#' @param fsf_syntax a character vector containing all fsf syntax
#' @param feat_args a named list 
#' @return a modified version of \code{fsf_syntax} that now includes all custom
#'   arguments from \code{feat_args}
#'
#' @details If a tag already exists in the fsf syntax, it is replaced by the custom argument
#'   in \code{feat_args}. If the tag does not exist, it is inserted.
#' @keywords internal
add_custom_feat_syntax <- function(fsf_syntax, feat_args, lg=NULL) {

  checkmate::assert_character(fsf_syntax)
  checkmate::assert_list(feat_args, null.ok = TRUE)

  if (is.null(feat_args)) {
    #no custom arguments, return fsf_syntax unchanged
    return(fsf_syntax)
  }

  if (is.null(lg)) {
    lg <- lgr::get_logger() #use root logger if not passed
  }

  # for now, this only supports single name/value pairs, even though the list may contain
  # things like list(arg=c(1,2,3,4)). The function skips these for now.
  # TODO: Provide better enforcement of expected structure in finalize_pipeline_configuration
  # or setup_glm_pipeline so that user rectifies things at setup.

  # handle additional custom feat level 1 fields in fsf syntax

  for (ii in seq_along(feat_args)) {
    this_name <- names(feat_args)[ii]
    this_value <- feat_args[[ii]]
    if (length(this_value) > 1L) {
      lg$warn("feat_args: Cannot handle settings with multiple values: %s. Skipping out.", this_name)
      next
    } else {
      lg$debug("Adding custom feat setting: %s = %s", this_name, this_value)
      if (any(grepl(paste0("set fmri(", this_name, ")"), fsf_syntax, fixed = TRUE))) {
        lg$debug("Substituting existing value of feat setting: %s", this_name)
        fsf_syntax <- gsub(paste0("(set fmri\\s*\\(", this_name, "\\))\\s*(.*)"), paste0("\\1 ", this_value),
          fsf_syntax,
          perl = TRUE
        )
      } else {
        fsf_syntax <- c(fsf_syntax, paste0("set fmri(", this_name, ") ", this_value))
      }
    }
  }

  return(fsf_syntax)

}

#' small helper function to look at feat outputs and files to determine if execution is complete
#'
#' @param gpa a \code{glm_pipeline_arguments object}
#' @param level the level of analysis to be refreshed (1, 2, or 3)
#' @return a modified copy of \code{gpa} with the feat columns of
#'   $l1_model_setup, $l2_model_setup, or $l3_model_setup refreshed
#'
#' @keywords internal
#' @importFrom dplyr select
#' @importFrom purrr pmap_dfr
refresh_feat_status <- function(gpa, level = 1L, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)
  setup_name <- paste0("l", level, "_model_setup")

  if ("fsl" %in% gpa$glm_software && setup_name %in% names(gpa)) {
    if (is.null(lg)) lg <- lgr::get_logger()
    lg$info("Found existing %s field. Refreshing status of L%d feat execution and outputs.", setup_name, level)
    refresh <- gpa[[setup_name]]$fsl %>%
      dplyr::select(feat_dir, feat_fsf) %>%
      purrr::pmap_dfr(get_feat_status, lg = lg)

    # copy back relevant columns into data structure
    gpa[[setup_name]]$fsl[, names(refresh)] <- refresh
  }

  return(gpa)
}