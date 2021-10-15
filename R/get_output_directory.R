# TODO: clean up the mess of the 'local' approach, but for now, prefer the glue expression
# consider whether we want a 'what' type argument, like subject directory, analysis directory,
# run directory, model-specific analysis directory, etc.

#' small helper function to return the location of an l1 directory based on
#'   id, session, and run number
#'
#' @param id The id of a participant
#' @param session The session number to lookup
#' @param run_number The run number to lookup
#' @param gpa A \code{glm_pipeline_arguments} object
#' @param glm_software which software is being used for the analysis (since directories may vary)
#' @param create_if_missing whether to create the directory if it does not exist
get_output_directory <- function(id = NULL, session = NULL, run_number = NULL,
  l1_model = NULL, l2_model = NULL, l3_model = NULL,
  l1_contrast = NULL, l2_contrast = NULL, l3_contrast = NULL,
  what="l1", gpa, glm_software = "fsl", create_if_missing = FALSE) {

  checkmate::assert_scalar(id, null.ok=TRUE)
  checkmate::assert_integerish(session, null.ok = TRUE)
  if (is.null(session)) session <- 1
  checkmate::assert_integerish(run_number, lower = 1, null.ok = TRUE)
  checkmate::assert_string(l1_model, null.ok=TRUE)
  checkmate::assert_string(l2_model, null.ok = TRUE)
  checkmate::assert_string(l3_model, null.ok = TRUE)
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  stopifnot("run_data" %in% names(gpa))
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_logical(create_if_missing, len=1L)

  if (!is.null(l1_model)) checkmate::assert_subset(l1_model, names(gpa$l1_models$models))
  if (!is.null(l2_model)) checkmate::assert_subset(l2_model, names(gpa$l2_models$models))
  if (!is.null(l3_model)) checkmate::assert_subset(l3_model, names(gpa$l3_models$models))

  lg <- lgr::get_logger("glm_pipeline/l1_setup")

  # create local synonyms for these in case they are used the glue() expression.
  l1_cope_name <- l1_contrast
  l2_cope_name <- l2_contrast
  l3_cope_name <- l3_contrast

  # helper subfunction to select rows where all conditions in filter_list match
  subset_run_data <- function(run_data, filter_list) {
    res <- sapply(seq_along(filter_list), function(xx) {
      run_data[[ names(filter_list)[xx] ]] == filter_list[[xx]]
    })
    all_match <- which(apply(res, 1, all))
    rinfo <- run_data %>% dplyr::slice(all_match)
    if (nrow(rinfo) == 0L) {
      lg$error("Unable to locate a record in gpa$run_data for id %s, session %s, run_number %s.", id, session, run_number)
      return(NULL)
    } else {
      return(rinfo)
    }
  }

  if (what == "l1") {
    if (gpa$output_locations$feat_l1_directory == "local") {
      checkmate::assert_integerish(run_number, lower = 1, null.ok = FALSE) # enforce run_number for l1 lookup
      rinfo <- subset_run_data(gpa$run_data, list(id=id, session=session, run_number=run_number))
      if ("mr_dir" %in% names(rinfo)) {
        lg$debug("Using mr_dir l1 directory lookup: %s", rinfo$mr_dir[1L])
        if (is.null(l1_model)) {
          lg$debug("Lookup from analysis_name: %s", gpa$analysis_name)
          out_dir <- file.path(normalizePath(file.path(rinfo$mr_dir[1L], "..")), gpa$analysis_name)
        } else {
          lg$debug("Lookup from model outdir: %s", gpa$l1_models$models[[l1_model]]$outdir)
          out_dir <- file.path(normalizePath(file.path(rinfo$mr_dir[1L], "..")), gpa$l1_models$models[[l1_model]]$outdir)
        }
      } else {
        # look in parent folder of relevant run nifti and place l1 model there
        rn <- get_mr_abspath(rinfo[1, , drop = F], "run_nifti")
        lg$debug("Using run_nifti l1 directory lookup: %s", rn)

        if (is.null(l1_model)) {
          lg$debug("Lookup from analysis_name: %s", gpa$analysis_name)
          out_dir <- file.path(
            normalizePath(file.path(dirname(rn), "..")),
            gpa$analysis_name
          )
        } else {
          lg$debug("Lookup from model outdir: %s", gpa$l1_models$models[[l1_model]]$outdir)
          out_dir <- file.path(
            normalizePath(file.path(dirname(rn), "..")),
            gpa$l1_models$models[[l1_model]]$outdir
          )
        }
      }
    } else {
      out_dir <- glue::glue(gpa$output_locations$feat_l1_directory) # evaluate glue expression
    }
  } else if (what == "l2") {
    out_dir <- glue::glue(gpa$output_locations$feat_l2_directory) # evaluate glue expression
  } else if (what == "l3") {
    out_dir <- glue::glue(gpa$output_locations$feat_l3_directory) # evaluate glue expression
  } else if (what == "sub") {
    if (gpa$output_locations$feat_sub_directory == "local") {
      rinfo <- subset_run_data(gpa$run_data, list(id = id, session = session, run_number = run_number))
    } else {
      out_dir <- glue::glue(gpa$output_locations$feat_sub_directory) # evaluate glue expression
    }
  } else if (what == "ses") {
    if (gpa$output_locations$feat_ses_directory == "local") {
      rinfo <- subset_run_data(gpa$run_data, list(id = id, session = session))
    } else {
      out_dir <- glue::glue(gpa$output_locations$feat_ses_directory) # evaluate glue expression
    }
  }

  #glue returns empty strings when some variables are NULL
  if (identical(glue("{}"), out_dir) || nchar(out_dir) == 0L) {
    lg$error("Unable to sort out %s output directory", what)
    return(NULL)
  } else {
    #out_dir <- make.names(out_dir) #remove weird characters
  }

  if (isTRUE(create_if_missing) && !dir.exists(out_dir)) {
    lg$debug("Create output directory: %s", out_dir)
    dir.create(out_dir, recursive = TRUE)
  }

  return(out_dir)
}