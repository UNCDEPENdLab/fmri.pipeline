create_fwe_spec <- function(gpa, level = level) {
  checkmate::assert_integerish(level, lower=1, upper=3)
  cat("Welcome to the FWE setup menu\n")
  cat("At present, we offer 3 forms of FWE correction: pTFCE, AFNI 3dFWHMx + 3dClustSim, and PALM permutation\n")

  which_fwe <- menu(c("pTFCE", "3dFWHMx + 3dClustSim", "PALM"), title = c("Which correction would you like to add?"))

}

#' R6 class for an FWE correction method that can apply to one or more models
#'
#' @importFrom R6 R6Class
#' @importFrom data.table data.table
#' @importFrom gtools permutations
#' @importFrom checkmate assert_data_frame assert_logical
#' @keywords internal
fwe_spec <- R6::R6Class("fwe_spec",
  private = list(
    # data: keyed data.table object
    data = NULL,
    fwe_type = NULL
  ), public = list(
    initialize = function(fwe_type = NULL) {
      checkmate::assert_subset(fwe_type, c("ptfce, 3dclustsim", "palm", "randomise"), empty.ok = FALSE)
    }
  )
)



build_fwe_correction <- function(gpa, lg=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_class(lg, "Logger")

  level <- 3L # for now, only supporting FWE at group analysis level

  if (!checkmate::test_class(gpa$l3_models, "hi_model_set")) {
    msg <- "Cannot setup FWE correction until build_l3_models is complete and $l3_models is populated."
    lg$error(msg)
    stop(msg)
  }

  # to setup an FWE specification, we only need the models in place, but they need not be complete.
  # cat(glue("\n\n---\nCurrent {field_desc} columns: {c_string(chosen_cols)}\n\n", .trim=FALSE))
  
  action <- 0L
  while (action != 4L) {
    action <- menu(c("Add FWE correction", "Modify FWE correction", "Delete FWE correction", "Done with FWE setup"),
      title = "What would you like to like to do?"
    )

    if (action == 1L) {
      gpa <- create_fwe_spec(gpa, level = level)

    } else if (action == 2L) {

    } else if (action == 3L) {

    } else if (action == 4L) {
      cat("Finishing FWE correction setup\n")
    }
  }
  


  if ("fsl" %in% gpa$glm_software) {
    # enforce that at least some L3 models are complete
    enforce_glms_complete(gpa, level = 3L, lg)

    #complete_l3 <- #gpa$
  } else if ("spm" %in% gpa$glm_software) {
    stop("Not supported yet")
  } else if ("afni" %in% gpa$glm_software) {
    stop("Not supported yet")
  }

  
}