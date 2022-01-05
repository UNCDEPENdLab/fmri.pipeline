create_fwe_spec <- function(gpa, level = level, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1, upper = 3)
  checkmate::assert_class(lg, "Logger")

  cat("Welcome to the FWE setup menu\n")
  cat("At present, we offer 4 forms of FWE correction: pTFCE, AFNI 3dFWHMx + 3dClustSim, 3dttest++ residual permutation + 3dClustSim, and PALM permutation\n")

  which_fwe <- select.list(
    choices = c("pTFCE", "3dFWHMx + 3dClustSim", "3dtest++ -randomsign + 3dClustSim", "PALM"),
    title = c("Which correction would you like to add?")
  )

  stopifnot(level == 3L) # only support 3rd-level FWE for now

  # need to choose to which models this applies...
  # for now, skip the model lookup from the run step and just rely on l3_model_setup
  # if (level == 3L) {
  #   l3_set <- choose_glm_models(gpa, "prompt", level = 3L, lg = lg)
  #   if (is.null(l3_set)) {
  #     message("No l3 models were selected for this FWE.")
  #     return(NULL) #maybe it's okay not to have any models?
  #   }
  # }

  fwe_set <- list()
  fields <- c("l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model")
  for (ff in fields) {
    fwe_set[[ff]] <- select.list(unique(gpa$l3_model_setup$fsl[[ff]]),
      multiple = TRUE,
      title = paste("For this FWE, choose which", ff, "to include in this FWE.")
    )
  }

  # build filter statement
  filter_expr <- paste(lapply(seq_along(fwe_set), function(x) {
    paste(names(fwe_set)[x], "%in%", paste(deparse(fwe_set[[x]]), collapse = ""))
  }), collapse = " & ")

  # https://www.r-bloggers.com/2020/05/filtering-with-string-statements-in-dplyr/
  to_fwe <- gpa$l3_model_setup$fsl %>% dplyr::filter(rlang::eval_tidy(rlang::parse_expr(filter_expr)))
  
  if (nrow(to_fwe) == 0L) {
    msg <- "No model outputs match this combination of selections."
    lg$error(msg)
    stop(msg)
  }

  if (any(to_fwe$feat_complete == FALSE)) {
    lg$warn("Some L3 models requested are currently incompleted, preventing FWE correction.")
    lg$warn("%s", capture.output(print(to_fwe %>% filter(feat_complete == FALSE))))
    to_fwe <- to_fwe %>% dplyr::filter(feat_complete == TRUE)
    if (nrow(to_fwe) == 0L) {
      msg <-"All L3 models requested are incomplete and this function cannot continue."
      lg$error(msg)
      stop(msg)
    }
  }


  # I guess we could have lower level subsetting, too.
  # gfeat_set <- gpa$l3_model_setup$fsl %>% dplyr::filter(l3_model == !!l3_set & feat_complete == TRUE)
  
  #choices = c("pTFCE", "3dFWHMx + 3dClustSim", "3dtest++ -randomsign + 3dClustSim", "PALM"),
  if (which_fwe == "pTFCE") {
    ptfce_objs <- lapply(to_fwe$feat_dir, function(gg) {
      oo <- ptfce_spec$new(gfeat_dir = gg)
    })
  }

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