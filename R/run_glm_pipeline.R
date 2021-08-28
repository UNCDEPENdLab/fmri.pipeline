

run_glm_pipeline <- function(gpa, l1_models = NULL, l2_models = NULL, l3_models = NULL, glm_software = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_string(l1_models, null.ok = TRUE)
  checkmate::assert_string(l2_models, null.ok = TRUE)
  checkmate::assert_string(l3_models, null.ok = TRUE)

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")
  batch_id <- uuid::UUIDgenerate()
  if (is.null(gpa$batch_run)) gpa$batch_run <- list()


  if (is.null(gpa$finalize_complete) || isFALSE(gpa$finalize_complete)) {
    lg$info("finalize_pipeline_configuration has not been run on this object. We will start with this step.")
    # send this to slurm -- get rslurm jobid

    gpa <- finalize_pipeline_configuration(gpa)
  }

  # solution to old run_glm_pipeline problem
  # submit job that does calls run_featsep etc., then internally does the wait_for_job for all children that are launched. Then make
  # jobs after that dependent on the job that lauches run_featsep.

  # create all FSF files for level one runs
  gpa <- setup_l1_models(gpa)

  #### execution

  jobs <- run_feat_sepjobs(gpa, level = 1L)

  # todo
  # gpa <- verify_lv1_runs(gpa)

  # load("test_gpa.RData")

  # new nomenclature
  # gpa$l1_model_setup$fsl <- gpa$l1_model_setup$fsl %>%
  #   dplyr::rename(l1_feat_fsf=feat_file) %>%
  #   dplyr::mutate(l1_feat_dir=sub(".fsf", ".feat", l1_feat_fsf, fixed=TRUE))

  # save(gpa, file="test_gpa.RData")

  # gpa$parallel$fsl$l2_feat_time <- "1:00:00"
  # gpa$parallel$fsl$l2_feat_memgb <- "20"
  # gpa$parallel$fsl$l3_feat_time <- "1:00:00"
  # gpa$parallel$fsl$l3_feat_memgb <- "20"

  # setup of l2 models (should follow l1)

  # this should be run *after* level 1 runs
  gpa <- setup_l2_models(gpa)

  jobs <- run_feat_sepjobs(gpa, level = 2)

  gpa$parallel$fsl$l3_feat_cpusperjob <- 16

  gpa <- setup_l3_models(gpa)

  jobs <- run_feat_sepjobs(gpa, level = 3)
}

#' small helper function to parse duration syntax of days-hours:minutes:seconds
#'   into lubridate duration object
#' 
#' @param str string containing a duration that may include a days specification
#' @importFrom lubridate hms
#' @keywords internal
dhms <- function(str) {
  checkmate::assert_string(str)
  if (grepl("^\\d+-", str, perl=TRUE)) {
    split_hyphen <- strsplit(str, "-", fixed = TRUE)[[1]]
    days <- as.numeric(split_hyphen[1])
    period <- lubridate::hms(split_hyphen[2:length(split_hyphen)])
    period@day <- days
  } else {
    period <- lubridate::hms(str)
  }
  return(period)
}


#
# conceptual sketch:
#
# initiate_batch(
#  expression(gpa <- finalize_pipeline_configuration(gpa)),
#
#
#
# )
