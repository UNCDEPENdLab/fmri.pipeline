

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

R_batch_job <- R6::R6Class("batch_job",
  private = list(),
  public = list(
    parent_jobs = NULL,
    child_jobs = NULL,
    cpu_time = c("4:00:00"), # 4 hours default
    n_nodes = 1,
    n_cpus = c("4"),
    scheduler = "slurm",
    input_environment = "job_environment.RData",
    output_environment = "job_environment.RData",
    sqlite_db = NULL,
    batch_id = NULL,
    batch_directory = "~/", # needs to be somewhere that both compute nodes and login nodes can access, so /tmp is not good.
    job_id = NULL,
    code_to_run = NULL,
    r_packages = NULL,
    initialize = function(
      parent_jobs = NULL, job_id=NULL, n_nodes = NULL, n_cpus = NULL, 
      cpu_time = NULL, batch_id = NULL, code_to_run = NULL,
      r_packages = NULL, scheduler = NULL, scheduler_args = NULL) {

      if (is.null(job_id)) public$job_id <- uuid::UUIDgenerate()
      if (is.null(cpud_time)) {
        message("Using default cpu_time of: ", public$cpu_time)
      } else {
        public$cpu_time <- as.character(cpu_time)
      }
      
      if (is.null(code_to_run)) {
        stop("Unable to initialize R_batch_job object without code_to_run")
      } else {
        checkmate::assert_multi_class(code_to_run, c("expression", "character"))
        public$code_to_run <- as.character(code_to_run)
      }

      checkmate::assert_subset(scheduler, c("torque", "qsub", "slurm", "sbatch", "sh", "local"))
      if (!is.null(scheduler)) {
        public$scheduler <- scheduler
      }

      # lookup existing records in database
      #

      public$parent_jobs <- parent_jobs
    },

    #' write the sbatch, pbs, or local bash script
    write_batch_file = function() {
      syntax <- c("#!/bin/sh")
      if (public$scheduler %in% c("slurm", "sbatch")) {
        syntax <- c(
          syntax, paste("#SBATCH", public$scheduler_args),
          paste("#SBATCH -n", public$n_cpus), paste("-N", public$n_nodes)
        )
        # SBATCH -p general
        # SBATCH -N 1
        # SBATCH -t 72:00:00
        # SBATCH --mem=6g
        # SBATCH -n 1
      }
    },

    #' write code to be executed
    write_compute_file = function() {
      syntax <- c()
      if (!is.null(public$r_packages)) {
        syntax <- c(
          syntax,
          "if (!require(pacman)) { install.packages('pacman'); library(pacman) }",
          paste0("pacman::p_load(", paste(public$r_packages, collapse = ", "), ")")
        )
      }

      if (!is.null(public$input_environment)) {
        if (!file.exists(public$input_environment)) {
          warning(sprintf("Cannot add load input_environment statement %s to compute script.", public$input_environment))
        } else {
          syntax <- c(syntax, paste0("load(", public$input_environment, ")"))
        }
      }

      syntax <- c(syntax, code_to_run)

      if (!is.null(public$output_environment)) {
        syntax <- c(syntax, paste0("save.image(file=", public$output_environment, ")"))
      }
      writeLines(syntax, con = )
    },
    submit = function() {
      if (public$scheduler %in% c("slurm", "sbatch")) {

      } else if (public$scheduler %in% c("torque", "qsub")) {

      } else if (public$schedule %in% c("local", "sh")) {

      }


      # save environment to file
      # write R script that
      # a) loads environment
      # b) load packages
      # c) evaluates expression
      # call cluster_submit_job
    }
  ),
)

#
# conceptual sketch:
#
# initiate_batch(
#  expression(gpa <- finalize_pipeline_configuration(gpa)),
#
#
#
# )
