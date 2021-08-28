
R_batch_job <- R6::R6Class("batch_job",
  private = list(
    batch_generated = FALSE,
    #' write the sbatch, pbs, or local bash script
    write_batch_file = function() {
      syntax <- c("#!/bin/sh")
      
      if (self$scheduler %in% c("slurm", "sbatch")) {
        mem_string <- NULL
        if (!is.null(self$mem_per_cpu)) {
          mem_string <- paste0("#SBATCH --mem-per-cpu=", self$mem_per_cpu)
        } else if (!is.null(self$mem_total)) {
          mem_string <- paste0("#SBATCH --mem-per-cpu=", self$mem_total)
        }

        job_string <- if (is.null(self$job_name)) NULL else paste("#SBATCH -J", self$job_name)
        sched_string <- if (is.null(self$scheduler_options)) NULL else paste("#SBATCH", self$scheduler_options)

        syntax <- c(
          syntax,
          paste("#SBATCH -N", self$n_nodes),
          paste("#SBATCH -n", self$n_cpus),
          paste("#SBATCH -t", self$cpu_time),
          mem_string,
          job_string,
          sched_string,
          "",
          "cd $SLURM_SUBMIT_DIR",
          self$batch_code
        )

      } else if (self$scheduler %in% c("torque", "qsub")) {
        mem_string <- NULL
        if (!is.null(self$mem_per_cpu)) {
          mem_string <- paste0("#PBS -l pmem=", self$mem_per_cpu)
        } else if (!is.null(self$mem_total)) {
          mem_string <- paste0("#PBS -l mem=", self$mem_total)
        }

        syntax <- c(
          syntax,
          paste0("#PBS -l nodes=", self$n_nodes, ":ppn=", self$n_cpus),
          paste("#PBS -l walltime=", self$cpu_time),
          mem_string,
          ifelse(is.null(self$job_name), NULL, paste("#PBS -N", self$job_name)),
          paste("#PBS", self$scheduler_options),
          "",
          "cd $PBS_O_WORKDIR"
        )
      } else {
        #local scheduler considerations here
      }

      syntax <- c(syntax, paste("Rscript --vanilla batch_run.R"))

      message("Writing batch script to: ", file.path(self$batch_directory, "submit_batch.sh"))
      writeLines(syntax, con = file.path(self$batch_directory, "submit_batch.sh"))
    },

    #' write code to be executed
    write_compute_file = function() {
      syntax <- c()
      if (!is.null(self$r_packages)) {
        syntax <- c(
          syntax,
          "if (!require(pacman)) { install.packages('pacman'); library(pacman) }",
          paste0("pacman::p_load(", paste(self$r_packages, collapse = ", "), ")")
        )
      }

      if (!is.null(self$input_environment)) {
        if (!file.exists(self$input_environment)) {
          warning(sprintf("input_environment %s did not exist at the time of script generation.", self$input_environment))
        }
        syntax <- c(syntax, paste0("if (file.exists('", self$input_environment, "')) load('", self$input_environment, "')"))
      }

      syntax <- c(syntax, self$r_code)

      if (!is.null(self$output_environment)) {
        syntax <- c(syntax, paste0("save.image(file='", self$output_environment, "')"))
      }
      message("Writing R script to: ", file.path(self$batch_directory, "batch_run.R"))
      writeLines(syntax, con = file.path(self$batch_directory, "batch_run.R"))
    }
  ),
  public = list(
    parent_jobs = NULL,
    child_jobs = NULL,
    job_id = NULL,
    job_name = NULL,
    cpu_time = "4:00:00", # 4 hours default
    n_nodes = "1",
    n_cpus = "4",
    mem_total = "4g",
    mem_per_cpu = NULL,
    input_environment = "job_environment.RData",
    output_environment = "job_environment.RData",
    sqlite_db = NULL,
    batch_id = NULL,
    batch_directory = "~/", # needs to be somewhere that both compute nodes and login nodes can access, so /tmp is not good.
    batch_code=NULL,
    r_code = NULL,
    r_packages = NULL,
    scheduler = "slurm",
    scheduler_options = NULL,

    #' initialization routine for object
    initialize = function(
      parent_jobs = NULL, job_id=NULL, job_name=NULL, n_nodes = NULL, n_cpus = NULL,
      cpu_time = NULL, mem_per_cpu = NULL, mem_total = NULL, batch_id = NULL,
      r_code = NULL, batch_code, r_packages = NULL, scheduler = NULL, scheduler_options = NULL) {

      if (!is.null(parent_jobs)) self$parent_jobs <- parent_jobs
      if (is.null(job_id)) self$job_id <- uuid::UUIDgenerate()
      if (!is.null(job_name)) self$job_name <- as.character(job_name)
      if (!is.null(n_nodes)) self$n_nodes <- as.character(n_nodes)
      if (!is.null(n_cpus)) self$n_cpus <- as.character(n_cpus)
      if (is.null(cpu_time)) {
        message("Using default cpu_time of: ", self$cpu_time)
      } else {
        self$cpu_time <- as.character(cpu_time)
      }

      if (!is.null(mem_per_cpu)) {
        checkmate::assert_string(mem_per_cpu)
        self$mem_per_cpu <- mem_per_cpu
      }

      if (!is.null(mem_total)) {
        checkmate::assert_string(mem_total)
        self$mem_total <- mem_total
      }

      if (!is.null(batch_id)) self$batch_id <- as.character(batch_id)

      if (is.null(r_code)) {
        stop("Unable to initialize R_batch_job object without r_code")
      } else {
        checkmate::assert_multi_class(r_code, c("expression", "character"))
        self$r_code <- as.character(r_code)
      }

      if (!is.null(batch_code)) {
        checkmate::assert_character(batch_code)
        self$batch_code <- batch_code
      }
      
      if (!is.null(r_packages)) {
        checkmate::assert_character(r_packages)
        self$r_packages <- r_packages
      }

      checkmate::assert_subset(scheduler, c("torque", "qsub", "slurm", "sbatch", "sh", "local"))
      if (!is.null(scheduler)) self$scheduler <- scheduler

      if (!is.null(scheduler_options)) {
        checkmate::assert_character(scheduler_options)
        self$scheduler_options <- scheduler_options
      }
    },

    #' Helper function that generates the batch and compute files for a job
    #' 
    #' @details this is called by \code{$submit} when a job is submitted and is provided
    #'   here in case the user wants to generate the batch files without executing them
    generate = function() {
      private$write_batch_file()
      private$write_compute_file()
      private$batch_generated <- TRUE
    },
    submit = function() {
      if (isFALSE(private$batch_generated)) self$generate()

      batch_script <- file.path(self$batch_directory, "submit_batch.sh")
      self$job_id <- fmri.pipeline::cluster_job_submit(batch_script, scheduler=self$scheduler, wait_jobs=self$parent_jobs)

      # save environment to file
      # write R script that
      # a) loads environment
      # b) load packages
      # c) evaluates expression
      # call cluster_job_submit
    }
  ),
)