#' Description of R_batch_job R6 class
#' @export
R_batch_job <- R6::R6Class("batch_job",
  private = list(
    # whether the batch and compute files have already been created
    batch_generated = FALSE,

    # absolute path to batch file
    batch_file_name = NULL,

    # absolute path to compute file
    compute_file_name = NULL,

    # field job_id the process or scheduler job id that uniquely identifies this job
    job_id = NULL,
    get_unique_file_name = function(script_name = "") {
      base <- tools::file_path_sans_ext(script_name)
      ext <- tools::file_ext(script_name)

      script_name <- paste0(
        base,
        ifelse(is.null(self$job_name), "", paste0("_", self$job_name)),
        ".", ext
      )

      if (file.exists(normalizePath(file.path(self$batch_directory, script_name), mustWork = FALSE))) {
        # add unique random string to script name if the file already exists
        script_name <- sub(
          paste0("\\.", ext, "$"),
          paste0("_", sub("^/", "", tempfile(pattern = "", tmpdir = "")), ".", ext),
          script_name
        )
      }

      return(normalizePath(file.path(self$batch_directory, script_name), mustWork = FALSE))
    },

    # helper to generate and return the full path to the batch script
    get_batch_file_name = function() {
      if (is.null(private$batch_file_name)) {
        private$batch_file_name <- private$get_unique_file_name("submit_batch.sh")
      }
      return(private$batch_file_name)
    },

    # helper to generate and return the full path to the compute script
    get_compute_file_name = function() {
      if (is.null(private$compute_file_name)) {
        private$compute_file_name <- private$get_unique_file_name("batch_run.R")
      }
      return(private$compute_file_name)
    },

    # write the sbatch, pbs, or local bash script
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
        # local scheduler considerations here
      }

      # syntax <- c(syntax, paste("Rscript --vanilla", private$get_compute_file_name()))
      syntax <- c(syntax, paste("R CMD BATCH --no-save --no-restore", private$get_compute_file_name()))

      message("Writing batch script to: ", private$get_batch_file_name())
      writeLines(syntax, con = private$get_batch_file_name())
    },

    # write code to be executed
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
          # warning(sprintf("input_environment %s did not exist at the time of script generation.", self$input_environment))
        }
        syntax <- c(syntax, paste0("if (file.exists('", self$input_environment, "')) load('", self$input_environment, "')"))
      }

      syntax <- c(syntax, self$r_code)

      if (!is.null(self$output_environment)) {
        syntax <- c(syntax, paste0("save.image(file='", self$output_environment, "')"))
      }
      message("Writing R script to: ", private$get_compute_file_name())
      writeLines(syntax, con = private$get_compute_file_name())
    }
  ),
  public = list(
    #' @field parent_jobs a vector of parent_jobs that are upstream of this job and may influence its execution
    parent_jobs = NULL,

    #' @field depends_on_parents logical or character string indicating whether this job should wait until
    #'    \code{parent_jobs} complete. If a character string is passed, it indicates which named elements of
    #'    \code{parent_jobs} are must complete before this job begins.
    depends_on_parents = FALSE,

    #' @field child_jobs a vector of job ids generated by the submission of this job. Only relevant if this job
    #'    submits additional jobs through the scheduler directly
    child_jobs = NULL,

    #' @field job_name a user-defined name for the job used for specifying job dependencies and informative job
    #'    status queries on a job scheduler
    job_name = NULL,

    #' @field cpu_time The amount of time requested on the job scheduler, following d-hh:mm:ss format. Defaults to
    #'    "4:00:00", which is 4 hours.
    cpu_time = "4:00:00",

    #' @field n_nodes The number of nodes to be requested on the job scheduler
    n_nodes = "1",

    #' @field n_cpus The number of cores (aka 'cpus', ignoring hyperthreading) to be requested on the job scheduler
    n_cpus = "4",

    #' @field mem_total The total amount of memory (RAM) requested by the job
    mem_total = "4g",

    #' @field mem_per_cpu The amount of memory (RAM) requested per cpu (total = mem_per_cpu * n_cpus)
    mem_per_cpu = NULL,

    #' @field input_environment The name of the environment to be loaded at the beginning of the R batch prior to executing
    #     other code. Used to setup any local objects needed to begin computation.
    input_environment = "job_environment.RData",

    #' @field output_environment The name of the environment to be saved at the end of the R batch execution, which can then
    #'    be loaded by subsequent jobs.
    output_environment = "job_environment.RData",

    #' @field sqlite_db Not used yet, but will be used for job tracking in future
    sqlite_db = NULL,

    #' @field batch_id Not currently used, but intended for job sequence tracking
    batch_id = NULL,

    #' @field batch_directory Location of batch scripts to be written
    batch_directory = "~/", # needs to be somewhere that both compute nodes and login nodes can access, so /tmp is not good.

    #' @field batch_code Shell code to be included in the batch script prior to the R code to be run. This can include
    #'    module load statements, environment variable exports, etc.
    batch_code = NULL,

    #' @field r_code The R code to be executed by the scheduler. This can be a character vector that includes multiple
    #'    R statements or an expression object containing the R code to be evaluated
    r_code = NULL,

    #' @field r_packages The R packages to be loaded into the environment before job execution. These are loaded by
    #'    pacman::p_load, which will install any missing packages before attempting to load
    r_packages = NULL,

    #' @field scheduler The job scheduler to be used for this batch. Options are: "slurm", "torque", or "local".
    scheduler = "slurm",

    #' @field scheduler_options An optional character vector of scheduler arguments to be included in the batch script header
    #'    that control additional features such as job emails or group permissions. These directives are added with #SBATCH
    #'    or #PBS headings, depending on the scheduler, and are ignored if the scheduler is "local".
    scheduler_options = NULL,

    #' @field repolling_interval The number of seconds to wait between successive checks on whether parent jobs have completed.
    #'    This is mostly relevant to the 'local' scheduler.
    repolling_interval = 60, # seconds

    #' @description Create a new R_batch_job object for execution on an HPC cluster
    #' @param parent_jobs A vector of one or more job ids that are parents to this job. This can be a named vector, to
    #'    be used in conjunction with \code{depends_on_parents} to specify which parent jobs must be completed before this
    #'    job begins.
    #' @param job_name The name of this job
    initialize = function(batch_directory = NULL, parent_jobs = NULL, job_name = NULL, n_nodes = NULL, n_cpus = NULL,
                          cpu_time = NULL, mem_per_cpu = NULL, mem_total = NULL, batch_id = NULL, r_code = NULL,
                          batch_code = NULL, r_packages = NULL, scheduler = NULL, scheduler_options = NULL, repolling_interval = NULL) {
      if (!is.null(batch_directory)) self$batch_directory <- batch_directory
      if (!is.null(parent_jobs)) self$parent_jobs <- parent_jobs
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

      if (!is.null(repolling_interval)) {
        checkmate::assert_number(repolling_interval, lower = 0.1, upper = 2e5)
        self$repolling_interval <- repolling_interval
      }
    },

    #' @description Helper function that generates the batch and compute files for a job
    #'
    #' @details this is called by \code{$submit} when a job is submitted and is provided
    #'   here in case the user wants to generate the batch files without executing them
    generate = function() {
      # create batch_directory, if missing
      if (!dir.exists(self$batch_directory)) dir.create(self$batch_directory, recursive = TRUE)
      private$write_batch_file()
      private$write_compute_file()
      private$batch_generated <- TRUE
    },

    #' @description Submit job to scheduler or local compute
    submit = function() {
      if (isFALSE(private$batch_generated)) self$generate()
      cd <- getwd()
      setwd(self$batch_directory)

      # For local compute of an R job, use Rscript --vanilla to get proper job_id in return,
      # rather than job_id of shell script that called Rscript.
      # if (scheduler %in% c("local", "sh")) {
      #  batch_script <- private$get_compute_file_name()
      # } else {
      batch_script <- private$get_batch_file_name()
      # }

      # submit job for computation, waiting for parents if requested
      wait_jobs <- NULL
      if (checkmate::test_logical(self$depends_on_parents, max.len = 1) && isTRUE(self$depends_on_parents)) {
        # wait on any/all parent jobs
        wait_jobs <- self$parent_jobs
      } else if (checkmate::test_character(self$depends_on_parents)) {
        # enforce that all parent job dependencies need to be named elements of $parent_jobs
        if (checkmate::test_subset(self$depends_on_parents, names(self$parent_jobs))) {
          wait_jobs <- self$parent_jobs[self$depends_on_parents]
        } else {
          # could make this more informative
          warning("Could not add parent job dependency because one or more depend_on_parents elements were not in parent_jobs")
        }
      }

      # TODO: if a job_id already exists and submit is called again, do we check job status, insist a 'forced' submission?
      private$job_id <- fmri.pipeline::cluster_job_submit(batch_script,
        scheduler = self$scheduler,
        wait_jobs = wait_jobs, repolling_interval = self$repolling_interval
      )

      if (!is.null(cd) && dir.exists(cd)) setwd(cd) # reset working directory (don't attempt if that directory is absent)
    },

    #' @description Function to create a deep copy of a batch job
    #' @details Note that this also resets the compute_file_name and batch_file_name fields so that the
    #'    copied object doesn't create files that collide with the original
    copy = function() {
      cloned <- self$clone(deep = TRUE)
      cloned$reset_file_names()
      return(cloned)
    },

    #' @description helper function to reset names of compute and batch files that will be generated by this job.
    #' @details This needs to be exposed as a public method for copied objects to be able to reset these private fields.
    reset_file_names = function() {
      private$batch_generated <- FALSE
      private$batch_file_name <- NULL
      private$compute_file_name <- NULL
    },

    #' @description Return the job id of this job (populated by job submission)
    get_job_id = function() {
      private$job_id
    }
  ),
)

#' Description of R_batch_sequence R6 class
#' @export
R_batch_sequence <- R6::R6Class("batch_sequence",
  private = list(
    # list of R_batch_job classes to be run sequentially
    sequence_jobs = list()
  ),
  public = list(
    #' @description create a new R_batch_sequence object
    initialize = function(..., joblist = NULL) {
      others <- list(...)
      if (!is.null(joblist) && is.list(joblist)) {
        sapply(joblist, checkmate::assert_class, "batch_job")
        private$sequence_jobs <- joblist
      } else {
        private$sequence_jobs <- others
      }
    },

    #' @description submit the job sequence to the scheduler or local compute
    submit = function() {
      job_list <- private$sequence_jobs
      njobs <- length(job_list)
      for (ii in seq_len(njobs)) {
        this_job <- job_list[[ii]]
        this_job$submit()
        if (ii < njobs) {
          dependent_children <- sapply(seq_along(job_list), function(jnum) {
            this_job$job_name %in% job_list[[jnum]]$depends_on_parents && jnum > ii
          })

          if (any(dependent_children)) {
            job_list[dependent_children] <- lapply(job_list[dependent_children], function(job) {
              # Always add parent job id to any downstream jobs in sequence that depend on this job.
              job$parent_jobs[this_job$job_name] <- this_job$get_job_id()
              return(job)
            })
          }
        }
      }
    }
  )
)
