#' R6 class for running 3dFWHMx on a single input file based on user specification
#'
#' @importFrom R6 R6Class
#' @details N.B. This class doesn't even expose the Gaussian ACF options given false positive problems
#' @export
fwhmx <- R6::R6Class("fwhmx",
  private = list(
    input_file = NULL,
    mask_file = NULL,
    out_dir = NULL, # where to put calculations
    out_summary = NULL,
    out_detailed = NULL,
    out_by_subbrik = NULL,
    demed = FALSE,
    unif = FALSE,
    average = "-geom",
    call = NULL,
    acf_params = rep(NA_real_, 4), # should always be 4 in length
    ret_code = NULL,
    fwhmx_complete = FALSE,
    build_call = function() {
      if (is.null(private$mask_file)) {
        mask_string <- "-automask"
      } else {
        mask_string <- glue::glue("-mask {private$mask_file}")
      }
      private$call <- glue::glue(
        "3dFWHMx -overwrite -acf {private$out_detailed}",
        " -out {private$out_by_subbrik}",
        " -input {private$input_file} {mask_string} {private$average}",
        " > {private$out_summary}"
      )
    },
    populate_params = function(force=FALSE) {
      if (isTRUE(private$fwhmx_complete) && isFALSE(force)) {
        # don't re-read summary object unless requested explicitly
        return(invisible(NULL))
      }

      if (checkmate::test_file_exists(private$out_summary)) {
        rr <- readLines(private$out_summary)
        if (length(rr) == 2L) {
          private$acf_params <- as.numeric(strsplit(trimws(rr[2L]), "\\s+")[[1L]])
          private$fwhmx_complete <- TRUE
        } else {
          warning("Cannot parse ACF params file: ", private$out_summary)
          print(rr)
        }
      } else {
        warning("Cannot locate ACF summary statistics in file: ", private$out_summary)
      }

      # always return self for side-effect methods
      return(invisible(self))
    }
  ),
  public = list(
    initialize = function(input_file = NULL, mask_file = NULL, out_dir = NULL, demed = NULL, unif = NULL, average = "geometric") {
      if (is.null(input_file)) {
        stop("Cannot run 3dFWHMx without -input dataset.")
      } else {
        checkmate::assert_file_exists(input_file)
        private$input_file <- normalizePath(input_file)
      }

      if (is.null(mask_file)) {
        message("No mask_file provided. Will default to 3dFWHMx -automask")
      } else {
        checkmate::assert_file_exists(mask_file)
        private$mask_file <- mask_file
      }

      if (is.null(out_dir)) {
        private$out_dir <- dirname(private$input_file)
      } else {
        private$out_dir <- normalizePath(out_dir)
      }

      private$out_summary <- file.path(private$out_dir, "3dFWHMx_acf.txt")
      private$out_by_subbrik <- file.path(private$out_dir, "3dFWHMx_acf_bysubbrik.txt")
      private$out_detailed <- file.path(private$out_dir, "3dFWHMx_acf_radius.txt")

      if (file.exists(private$out_summary)) {
        # populate from cached txt file (rather than forcing the command to re-run)
        private$populate_params()
      }

      # demedian
      if (!is.null(demed)) {
        checkmate::assert_logical(demed, len = 1L)
        private$demed <- demed
      }

      # normalize sub-briks to the same MAD
      if (!is.null(unif)) {
        checkmate::assert_logical(unif, len = 1L)
        private$unif <- unif
      }

      if (!is.null(average)) {
        checkmate::assert_string(average)
        checkmate::assert_subset(average, c("arithmetic", "geometric"))
        if (average == "geometric") {
          private$average <- "-geom"
        } else if (average == "arithmetic") {
          private$average <- "-arith"
        }
      }
    },
    run = function(force=FALSE) {
      if (!checkmate::test_directory_exists(private$out_dir)) {
        message("Creating output directory: ", private$out_dir)
        dir.create(private$out_dir, showWarnings = FALSE, recursive = TRUE)
      }

      if (isTRUE(private$fwhmx_complete) && isFALSE(force)) {
        message("3dFWHMx has already completed for this input.")
        return(private$ret_code)
      }

      private$build_call()
      private$ret_code <- run_afni_command(private$call)
      if (private$ret_code != 0) {
        warning("3dFWHMx call returned error code.")
      } else {

      }
    },
    get_call = function() {
      private$build_call()
      return(private$call)
    },
    get_acf_params = function() {
      private$populate_params()
      if (isFALSE(private$fwhmx_complete)) {
        warning("3dFWHMx has not run successfully for this subject. Please use the $run() method first!")
      }
      private$acf_params
    },
    get_input_file = function() {
      private$input_file
    },
    get_outputs = function() {
      c(
        out_dir = private$out_dir,
        out_summary = private$out_summary,
        out_detailed = private$out_detailed,
        out_by_subbrik = private$out_by_subbrik
      )
    },
    is_fwhmx_complete = function() {
      private$fwhmx_complete
    }
  )
)

#' R6 class for running 3dFWHMx on a group of input files using a scheduler/cluster
#'
#' @importFrom R6 R6Class
#' @export
fwhmx_set_spec <- R6::R6Class("fwhmx_set_spec",
  private = list(
    fwhmx_batch = NULL, # batch object
    walltime = "10:00:00", # 10 hours default for all jobs
    fwhmx_objs = list(),
    fwhmx_df = data.frame(),
    fwhmx_acf_avg = setNames(rep(NA_real_, 4), c("a", "b", "c")),
    fwhmx_effective_fwhm = NA_real_,
    all_fwhmx_complete = FALSE,
    scheduler = "slurm",
    input_files = NULL,
    mask_files = NULL,

    populate_acf_df = function() {
      res <- lapply(private$fwhmx_objs, function(x) x$get_acf_params())
      inp_files <- sapply(private$fwhmx_objs, function(x) x$get_input_file())
      lens <- sapply(res, length)
      if (!all(lens == 4)) {
        stop("The length of some ACF outputs from $get_acf_params is not 4. Don't know how to proceed!")
      }
      acf_mat <- do.call(rbind, res)
      colnames(acf_mat) <- c("a", "b", "c", "effective_fwhm")

      private$fwhmx_df <- data.frame(input = inp_files, acf_mat)

      # always return self for side-effect methods
      return(invisible(self))
    },
    get_fwhmx_calls = function(include_complete = FALSE) {
      if (isTRUE(include_complete)) {
        fwhmx_return <- rep(TRUE, length(private$fwhmx_objs))
      } else {
        fwhmx_return <- sapply(private$fwhmx_objs, function(x) !x$is_fwhmx_complete())
      }

      return(sapply(private$fwhmx_objs[fwhmx_return], function(x) x$get_call()))
    },
    generate_batch = function(include_complete = FALSE) {
      # this will return all 3dfwhmx calls

      fwhmx_calls <- private$get_fwhmx_calls(include_complete = include_complete)
      if (length(fwhmx_calls) == 0L) {
        message("All 3dFWHMx runs have already completed. If you want to force a re-run, use $run(force=TRUE)")
        private$fwhmx_batch <- NULL # reset to NULL for consistency
        return(invisible(NULL))
      }

      # create parent batch job to run all 3dFWHMx scripts
      private$fwhmx_batch <- R_batch_job$new(
        job_name = "run_3dfwhmx", n_cpus = 1,
        cpu_time = private$walltime, scheduler = private$scheduler,
        input_objects = list(fwhmx_calls = fwhmx_calls), # export this object to the job
        wait_for_children = TRUE, r_packages = "fmri.pipeline",
        r_code = glue(
          "child_job_ids <- cluster_submit_shell_jobs(fwhmx_calls, memgb_per_command=8, fork_jobs=TRUE, scheduler='{private$scheduler}')"
        )
        # cleanup_r_code = glue(
        #   "batch_obj$refresh()" #update fwhmx objs/params, 
        # )
      )
    }
  ),
  public = list(
    initialize = function(input_files = NULL, mask_files = NULL, scheduler = NULL, walltime="10:00:00", ...) {
      if (is.null(input_files)) {
        stop("fwhmx_set_spec requires at least one input file")
      }

      checkmate::assert_character(input_files)
      checkmate::assert_file_exists(input_files)
      private$input_files <- input_files

      if (!is.null(mask_files)) {
        checkmate::assert_file_exists(mask_files)
        stopifnot(length(mask_files) == length(input_files)) # enforce length match for inputs and masks (otherwise, how do they line up?)
        private$mask_files <- mask_files
      }

      self$refresh() # populates fwhmx_objs and all_fwhmx_complete based on inputs

      if (!is.null(scheduler)) {
        checkmate::assert_string(scheduler)
        checkmate::assert_subset(scheduler, c("torque", "slurm", "local"))
        private$scheduler <- scheduler
      }

      if (!is.null(walltime)) {
        private$walltime <- validate_dhms(walltime)
      }
    },

    #' @description recreate fwhmx objects and completion status
    #' @details useful for updating job status and ACF params after a run completes
    refresh = function(...) {
      #N.B. This should refresh from the object generated by the batch... basically the batch runs, then saves the modified copy in an RData object. Then we
      # pull this in from that object and refresh the parent/calling object.
      # so if output_rdata_file is /tmp/test.RData, we need to pull cobj out of that, then get its results into this object.
      # basically, the parent object should always know where to look for the output RData of the batch job

      # create 3dFWHMx object for each input file
      private$fwhmx_objs <- lapply(seq_along(private$input_files), function(ii) {
        fwhmx$new(input_file = private$input_files[ii], mask_file = private$mask_files[ii], ...)
      })

      private$all_fwhmx_complete <- all(sapply(private$fwhmx_objs, function(x) x$is_fwhmx_complete()))
    },

    #' @description Run the 3dFWHMx batch for all inputs
    #' @param force If \code{TRUE}, completed 3dFWHMx runs will be included in the batch. Default: FALSE
    run = function(force = FALSE) {
      private$generate_batch(include_complete = force)
      if (!is.null(private$fwhmx_batch)) {
        private$fwhmx_batch$submit()
      } else {
        message("Nothing to run")
      }
    },

    #' @description Return R_batch_job objects used to run 3dFWHMx for all inputs
    #' @param include_complete If TRUE, complete 3dFWHMx jobs will be included in the batch, leading these to be re-run
    get_batch = function(include_complete = FALSE) {
      private$generate_batch(include_complete = include_complete)
      return(private$fwhmx_batch)
    },
    get_acf_average = function(allow_incomplete_fwhmx=FALSE) {
      private$populate_acf_df()
      if (any(is.na(private$fwhmx_df)) && isFALSE(allow_incomplete_fwhmx)) {
        stop("Unable to get ACF average because parameters are missing for some datasets. Make sure you've run 3dFWHMx for all inputs provided!")
      }

      cmeans <- colMeans(private$fwhmx_df[, c("a", "b", "c", "effective_fwhm")], na.rm = TRUE)
      private$fwhmx_acf_avg <- cmeans[c("a", "b", "c")]
      private$fwhmx_effective_fwhm <- cmeans["effective_fwhm"]
      return(private$fwhmx_acf_avg)
    },
    get_acf_df = function() {
      private$populate_acf_df()
      private$fwhmx_df
    },
    get_effective_fwhm = function() {
      if (is.na(private$fwhmx_effective_fwhm)) {
        private$populate_acf_df()
      }
      return(private$fwhmx_effective_fwhm)
    },

    #' @description Simple method to return whether 3dFWHMx is complete for all input files
    is_fwhmx_complete = function() {
      private$all_fwhmx_complete
    }
  )
)


#' R6 class for 3dClustSim automation
#' @keywords internal
#' @importFrom tibble tibble
clustsim_spec <- R6::R6Class("clustsim_spec",
  private = list(
    out_dir = getwd(),
    clustsim_df = tibble::tibble(),
    out_files = NULL,
    prefix = "clustsim_",
    use_fwhmx_acf = FALSE,
    acf_params = setNames(rep(NA_real_, 3), c("a", "b", "c")),
    clustsim_batch = NULL,
    clustsim_call = NULL,
    clustsim_mask = NULL,
    clustsim_complete = FALSE,
    method = "fwhmx", # or 'inset' for 3dttest++ approach or 'xyz' for voxel + matrix approach
    nopad = "",
    pthr = ".02 .01 .005 .002 .001 .0005 .0002 .0001",
    athr = ".05 .02 .01 .005 .002 .001 .0005 .0002 .0001",
    iter = 10000L,
    nodec = "",
    seed = 0,
    sim_string = "",
    scheduler = "slurm",
    ncpus = 1L,
    which_sim = function() {
      if (!is.null(self$inset_files)) {
        private$sim_string <- glue("-inset {paste(self$inset_files, collapse=' ')}")
        if (!is.null(private$clustsim_mask)) {
          private$sim_string <- paste(private$sim_string, glue::glue("-mask {private$clustsim_mask}")) # -mask compatible with -inset
        }
      } else if (!is.null(private$clustsim_mask)) {
        private$sim_string <- glue::glue("-mask {private$clustsim_mask}")
      } else {
        # use voxel size and matrix size
        private$sim_string <- glue::glue("-nxyz {paste(private$nxyz, collapse=' ')} -dxyz {paste(private$dxyz, collapse=' ')}")
      }
    },
    build_call = function() {
      if (isTRUE(private$use_fwhmx_acf)) {
        if (isFALSE(self$fwhmx_set$is_fwhmx_complete())) {
          warning("Cannot build 3dClustSim call because some 3dFWHMx runs are incomplete.")
          return(invisible(NULL))
        }
        acf_string <- glue::glue("-acf {paste(self$fwhmx_set$get_acf_average(), collapse=' ')}")
      } else {
        acf_string <- glue::glue("-acf {paste(private$acf_params, collapse=' ')}")
      }

      private$which_sim()
      private$clustsim_call <- glue::glue(
        "3dClustSim {private$sim_string} -iter {private$iter} {acf_string} -prefix {private$prefix}",
        " -pthr {private$pthr} -athr {private$athr} -seed {private$seed}{private$nopad}{private$nodec}"
      )
    },
    build_df = function() {
      if (isTRUE(private$clustsim_complete) && nrow(private$clustsim_df) > 0L) {
        return(invisible(NULL)) # data frame has already been compiled
      }

      clustsim_files <- list.files(path = private$out_dir, pattern = glue::glue("{private$prefix}.*sided\\.1D"), full.names = TRUE)
      if (length(clustsim_files) == 0L) {
        warning("Cannot find any clustsim output files in: ", private$out_dir)
        return(invisible(NULL))
      } else if (length(clustsim_files) == 9L) {
        private$clustsim_complete <- TRUE
      } else {
        warning("Fewer than 9 clustsim files output files detected: ", paste(clustsim_files, collapse=", "))
      }

      read_clustsim_1d <- function(file) {
        checkmate::assert_file_exists(file)
        nn <- as.integer(sub("^.*NN([1-3])_.*", "\\1", file, perl = TRUE))
        sided <- sub("^.*NN[1-3]_(1|2|bi)sided.*", "\\1", file, perl = TRUE)
        data <- read.table(file, header = FALSE)
        header <- grep(pattern = "^\\s*#.*", readLines(file), value = TRUE)
        p_line <- gsub("[#|]", "", grep(pattern = "^\\s*#\\s*pthr.*", header, value = TRUE)) # strip # and | from header line
        c_names <- strsplit(trimws(p_line), "\\s+")[[1]]
        stopifnot(length(c_names) == ncol(data))
        names(data) <- c_names
        data_long <- data %>%
          tidyr::pivot_longer(cols = !pthr, names_to = "athr", values_to = "nvoxels") %>%
          dplyr::mutate(athr = as.numeric(athr), nn = !!nn, sided = !!sided) %>%
          dplyr::select(nn, sided, pthr, athr, nvoxels)
        return(data_long)
      }

      # read each clustsim output and compile into single indexed data.frame
      private$clustsim_df <- dplyr::bind_rows(lapply(clustsim_files, read_clustsim_1d))
    }
  ),
  public = list(
    #' @field fwhmx_set A fwhmx_set_spec object containing 3dFWHMx information for all fwhmx_input_files
    fwhmx_set = NULL,

    #' @field inset_files A character vector of files to use directly as volumes to threshold and clusterize
    inset_files = NULL,

    initialize = function(out_dir = NULL, prefix=NULL,
                          fwhmx_input_files = NULL, fwhmx_mask_files = NULL, inset_files = NULL, dxyz = NULL, nxyz = NULL,
                          clustsim_mask = NULL, acf_params = NULL, nopad = NULL, pthr = NULL, athr = NULL, iter = NULL, nodec = NULL, 
                          seed = NULL, scheduler = NULL, ncpus = NULL) {

      # only one method can be used for the volume over which 3dClustSim simulates. (the first two are parts of the same method)
      n_passed <- sapply(list(dxyz, nxyz, clustsim_mask, inset_files), function(x) !is.null(x))
      methods_passed <- sum(n_passed[1] || n_passed[2], n_passed[3], n_passed[4])

      if (methods_passed != 1L) {
        stop("Only one method can be passed for the 3dClustSim volume settings. Use dxyz + nxyz, clustsim_mask, or inset_files")
      }

      if (!is.null(scheduler)) {
        checkmate::assert_string(scheduler)
        checkmate::assert_subset(scheduler, c("torque", "slurm", "local"))
        private$scheduler <- scheduler
      }

      if (!is.null(out_dir)) {
        checkmate::assert_directory_exists(out_dir)
        private$out_dir <- out_dir
      }

      if (!is.null(prefix)) {
        checkmate::assert_string(prefix)
        private$prefix <- prefix
      }

      if (!is.null(fwhmx_input_files)) {
        private$use_fwhmx_acf <- TRUE
        self$fwhmx_set <- fwhmx_set_spec$new(input_files = fwhmx_input_files, mask_files = fwhmx_mask_files, scheduler=scheduler)
        if (!is.null(acf_params)) {
          stop("Cannot pass fwhmx_input_files and acf_params since the ACF parameters are calculated by 3dFWHMx!")
        }
      }

      if (!is.null(inset_files)) {
        if (!is.null(fwhmx_input_files)) {
          stop("Cannot specify both inset_files and fwhmx_input_files. See 3dClustSim documentation!")
        }

        private$method <- "inset"
        checkmate::assert_file_exists(inset_files)
        self$inset_files <- inset_files
      }

      if (!is.null(dxyz) || !is.null(nxyz)) {
        private$method <- "xyz"
        if (!is.null(dxyz)) {
          private$dxyz <- check_nums(dxyz, lower = .01, upper = 100) # millimeters
        } else {
          stop("Spatial domain 3dClustSim specification requires both dxyz and nxyz. I only received nxyz.")
        }

        if (!is.null(nxyz)) {
          private$nxyz <- check_nums(nxyz, lower = 1, upper = 10000) # voxels
        } else {
          stop("Spatial domain 3dClustSim specification requires both dxyz and nxyz. I only received dxyz.")
        }
      }

      if (!is.null(clustsim_mask)) {
        checkmate::assert_string(clustsim_mask)
        checkmate::assert_file_exists(clustsim_mask)
        private$clustsim_mask <- clustsim_mask
      } else {
        message("No clustsim_mask argument received. All voxels in volume will be used for cluster simulations!")
      }

      if (!is.null(acf_params)) {
        checkmate::assert_numeric(acf_params, len=3)
      }

      if (!is.null(nopad)) {
        checkmate::assert_logical(nopad, len = 1L)
        if (isTRUE(nopad)) private$nopad <- " -nopad"
      }

      check_nums <- function(inp, lower = 0, upper = 1e10) {
        inp_name <- deparse(substitute(inp)) # get name of object passed in
        if (checkmate::test_string(inp)) {
          inp <- suppressWarnings(as.numeric(strsplit(inp, "\\s+")[[1]]))
        } else if (checkmate::test_character(inp)) {
          inp <- suppressWarnings(as.numeric(inp))
        }

        if (checkmate::test_numeric(inp, lower = lower, upper = upper, any.missing = FALSE)) {
          return(paste(inp, collapse = " "))
        } else {
          stop("Problem with ", inp_name, " specification: ", paste(inp, collapse = " "))
        }
      }

      if (!is.null(pthr)) {
        private$pthr <- check_nums(pthr, lower = 1e-10, upper = .999)
      }

      if (!is.null(athr)) {
        private$athr <- check_nums(athr, lower = 1e-10, upper = .999)
      }

      if (!is.null(iter)) {
        checkmate::assert_integerish(iter, lower = 10, upper = 1e8)
        private$iter <- as.integer(iter)
      }

      if (!is.null(nodec)) {
        checkmate::assert_logical(nodec, len = 1L)
        if (isTRUE(nodec)) private$nodec <- " -nodec"
      }

      if (!is.null(seed)) {
        checkmate::assert_integerish(seed, len = 1L, lower=0)
        private$seed <- as.integer(seed)
      }

      if (!is.null(ncpus)) {
        checkmate::assert_integerish(ncpus, lower = 1, len = 1)
        private$ncpus <- ncpus
      }
    },
    run = function(force=FALSE) {
      if (isTRUE(private$clustsim_complete) && isFALSE(force)) {
        message("3dClustSim already finished for these inputs. Use $run(force=TRUE) to re-run or $get_clustsim_df() to retrieve results.")
        return(invisible(NULL))
      }

      # create batch job
      private$clustsim_batch <- R_batch_job$new(
        job_name = "run_3dclustsim", n_cpus = private$ncpus,
        cpu_time = private$walltime, scheduler = private$scheduler,
        input_objects = list(csim_obj = self), # pass the current object as input to the batch
        r_packages = "fmri.pipeline",
        r_code = c(
          "# refresh completeness of 3dFWHMx runs in case these were run by a preceding batch job",
          "if (isTRUE(csim_obj$get_use_fwhmx_acf())) csim_obj$fwhmx_set$refresh()",
          "setwd(csim_obj$get_out_dir())",
          "run_afni_command(csim_obj$get_call(), omp_num_threads = csim_obj$get_ncpus())"
        )
      )

      # need to run 3dFWHMx before 3dClustSim can run
      if (isTRUE(private$use_fwhmx_acf) && isFALSE(self$fwhmx_set$is_fwhmx_complete())) {
        fwhmx_batch <- self$fwhmx_set$get_batch()
        batch_obj <- private$clustsim_batch
        batch_obj$depends_on_parents <- "run_3dfwhmx"
        clust_seq <- R_batch_sequence$new(fwhmx_batch, batch_obj)
        clust_seq$submit()
        # clust_seq$generate()
      } else {
        private$clustsim_batch$submit()
      }
      
      return(invisible(self))
    },
    get_clustsim_df = function() {
      private$build_df()
      private$clustsim_df
    },
    get_call = function() {
      private$build_call()
      private$clustsim_call
    },
    get_out_dir = function() {
      private$out_dir
    },
    get_ncpus = function() {
      private$ncpus
    },
    get_use_fwhmx_acf = function() {
      private$use_fwhmx_acf
    }

  )
)
