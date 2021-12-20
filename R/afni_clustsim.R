# N.B. This class doesn't even expose the Gaussian ACF options given false positive problems
fwhmx_spec <- R6::R6Class("fwhmx_spec",
  private = list(
    data_file = NULL,
    out_dir = NULL, # where to put calculations
    out_summary = NULL,
    out_detailed = NULL,
    out_by_subbrik = NULL,
    mask_file = NULL,
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
        mask_string <- glue("-mask {private$mask_file}")
      }
      private$call <- glue(
        "3dFWHMx -overwrite -acf {private$out_detailed}",
        " -out {private$out_by_subbrik}",
        " -input {private$data_file} {mask_string} {private$average}",
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
    }
  ),
  public = list(
    initialize = function(data_file = NULL, out_dir = NULL, mask_file = NULL, demed = NULL, unif = NULL, average = "geometric") {
      if (is.null(data_file)) {
        stop("Cannot run 3dFWHMx without -input dataset.")
      } else {
        checkmate::assert_file_exists(data_file)
        private$data_file <- normalizePath(data_file)
      }

      if (is.null(out_dir)) {
        private$out_dir <- dirname(private$data_file)
      } else {
        private$out_dir <- normalizePath(out_dir)
      }

      private$out_summary <- file.path(private$out_dir, "3dFWHMx_acf.txt")
      private$out_by_subbrik <- file.path(private$out_dir, "3dFWHMx_acf_bysubbrik.txt")
      private$out_detailed <- file.path(private$out_dir, "3dFWHMx_acf_radius.txt")

      if (is.null(mask_file)) {
        message("No mask_file provided. Will default to 3dFWHMx -automask")
      } else {
        checkmate::assert_file_exists(mask_file)
        private$mask_file <- mask_file
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
      private$data_file
    },
    get_outputs = function() {
      c(
        out_dir = private$out_dir,
        out_summary = private$out_summary,
        out_detailed = private$out_detailed,
        out_by_subbrik = private$out_by_subbrik
      )
    }
  )
)



clustsim_spec <- R6::R6Class("clustsim_spec",
  private = list(
    fwhmx_objs = NULL,
    fwhmx_df = data.frame(),
    fwhmx_acf_avg = setNames(rep(NA_real_, 4), c("a", "b", "c")),
    fwhmx_effective_fwhm = NA_real_,
    clustsim_call = NULL,
    method = "fwhmx", # or 'inset' for 3dttest++ approach or 'xyz' for voxel + matrix approach
    nopad = "",
    pthr = ".02 .01 .005 .002 .001 .0005 .0002 .0001",
    athr = ".05 .02 .01 .005 .002 .001 .0005 .0002 .0001",
    clustsim_mask = NULL,
    iter = 10000L,
    sim_string = "",
    scheduler = "slurm",
    which_sim = function() {
      if (!is.null(self$inset_files)) {
        private$sim_string <- glue("-inset {paste(self$inset_files, collapse=' ')}")
        if (!is.null(private$clustsim_mask)) {
          private$sim_string <- paste(private$sim_string, glue("-mask {private$clustsim_mask}")) # -mask compatible with -inset
        }
      } else if (!is.null(private$clustsim_mask)) {
        private$sim_string <- glue("-mask {private$clustsim_mask}")
      } else {
        # use voxel size and matrix size
        private$sim_string <- glue("-nxyz {paste(private$nxyz, collapse=' ')} -dxyz {paste(private$dxyz, collapse=' ')}")
      }
    },
    build_call = function() {
      if (private$method == "fwhmx") {

      }
      private$which_sim()
      private$get_acf_average()
      private$clustsim_call <- glue(
        "3dClustSim {private$sim_string} -iter {private$iter} -acf {private$fwhmx_acf_avg}",
        " -pthr {private$pthr} -athr {private$athr} {private$nopad}"
      )
      browser()
    },
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
    }
  ),
  public = list(
    #' @field fwhmx_files A character vector of files from run-level data to be passed to 3dFWHMx (these should usually be residuals)
    fwhmx_files = NULL,

    #' @field A character vector of masks for each element of \code{fwhmx_files} -- passed as -mask to 3dFWHMx
    mask_files = NULL,

    #' @field inset_files A character vector of files to use directly as volumes to threshold and clusterize
    inset_files = NULL,
    initialize = function(fwhmx_files = NULL, mask_files = NULL, inset_files = NULL, dxyz = NULL, nxyz = NULL,
                          clustsim_mask = NULL, nopad = NULL, pthr = NULL, athr = NULL, iter = NULL, scheduler = NULL) {
      if (!is.null(fwhmx_files)) {
        checkmate::assert_character(fwhmx_files)
        checkmate::assert_file_exists(fwhmx_files)
        private$method <- "fwhmx"
        self$fwhmx_files <- fwhmx_files
      }

      if (!is.null(mask_files)) {
        checkmate::assert_file_exists(mask_files)
        stopifnot(length(mask_files) == length(fwhmx_files)) # enforce length match for inputs and masks (otherwise, how do they line up?)
        self$mask_files <- mask_files
      }

      if (!is.null(inset_files)) {
        if (!is.null(fwhmx_files)) {
          stop("Cannot specify both inset_files and fwhmx_files. See 3dClustSim documentation!")
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

      if (!is.null(nopad)) {
        checkmate::assert_logical(nopad, len = 1L)
        if (isTRUE(nopad)) private$nopad <- "-nopad"
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

      if (!is.null(pthr)) {
        private$athr <- check_nums(athr, lower = 1e-10, upper = .999)
      }

      if (!is.null(iter)) {
        checkmate::assert_integerish(iter, lower = 10, upper = 1e8)
        private$iter <- as.integer(iter)
      }

      # if (is.null(private$fwhmx_objs) || !inherits(private$fwhmx_objs, "list") || length(private$fwhx_objs) != length(self$fwhmx_files)) {
      if (private$method == "fwhmx") {
        # execute 3dFWHMx for each input file
        private$fwhmx_objs <- lapply(seq_along(self$fwhmx_files), function(ii) {
          fwhmx_spec$new(data_file = self$fwhmx_files[ii], mask_file = self$mask_files[ii])
        })
      }

      if (!is.null(scheduler)) {
        checkmate::assert_string(scheduler)
        checkmate::assert_subset(scheduler, c("torque", "slurm", "local"))
        private$scheduler <- scheduler
      }
    },
    run_3dfwhmx = function(force = FALSE) {
      if (private$method != "fwhmx") {
        stop("This clustsim object does not appear to be setup for 3dFHWMx inputs!")
      }

      # this will return all 3dfwhmx calls
      fwhmx_calls <- sapply(private$fwhmx_objs, function(x) x$get_call())

      # create parent batch job to run all 3dFWHMx scripts
      fwhmx_batch <- R_batch_job$new(
        job_name = "run_3dfwhmx", n_cpus = 1,
        cpu_time = "10:00:00", scheduler = private$scheduler,
        input_objects = "fwhmx_calls", # export this object to the job
        wait_for_children = TRUE, r_packages="fmri.pipeline",
        r_code = glue(
          "child_job_ids <- cluster_submit_shell_jobs(fwhmx_calls, memgb_per_command=8, fork_jobs=TRUE, scheduler='{private$scheduler}')"
        )
      )

      fwhmx_batch$submit()
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
    run = function() {

    },
    get_call = function() {
      browser()
      private$build_call()
      private$clustsim_call
    }
  )
)

setwd("/proj/mnhallqlab/studies/MMClock/MR_Proc/10637_20140304/mni_5mm_aroma/sceptic_vchosen_ventropy_dauc_pemax_vtime_preconvolve")
res4d_files <- list.files(pattern = "res4d.nii.gz", getwd(), full.names = T, recursive = T)
mask_files <- list.files(pattern = "mask.nii.gz", getwd(), full.names = T, recursive = T)

mytest <- clustsim_spec$new(fwhmx_files = res4d_files, mask_files = mask_files)
