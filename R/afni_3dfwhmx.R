#' R6 class for running 3dFWHMx on a single input file based on user specification
#'
#' @importFrom R6 R6Class
#' @details N.B. This class doesn't even expose the Gaussian ACF options given false positive problems
#' @export
afni_3dfwhmx <- R6::R6Class("afni_3dfwhmx",
  private = list(
    input_file = NULL,
    mask_file = NULL,
    out_dir = NULL, # where to put calculations
    out_summary = NULL,
    out_by_radius = NULL,
    out_by_volume = NULL,
    demed = FALSE,
    unif = FALSE,
    average = "-geom",
    call = NULL,
    acf_params = rep(NA_real_, 4), # should always be 4 in length
    fwhm_by_volume = NULL,
    acf_by_radius = NULL,
    ret_code = NULL,
    fwhmx_complete = FALSE,
    ncpus = 1L,
    build_call = function() {
      if (is.null(private$mask_file)) {
        mask_string <- "-automask"
      } else {
        mask_string <- glue::glue("-mask {private$mask_file}")
      }

      demed_string <- ifelse(isTRUE(private$demed), "-demed", "")
      unif_string <- ifelse(isTRUE(private$unif), "-unif", "")

      private$call <- glue::glue(
        "3dFWHMx -overwrite -acf {private$out_by_radius}",
        " -out {private$out_by_volume}",
        " -input {private$input_file} {mask_string} {private$average} {demed_string} {unif_string}",
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
    },

    # read the detailed output (expanding radius) from the file
    populate_by_radius = function(force=FALSE) {
      if (isTRUE(private$fwhmx_complete) && isFALSE(force) && is.data.frame(private$acf_by_radius)) {
        # don't re-read detailed object unless requested explicitly
        return(invisible(NULL))
      }

      if (checkmate::test_file_exists(private$out_by_radius)) {
        private$acf_by_radius <- read.table(private$out_by_radius) %>%
          setNames(c("radius_mm", "acf_r", "model_r", "gaussian_newmodel_r"))
      } else {
        warning("Cannot locate ACF detailed params in file: ", private$out_by_radius)
      }
    },

    # read the sub-brik ACF outputs from the file
    populate_by_volume = function(force = FALSE) {
      if (isTRUE(private$fwhmx_complete) && isFALSE(force) && is.data.frame(private$fwhm_by_volume)) {
        # don't re-read detailed object unless requested explicitly
        return(invisible(NULL))
      }

      if (checkmate::test_file_exists(private$out_by_volume)) {
        private$fwhm_by_volume <- read.table(private$out_by_volume) %>%
          setNames(c("fwhm_x", "fwhm_y", "fwhm_z"))
      } else {
        warning("Cannot locate ACF sub-brik params in file: ", private$out_by_volume)
      }
    }
  ),
  public = list(
    #' @param input_file The input dataset whose smoothness should be calculated (often first-level GLM residuals)
    #' @param mask_file Only compute smoothness within this mask (if not provided, -automask will be used)
    #' @param out_dir The output directory for fwhmx files
    #' @param demed If \code{TRUE}, subtract the median of each voxels time series before calculating 
    #'   FWHM (wraps -demed). Default: FALSE.
    #' @param unif If \code{TRUE}, normalize each voxel's time series to have the same MAD before 
    #'   calculating FWHM (wraps -unif). Default: FALSE
    #' @param average For multi-volume data, compute the averaged ACF estimates by either the geometric
    #'   or arithmetic mean. Default: "geometric".
    #' @param ncpus The number of threads/cores to use for running 3dFWHMx. Default is 1. Controls OMP_NUM_THREADS.
    initialize = function(input_file = NULL, mask_file = NULL, out_dir = NULL, 
    demed = NULL, unif = NULL, average = "geometric", ncpus = 1L) {
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
        private$mask_file <- normalizePath(mask_file)
      }

      if (is.null(out_dir)) {
        private$out_dir <- dirname(private$input_file)
      } else {
        private$out_dir <- normalizePath(out_dir)
      }

      private$out_summary <- file.path(private$out_dir, "3dFWHMx_acf.txt")
      private$out_by_volume <- file.path(private$out_dir, "3dFWHMx_acf_bysubbrik.txt")
      private$out_by_radius <- file.path(private$out_dir, "3dFWHMx_acf_radius.txt")

      if (file.exists(private$out_summary)) {
        # populate from cached txt file (rather than forcing the command to re-run)
        private$populate_params()
      }

      # at present, use JIT population of by-radius and by-subbrik outputs if user calls their $get methods
      # if (file.exists(private$out_by_radius)) {
      #   # populate from cached txt file (rather than forcing the command to re-run)
      #   private$populate_by_radius()
      # }

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

      if (!is.null(ncpus)) {
        checkmate::assert_integerish(ncpus, lower = 1, upper = 1e3, len = 1)
        private$ncpus <- ncpus
      }
    },

    #' @description runs 3dFWHMx on this input dataset
    #' @param force If \code{TRUE}, 3dFWHMx will be run even if the expected output
    #'   files already exist.
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
      private$ret_code <- run_afni_command(private$call, omp_num_threads = private$ncpus)
      if (private$ret_code != 0) {
        warning("3dFWHMx call returned error code.")
      } else {
        # Assume that we have gotten good returns. We should mark the object as complete and populate ACF params
        private$populate_params()
      }

      return(invisible(self))
    },

    #' @description return the 3dFWHMx call used for the specified input
    #' @details this is useful if you want to call 3dFWHMx yourself directly or if you want to debug
    #'   the 3dFWHMx call specification.
    #' @return a character string with the 3dFWHMx call
    get_call = function() {
      private$build_call()
      return(private$call)
    },

    #' @description return the estimated ACF parameters for this run of data, averaged over volumes
    #' @details Will issue a warning if 3dFWHMx has not run successfully on this dataset already
    #' @return a three-element vector containing the ACF estimates for the dataset, averaging
    #'   over volumes.
    get_acf_params = function() {
      private$populate_params()
      if (isFALSE(private$fwhmx_complete)) {
        warning("3dFWHMx has not run successfully for this dataset. Please use the $run() method first!")
      }
      private$acf_params
    },

    #' @description return the estimated ACF parameters for this run of data
    #' @details Will issue a warning if 3dFWHMx has not run successfully on this dataset already
    #' @return a three-element vector containing the ACF estimates for the dataset, averaging
    #'   over volumes.
    get_acf_by_radius = function() {
      private$populate_by_radius()
      if (isFALSE(private$fwhmx_complete)) {
        warning("3dFWHMx has not run successfully for this dataset. Please use the $run() method first!")
      }
      private$acf_by_radius
    },

    #' @description return the estimated ACF parameters for this run of data
    #' @details Will issue a warning if 3dFWHMx has not run successfully on this dataset already
    #' @return a three-element vector containing the ACF estimates for the dataset, averaging
    #'   over volumes.
    get_fwhm_by_volume = function() {
      private$populate_by_volume()
      if (isFALSE(private$fwhmx_complete)) {
        warning("3dFWHMx has not run successfully for this dataset. Please use the $run() method first!")
      }
      private$fwhm_by_volume
    },

    #' @description return the input file (NIfTI) used for 3dFWHMx
    #' @return a character string of the input file location
    get_input_file = function() {
      private$input_file
    },

    #' @description return the mask file used for 3dFWHMx estimation
    #' @return a character string of the mask file location
    get_mask_file = function() {
      private$mask_file
    },

    #' @description return the expected output files related to this 3dFWHMx object
    #' @return a character vector containing expected output files
    get_outputs = function() {
      c(
        out_dir = private$out_dir,
        out_summary = private$out_summary,
        out_by_radius = private$out_by_radius,
        out_by_radius_png = paste0(private$out_by_radius, ".png"), # output by program automatically
        out_by_volume = private$out_by_volume
      )
    },

    #' method to indicate whether 3dFWHMx has already run and completed for this input
    is_complete = function() {
      private$fwhmx_complete
    },

    remove_outputs = function(prompt = FALSE) {
      checkmate::assert_logical(prompt, len = 1L)
      f_list <- self$get_outputs()
      f_exists <- sapply(f_list, checkmate::test_file_exists)
      f_list <- f_list[f_exists]
      if (length(f_list) == 0L) {
        message("No output files found for removal.")
      } else {
        if (isTRUE(prompt)) {
          for (ff in f_list) {
            cat("File:", ff, "\n")
            remove <- askYesNo("Delete?")
            if (isTRUE(remove)) {
              cat(sprintf("Removing: %s", ff), sep = "\n")
              unlink(ff)
            } else if (is.na(remove)) {
              # cancel
              return(invisible(NULL))
            }
          }
        } else {
          cat(sprintf("Removing: %s", f_list), sep="\n")
          unlink(f_list)
        }

      }
      private$fwhmx_complete <- FALSE # will need to re-run
      

    }
  )
)

# TESTS

# path <- "/proj/mnhallqlab/studies/MMClock/MR_Proc/10637_20140304/mni_5mm_aroma/sceptic_vchosen_ventropy_dauc_pemax_vtime_preconvolve/FEAT_LVL1_run1.feat"
# res4d_file <- file.path(path, "stats", "res4d.nii.gz")
# fwhmx_mask_file <- file.path(path, "mask.nii.gz")

# test <- afni_3dfwhmx$new(input_file = res4d_file, mask_file = fwhmx_mask_file, average = "geometric", ncpus = 2)
# test$run()
# test$is_complete()
# test$get_acf_params() # average ACF params
# test$get_fwhm_by_volume() # FWHM for each sub-brik (volume)
# test$get_acf_by_radius() # FWHM for each sub-brik (volume)
# test$remove_outputs(prompt = FALSE) # delete files generated by this input
