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

      demed_string <- ifelse(isTRUE(private$demed), "-demed", "")
      unif_string <- ifelse(isTRUE(private$unif), "-unif", "")

      private$call <- glue::glue(
        "3dFWHMx -overwrite -acf {private$out_detailed}",
        " -out {private$out_by_subbrik}",
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

      # always return self for side-effect methods
      return(invisible(self))
    }
  ),
  public = list(
    #' @param input_file The input dataset whose smoothness should be calculated (often first-level GLM residuals)
    #' @param mask_file Only compute smoothness within this mask (if not provided, -automask will be used)
    #' @param out_dir The output directory for fwhmx files
    #' @param demed
    #' 
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
    get_mask_file = function() {
      private$mask_file
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