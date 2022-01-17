#' wrapper class for 3dClusterize
#' @export
clusterize_spec <- R6::R6Class("clusterize_spec",
  private = list(

    # private fields
    pvt_inset_file = NULL,
    pvt_mask = NULL,
    pvt_sided = "bi", # default to bi-sided (i.e., adjacent positive and negative clusters cannot merge)
    pvt_lower_thresh = NULL, # lower threshold value for clusterizing (two and bi)
    pvt_upper_thresh = NULL, # upper threshold value for clusterizing (two and bi)
    pvt_one_thresh = NULL, # single threshold value for clusterizing (one)
    pvt_one_tail = NULL, # which tail to test in case of one-sided (LEFT or RIGHT)
    pvt_mask_from_hdr = NULL,
    pvt_out_mask = NULL,
    pvt_threshold_file = NULL, # two-part alternative to inset file -- stitched together by this class before 3dClusterize
    pvt_data_file = NULL, # two-part alternative to inset file -- stitched together by this class before 3dClusterize
    pvt_ithr = NULL,
    pvt_idat = NULL,
    pvt_NN = 1L, # faces + edges + corners must touch
    pvt_clust_nvox = 10, # very loose lower bound by default (tiny clusters)
    pvt_pref_map = NULL,
    pvt_pref_dat = NULL,
    pvt_1Dformat = TRUE,
    pvt_quiet = NULL,
    pvt_orient = "LPI", #AFNI default is RAI, which seems like a silly default to me
    pvt_binary = FALSE,
    pvt_clusterize_call = NULL,
    pvt_clusterize_output_file = "clusters.1D",
    pvt_dirty = FALSE, # whether object fields have changed since last build_call

    # private methods
    set_sided = function(val) {
      checkmate::assert_string(val)
      checkmate::assert_subset(val, c("1", "one", "2", "two", "bi"))
      if (val == "1") {
        val <- "one" # convert to standard nomenclature
      } else if (val == "2") {
        val <- "two"
      }
      private$pvt_sided <- val
      private$pvt_dirty <- TRUE # will rebuild call
    },
    set_NN = function(val) {
      if (checkmate::test_string(val)) {
        val <- as.integer(val)
      }

      checkmate::assert_integerish(val, lower = 1, upper = 3, len = 1L, any.missing = FALSE)
      private$pvt_NN <- as.integer(val)
      private$pvt_dirty <- TRUE # will rebuild call
    },
    set_clust_nvox = function(val) {
      if (checkmate::test_string(val)) {
        val <- as.integer(val)
      }

      checkmate::assert_integerish(val, lower = 1, upper = 1e5, len = 1L, any.missing = FALSE)
      private$pvt_clust_nvox <- as.integer(val)
      if (!is.null(private$pvt_clust_nvol)) {
        message("Resetting clust_nvol to NULL because clust_nvox was set.")
        private$pvt_clust_nvol <- NULL
      }
      private$pvt_dirty <- TRUE # will rebuild call
    },
    set_clust_nvol = function(val) {
      if (checkmate::test_string(val)) {
        val <- as.numeric(val)
      }

      checkmate::assert_numeric(val, lower = 0.5, upper = 1e6, len = 1L, any.missing = FALSE)
      private$pvt_clust_nvol <- as.numeric(val)
      if (!is.null(private$pvt_clust_nvox)) {
        message("Resetting clust_nvox to NULL because clust_nvol was set.")
        private$pvt_clust_nvox <- NULL
      }
      private$pvt_dirty <- TRUE # will rebuild call
    },
    build_call = function() {
      if (isFALSE(private$pvt_dirty)) return(invisible(NULL)) # no need to rebuild call
      combine_call <- NULL
      if (!is.null(private$pvt_threshold_file)) {
        if (is.null(private$pvt_data_file)) {
          in_file <- private$pvt_threshold_file # 3D input to clusterize
        } else {
          # build tmp dataset that combined threshold and data file
          in_file <- tempfile(pattern = "clustinput", fileext = ".nii.gz")
          combine_call <- glue("3dTcat -output {in_file} {private$pvt_threshold_file} {private$pvt_data_file}")
        }
      } else {
        in_file <- private$pvt_inset_file
      }

      # if output file is provided without any path specification, use the directory of the input file
      private$pvt_clusterize_output_file <- R.utils::getAbsolutePath(private$pvt_clusterize_output_file, workDirectory = dirname(in_file))

      str <- glue("3dClusterize -inset {in_file} -ithr {private$pvt_ithr} -NN {private$pvt_NN}")
      if (!is.null(private$pvt_mask)) str <- glue("{str} -mask {private$pvt_mask}")
      if (isTRUE(private$pvt_mask_from_hdr)) str <- glue("{str} -mask_from_hdr")
      if (!is.null(private$pvt_out_mask)) str <- glue("{str} -out_mask {private$pvt_out_mask}")
      
      # I do not support -within_range yet...

      if (!is.null(private$pvt_idat)) str <- glue("{str} -idat {private$pvt_idat}")
      if (private$pvt_sided == "one") {
        str <- glue("{str} -1sided {private$pvt_one_thresh} {private$pvt_one_tail}")
      } else if (private$pvt_sided == "two") {
        str <- glue("{str} -2sided {private$pvt_lower_thresh} {private$pvt_upper_thresh}")
      } else if (private$pvt_sided == "bi") {
        str <- glue("{str} -bisided {private$pvt_lower_thresh} {private$pvt_upper_thresh}")
      }

      if (!is.null(private$pvt_clust_nvox)) {
        str <- glue("{str} -clust_nvox {private$pvt_clust_nvox}")
      } else if (!is.null(private$pvt_clust_nvol)) {
        str <- glue("{str} -clust_nvol {private$pvt_clust_nvol}")
      } else {
        stop("Both clust_nvox and clust_nvol are missing... This shouldn't happen!")
      }

      if (!is.null(private$pvt_pref_map)) str <- glue("{str} -pref_map {private$pvt_pref_map}")
      if (!is.null(private$pvt_pref_dat)) str <- glue("{str} -pref_dat {private$pvt_pref_dat}")

      # I may need to just fix this to TRUE for parsing to go smoothly.
      if (isTRUE(private$pvt_1Dformat)) {
        str <- glue("{str} -1Dformat")
      } else {
        str <- glue("{str} -no_1Dformat")
      }

      if (isTRUE(private$pvt_quiet)) str <- glue("{str} -quiet")
      if (!is.null(private$pvt_orient)) str <- glue("{str} -orient {private$pvt_orient}")
      if (isTRUE(private$pvt_binary)) str <- glue("{str} -binary")

      str <- glue(str, " > {private$pvt_clusterize_output_file}")

      private$pvt_clusterize_call <- c(combine = combine_call, clusterize = str)

    }
  ),

  # active bindings (these are the only things that can change after initialization)
  active = list(
    #' @field fwe_p a vector of p-values used for familywise error (FWE) z-statistic threshold calculations in pTFCE.
    sided = function(val) {
      if (missing(val)) {
        return(private$pvt_sided)
      } else {
        private$set_sided(val)
      }
    },
    NN = function(val) {
      if (missing(val)) {
        return(private$pvt_NN)
      } else {
        private$set_NN(val)
      }
    },
    clust_nvox = function(val) {
      if (missing(val)) {
        return(private$pvt_clust_nvox)
      } else {
        private$set_clust_nvox(val)
      }
    },
    clust_nvol = function(val) {
      if (missing(val)) {
        return(private$pvt_clust_nvol)
      } else {
        private$set_clust_nvol(val)
      }
    }
  ),
  #' @param inset A 4D dataset containing the statistic to use for thresholding (ithr) and, optionally, the data value to output/retain
  #' @param mask If specified, the volume will be masked by \code{mask} prior to clusterizing
  #' @param threshold_file A 3D dataset containing the statistic to use for thresholding. Mutually exclusive with \code{inset}
  #'   If passed, \code{ithr} and \code{idat} are ignored because the \code{inset} file is generated internally.
  #' @param data_file A 3D dataset containing the data value to be retained in clusters post-thresholding. 
  #'   Must be passed with \code{threshold_file} and will be stitched together with it internally. Mutually exclusive with \code{inset}.
  #' 
  public = list(
    initialize = function(inset = NULL, mask = NULL, threshold_file = NULL, data_file = NULL, mask_from_hdr = NULL, out_mask = NULL, 
      ithr = NULL, idat = NULL, onesided = NULL, twosided = NULL, bisided = NULL, 
      lower_thresh = NULL, upper_thresh = NULL, one_thresh = NULL, one_tail = NULL,
      NN = NULL, clust_nvox = NULL,
      clust_nvol = NULL, pref_map = NULL, pref_dat = NULL, x1Dformat = NULL, no_1Dformat = NULL, 
      quiet = NULL, orient = NULL, binary = NULL, clusterize_output_file = NULL) {
      
      if (!is.null(inset)) {
        checkmate::assert_file_exists(inset)
        private$pvt_inset_file <- normalizePath(inset)
      } else if (!is.null(threshold_file)) {
        checkmate::assert_file_exists(threshold_file)
        private$pvt_threshold_file <- normalizePath(threshold_file)

        if (!is.null(data_file)) {
          checkmate::assert_file_exists(data_file)
          private$pvt_data_file <- normalizePath(data_file)
          ithr <- 0 # always put threshold image as first sub-brik in the internal file
          idat <- 1
        } else {
          ithr <- 0
          idat <- NULL
        }
      } else {
        stop("You must pass either inset (combined threshold + data file) or data_file and threshold_file.")
      }

      if (!is.null(mask)) {
        checkmate::assert_file_exists(mask)
        private$pvt_mask <- mask
      }

      if (!is.null(mask_from_hdr)) {
        checkmate::assert_logical(mask_from_hdr, len = 1L)
        private$pvt_mask_from_hdr <- mask_from_hdr
      }

      if (!is.null(out_mask)) {
        checkmate::assert_string(out_mask, len = 1L)
        private$pvt_out_mask <- out_mask
      }

      checkmate::assert_integerish(ithr, lower=0, upper=1e4)
      private$pvt_ithr <- ithr

      if (!is.null(idat)) {
        checkmate::assert_integerish(idat, lower = 0, upper = 1e4)
        private$pvt_idat <- idat
      }

      n_sided <- sum(sapply(list(onesided, twosided, bisided), function(x) !is.null(x)))
      if (n_sided > 1L) {
        stop("Only one 'sided' option can be specified: onesided, twosided, or bisided.")
      } else if (!is.null(onesided)) {
        checkmate::assert_logical(onesided, len = 1L)
        private$set_sided("one")
        if (is.null(one_thresh)) {
          stop("Must pass in a one_thresh value used for thresholding image prior to clusterizing")
        } else {
          checkmate::assert_number(one_thresh)
          private$pvt_one_thresh <- one_thresh
        }

        if (is.null(one_tail)) {
          stop("Must pass in a one_tail value -- LEFT or RIGHT.")
        } else {
          checkmate::assert_subset(one_tail, c("LEFT", "RIGHT", "LEFT_TAIL", "RIGHT_TAIL"), empty.ok = FALSE)
          private$pvt_one_tail <- one_tail
        }
      } else if (!is.null(twosided)) {
        checkmate::assert_logical(twosided, len = 1L)
        private$set_sided("two")
      } else if (!is.null(bisided)) {
        checkmate::assert_logical(bisided, len = 1L)
        private$set_sided("bi")
      }

      if (private$pvt_sided %in% c("two", "bi")) {
        checkmate::assert_number(lower_thresh)
        checkmate::assert_number(upper_thresh)
        private$pvt_lower_thresh <- lower_thresh
        private$pvt_upper_thresh <- upper_thresh
      }

      # set NN active binding
      if (!is.null(NN)) {
        private$set_NN(NN)
      }

      # set clust_nvox active binding
      if (is.null(clust_nvox) && is.null(clust_nvol)) {
        stop("Both clust_nvox and clust_nvol were not provided. Cannot move forward with clusterize!")
      } else if (!is.null(clust_nvox)) {
        private$set_clust_nvox(clust_nvox)
      } else if (!is.null(clust_nvol)) {
        if (!is.null(clust_nvox)) {
          stop("Cannot pass both clust_nvox and clust_nvol!")
        }
        private$set_clust_nvol(clust_nvol)
      }

      # set the integer-valued cluster mask output name
      if (!is.null(pref_map)) {
        checkmate::assert_string(pref_map)
        ext <- file_ext(pref_map)
        if (is.na(ext)) {
          pref_map <- paste0(pref_map, ".nii.gz") # force nifti output
        } else if (grepl("(brik|head).*", ext)) {
          message("Forcing -pref_map extension to .nii.gz for consistency")
          pref_map <- paste0(file_sans_ext(pref_map), ".nii.gz") # force nifti output
        }
        private$pvt_pref_map <- pref_map
      }

      # set the data output, masked by the clusterize result
      if (!is.null(pref_dat)) {
        if (is.null(private$pvt_idat)) {
          stop("Asking for pref_dat requires that you also specify ithr")
        }

        checkmate::assert_string(pref_dat)
        ext <- file_ext(pref_dat)
        if (is.na(ext)) {
          pref_dat <- paste0(pref_dat, ".nii.gz") # force nifti output
        } else if (grepl("(brik|head).*", ext)) {
          message("Forcing -pref_dat extension to .nii.gz for consistency")
          pref_dat <- paste0(file_sans_ext(pref_dat), ".nii.gz") # force nifti output
        }

        private$pvt_pref_dat <- pref_dat
      }

      if (!is.null(x1Dformat)) {
        if (!is.null(no_1Dformat)) {
          stop("Cannot specify both x1Dformat and no_1Dformat!")
        }
        checkmate::assert_logical(x1Dformat, len = 1L)
        private$pvt_1Dformat <- TRUE
      } else if (!is.null(no_1Dformat)) {
        checkmate::assert_logical(no_1Dformat, len = 1L)
        private$pvt_1Dformat <- !no_1Dformat
      }

      if (!is.null(quiet)) {
        checkmate::assert_logical(quiet, len = 1L)
        private$pvt_quiet <- quiet
      }

      if (!is.null(binary)) {
        checkmate::assert_logical(binary, len = 1L)
        private$pvt_binary <- binary
      }

      if (!is.null(orient)) {
        checkmate::assert_string(orient)
        # should probably check LPI, RAI, etc.
        private$pvt_orient <- orient
      }

      if (!is.null(clusterize_output_file)) {
        checkmate::assert_string(clusterize_output_file)
        private$pvt_clusterize_output_file <- clusterize_output_file
      }

    },
    run = function() {
      private$build_call()
      run_afni_command(private$pvt_clusterize_call)
    },
    get_clust_df = function() {
      if (!checkmate::test_file_exists(private$pvt_clusterize_output_file)) {
        warning(glue("The expectected 3dClusterize output does not exist: {private$pvt_clusterize_output_file}"))
        return(invisible(data.frame()))
      }

      clust_df <- read.table(private$pvt_clusterize_output_file)
      clust_lines <- readLines(private$pvt_clusterize_output_file)
      comment_lines <- grep("^\\s*#", clust_lines, perl = TRUE)

      # header 
      cols <- 
      browser()

    },
    get_clust_nifti = function() {

    },
    get_call = function() {
      private$build_call()
      return(private$pvt_clusterize_call)
    }
  )
)

# x <- clusterize_spec$new(
#   threshold_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-pe/L2m-l2_l2c-emotion.happy/L3m-age_sex/FEAT_l1c-EV_pe.gfeat/cope1.feat/stats/zstat6.nii.gz",
#   lower_thresh = -3, upper_thresh = 3, bisided = TRUE, NN = 1, clust_nvox = 35, pref_map = "zstat_clusterize.nii.gz"
# )

# x$get_call()
# x$get_clust_df()
