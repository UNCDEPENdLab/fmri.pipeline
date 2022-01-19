#' wrapper class for 3dClusterize
#' @export
afni_3dclusterize <- R6::R6Class("afni_3dclusterize",
  private = list(

    # private fields
    pvt_input_file = NULL, # internal record of the input file to be passed as -inset to 3dClusterize
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
    pvt_1Dformat = TRUE, # currently not used since parser depends on 1D format
    pvt_quiet = NULL,
    pvt_orient = "LPI", #AFNI default is RAI, which seems like a silly default to me
    pvt_binary = FALSE,
    pvt_combine_call = NULL, # If needed, the 3dTcat call for combining pvt_threshold_file and pvt_data_file
    pvt_clusterize_call = NULL,
    pvt_clusterize_output_file = "clusters.1D",
    pvt_dirty_call = FALSE, # whether object fields have changed since last build_call
    pvt_whereami = NULL,

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
      private$pvt_dirty_call <- TRUE # will rebuild call
    },
    set_NN = function(val) {
      if (checkmate::test_string(val)) {
        val <- as.integer(val)
      }

      checkmate::assert_integerish(val, lower = 1, upper = 3, len = 1L, any.missing = FALSE)
      private$pvt_NN <- as.integer(val)
      private$pvt_dirty_call <- TRUE # will rebuild call
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
      private$pvt_dirty_call <- TRUE # will rebuild call
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
      private$pvt_dirty_call <- TRUE # will rebuild call
    },
    get_combine_call = function() {

    },
    build_call = function() {
      if (isFALSE(private$pvt_dirty_call)) return(invisible(NULL)) # no need to rebuild call
      
      str <- glue("3dClusterize -overwrite -inset {private$pvt_input_file} -ithr {private$pvt_ithr} -NN {private$pvt_NN} -1Dformat")
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
      if (isTRUE(private$pvt_quiet)) str <- glue("{str} -quiet")
      if (!is.null(private$pvt_orient)) str <- glue("{str} -orient {private$pvt_orient}")
      if (isTRUE(private$pvt_binary)) str <- glue("{str} -binary")

      str <- glue(str, " > {private$pvt_clusterize_output_file}")

      private$pvt_clusterize_call <- str
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

  #' @description initialization function for a new afni_3dclusterize object. Arguments largely mirror the 3dClusterize parameters.
  #' @param inset A 4D dataset containing the statistic to use for thresholding (ithr) and, optionally, the data value to output/retain
  #' @param mask If specified, the volume will be masked by \code{mask} prior to clusterizing
  #' @param threshold_file A 3D dataset containing the statistic to use for thresholding. Mutually exclusive with \code{inset}
  #'   If passed, \code{ithr} and \code{idat} are ignored because the \code{inset} file is generated internally.
  #' @param data_file A 3D dataset containing the data value to be retained in clusters post-thresholding. 
  #'   Must be passed with \code{threshold_file} and will be stitched together with it internally. Mutually exclusive with \code{inset}.
  #' @param mask_from_hdr passes through as -mask_from_hdr
  #' @param out_mask passes through as -out_mask
  #' @param ithr sub-brik number for the voxelwise threshold. Passes through as -ithr
  #' @param onesided if TRUE, clusterizing will be conducted on one tail of the statistic distribution (-ithr)
  #' @param twosided if TRUE, clusterizing will be conducted on both tails of the statistic distribution (-ithr)
  #' @param pref_map File name for the integer-valued mask containing each cluster, ordered by descending voxel size. 
  #'   Passes through as -pref_map.
  #' @param NN 1, 2, 3. Default: 1. Passes through as -NN.
  #' @param quiet passes through as -quiet.
  public = list(
    initialize = function(inset = NULL, mask = NULL, threshold_file = NULL, data_file = NULL, mask_from_hdr = NULL, out_mask = NULL, 
      ithr = NULL, idat = NULL, onesided = NULL, twosided = NULL, bisided = NULL, 
      lower_thresh = NULL, upper_thresh = NULL, one_thresh = NULL, one_tail = NULL,
      NN = NULL, clust_nvox = NULL, clust_nvol = NULL, pref_map = NULL, pref_dat = NULL,
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
        orient <- toupper(orient)
        checkmate::assert_subset(orient, c("LPI", "RAI"))
        # should probably check LPI, RAI, etc.
        private$pvt_orient <- orient
      }

      if (!is.null(clusterize_output_file)) {
        checkmate::assert_string(clusterize_output_file)
        private$pvt_clusterize_output_file <- clusterize_output_file
      }

      # sort out the -inset to be used in the 3dClusterize call (if a threshold and data file are passed)
      if (!is.null(private$pvt_threshold_file)) {
        default_wd <- dirname(private$pvt_threshold_file)
        if (is.null(private$pvt_data_file)) {
          private$pvt_input_file <- private$pvt_threshold_file # 3D input to clusterize
        } else {
          # build tmp dataset that combined threshold and data file
          private$pvt_input_file <- tempfile(pattern = "clustinput", fileext = ".nii.gz")
          private$pvt_combine_call <- glue("3dTcat -output {private$pvt_input_file} {private$pvt_threshold_file} {private$pvt_data_file}")
        }
      } else {
        default_wd <- dirname(private$pvt_inset_file)
        private$pvt_input_file <- private$pvt_inset_file
      }

      # if output files are provided without any path specification, use the directory of the input file
      private$pvt_clusterize_output_file <- R.utils::getAbsolutePath(private$pvt_clusterize_output_file,
        workDirectory = default_wd
      )

      if (!is.null(private$pvt_pref_map)) {
        private$pvt_pref_map <- R.utils::getAbsolutePath(private$pvt_pref_map, workDirectory = default_wd)
      }

      if (!is.null(private$pvt_pref_dat)) {
        private$pvt_pref_dat <- R.utils::getAbsolutePath(private$pvt_pref_dat, workDirectory = default_wd)
      }

    },

    #' @description run the 3dClusterize command relevant to this object
    run = function() {
      private$build_call()
      run_afni_command(private$pvt_clusterize_call)
    },

    #' @description return the 3dClusterize table of clusters as a data.frame
    #' @details This function will return an empty data.frame if the 3dClusterize output file cannot be found.
    get_clust_df = function() {
      if (!checkmate::test_file_exists(private$pvt_clusterize_output_file)) {
        warning(glue("The expected 3dClusterize output does not exist: {private$pvt_clusterize_output_file}"))
        return(invisible(data.frame()))
      }

      clust_df <- read.table(private$pvt_clusterize_output_file)
      clust_txt <- readLines(private$pvt_clusterize_output_file)
      comment_lines <- grep("^\\s*#", clust_txt, perl = TRUE)
      table_lines <- grep("^\\s*#", clust_txt, perl = TRUE, invert = TRUE)

      # header
      cols <- clust_txt[min(table_lines) - 2L] # header is 2 lines before table starts
      cols <- make.names(strsplit(sub("^#", "", cols), "\\s{2,}")[[1L]]) # columns are separated by 2+ spaces
      stopifnot(length(cols) == ncol(clust_df))
      names(clust_df) <- cols

      clust_df$roi_num <- seq_len(nrow(clust_df))
      clust_df <- clust_df %>% dplyr::select(roi_num, everything()) # place roi_num first

      parse_header_rows <- function(lines, clust_df) {
        attrib <- lapply(lines, function(x) {
          # don't allow errant lines that lack the comment #; also only parse lines with key - value pairing (separated by =)
          if (!grepl("^#", x) || !grepl("=", x, fixed=TRUE)) {
            return(NULL)
          } else {
            keyval_split <- strsplit(x, "=", fixed = TRUE)[[1L]]
            if (length(keyval_split) != 2) {
              warning("Cannot sort out how to parse 3dClusterize header line: ", x)
              return(NULL)
            } else {
              key <- trimws(sub("^\\s*#\\s*\\[\\s*(.+)", "\\1", keyval_split[[1]], perl = TRUE))
              value <- trimws(sub("\\s*(.*)\\]\\s*$", "\\1", keyval_split[[2]], perl = TRUE))
              return(value %>% setNames(key))
            }
          }
        })

        # add each attribute to the data.frame
        for (aa in attrib) {
          if (!is.null(aa)) {
            attr(clust_df, names(aa)) <- aa
          }
        }

        #attr(clust_df, key) <- value
        return(clust_df)
      }

      comment_txt <- clust_txt[comment_lines]

      # There are a couple of weird lines for threshold value and nvoxel threshold. These allow for multiple = signs in the values.
      # Here, we pre-parse these lines and make them into single key-value pairs per line by adding one line for each.
      eq_lens <- sapply(gregexpr("=", text = comment_txt, fixed = TRUE), length)
      if (any(eq_lens > 1L)) {
        mult_eq_txt <- comment_txt[eq_lens > 1]

        # remove these lines from comment_txt (will be replaced momentarily by expanded copies)
        comment_txt <- comment_txt[eq_lens == 1L]
        for (line in mult_eq_txt) {
          firsteq <- gregexpr("=", text = line, fixed = TRUE)[[1]][1]
          stem <- trimws(substr(line, 1, firsteq - 1))
          values <- substr(line, firsteq + 1, nchar(line))
          values_split <- trimws(strsplit(values, ";", fixed = TRUE)[[1]])
          has_subkeys <- all(grepl("=", values_split, fixed = TRUE)) # TRUE for lines like left-tail thr=-3.000;  right-tail thr=3.000
          if (isTRUE(has_subkeys)) {
            sub_split <- strsplit(values_split, "=", fixed=TRUE)
            for (ii in seq_along(sub_split)) {
              # combine overall stem with subkey-specific stem, separated by comma. Add " ]" for all items except the last (has it already)
              comment_txt <- c(
                comment_txt,
                paste0(stem, ", ", sub_split[[ii]][1], " = ", sub_split[[ii]][2], ifelse(ii == length(sub_split), "", " ]"))
              )
            }
          } else {
            # really just two lines have been piled into one, like #[ Nvoxel threshold    = 35;  Volume threshold = 425.845 ]
            # split these into two complete lines
            for (ii in seq_along(values_split)) {
              if (ii == 1L) { # for first split, key has been pulled into stem
                comment_txt <- c(comment_txt, paste0(stem, " = ", values_split[ii], " ]"))
              } else { # for other splits, key is still to the left of the equals sign
                comment_txt <- c(comment_txt, paste0("#[ ", values_split[ii]))
              }
            }
          }
        }
      }

      # add header row information as attributes to data.frame
      clust_df <- parse_header_rows(comment_txt, clust_df)
      return(clust_df)

    },

    #' @description method to read and return the integer-valued clusterized mask (aka -pref_map) as an oro.nifti object
    #' @return an oro.nifti object containing the clusterized mask from 3dClusterize
    get_cluster_map_nifti = function() {
      if (checkmate::test_file_exists(private$pvt_pref_map)) {
        oro.nifti::readNIfTI(private$pvt_pref_map, reorient=FALSE)
      } else {
        NULL # return NULL if file is missing
      }
    },

    #' @description returns the 3dClusterize call for this specification
    get_call = function() {
      private$build_call()
      return(private$pvt_clusterize_call)
    },

    #' @description Provides a vector of expected output files that correspond to this 3dClusterize setup
    #' @param exclude_missing if TRUE (default), any output file that cannot be found will be returned as NA.
    #' @return a named vector of output files related to this 3dClusterize setup
    get_output_files = function(exclude_missing = TRUE) {

      cluster_map <- private$pvt_pref_map
      cluster_masked_data <- private$pvt_pref_dat
      cluster_table <- private$pvt_clusterize_output_file

      if (isTRUE(exclude_missing)) {
        if (!checkmate::test_file_exists(cluster_map)) cluster_map <- NA_character_
        if (!checkmate::test_file_exists(cluster_masked_data)) cluster_masked_data <- NA_character_
        if (!checkmate::test_file_exists(cluster_table)) cluster_table <- NA_character_
      }

      named_vector(cluster_table, cluster_map, cluster_masked_data)

    },

    #' @description returns the orientation code for this 3dClusterize call (LPI or RAI)
    get_orient = function() {
      private$pvt_orient
    },

    #' @description Add's an afni_whereami object to this class in the $whereami slot. The corresponding
    #'   whereami command is also run when this is added so that coordinates and labels can be obtained
    #'   immediately. To access the whereami object and its methods, use $whereami()
    add_whereami = function() {
      if (!is.null(private$pvt_whereami) && inherits(private$pvt_whereami, "afni_whereami")) {
        message("whereami object already added to this clusterize object. Use $whereami() to access it.")
      } else if (!self$is_complete()) {
        message("Cannot add whereami to 3dclusterize object until clusterization is run!")
      } else {
        private$pvt_whereami <- afni_whereami$new(
          afni_3dclusterize_obj = self
        )
        private$pvt_whereami$run()
      }
    },

    #' @description returns TRUE if all expected output files exist for this 3dClusterize call
    is_complete = function() {
      expect_files <- self$get_output_files(exclude_missing = FALSE)
      all(sapply(expect_files, checkmate::test_file_exists))
    },

    #' @description passthrough access to whereami object if that has been 
    whereami = function() { # expose nested object
      if (is.null(private$pvt_whereami)) {
        message("No whereami has been added to this object yet. Use $add_whereami() to do so.")
        return(invisible(NULL))
      } else {
        private$pvt_whereami
      }      
    }
  )
)

# x <- afni_3dclusterize$new(
#   threshold_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-pe/L2m-l2_l2c-emotion.happy/L3m-age_sex/FEAT_l1c-EV_pe.gfeat/cope1.feat/stats/zstat6.nii.gz",
#   lower_thresh = -3, upper_thresh = 3, bisided = TRUE, NN = 1, clust_nvox = 35, pref_map = "zstat_clusterize.nii.gz"
# )

# x$get_call()
# x$is_complete()
# vv <- x$get_clust_df()
# x$add_whereami()
# x$whereami()$get_whereami_df()

# x$get_output_files()
# yy <- x$get_cluster_map_nifti()

