#' This is a small helper function that returns the inputs provided in the feat_files field for
#' a set of .fsf files.
#' 
#' @details 
#' One can also pass in .feat or .gfeat directories and the function will
#' use the design.fsf files within each of these to find the inputs
#' 
#' @param feat_list a vector of design.fsf filenames, and/or .feat/.gfeat directory names
#' @param recursive a boolean indicating whether to drill down and find lower-level inputs for .gfeat/.feat inputs
#' @author Michael Hallquist
#' @importFrom checkmate assert_file_exists assert_directory_exists
#' @export
read_feat_inputs <- function(feat_list, recursive=FALSE) {
  is_feat_dir <- grepl(".g?feat$", feat_list, perl=TRUE)

  if (any(!is_feat_dir)) {
    check_fsfs <- feat_list[!is_feat_dir]
    sapply(check_fsfs, checkmate::assert_file_exists)
    if (!all(goodfiles <- grepl(".fsf$", check_fsfs, perl=TRUE))) {
      error("The following inputs do not end in .fsf, .gfeat, or .feat: ", paste(check_fsfs[!goodfiles], collapse=", "))
    }
  }

  #add the /design.fsf suffix to directory inputs
  if (any(is_feat_dir)) { 
    sapply(feat_list[is_feat_dir], checkmate::assert_directory_exists)
    feat_list[is_feat_dir] <- file.path(feat_list[is_feat_dir], "design.fsf") 
  }

  inputs <- lapply(feat_list, function(ff) {
    fsf <- readLines(ff)
    infile <- grep("^set feat_files\\(\\d+\\)", fsf, perl=TRUE, value=TRUE)
    #input <- paste0(sub("set feat_files\\(d+\\) \"([^\"]+)\"", "\\1", nifti, perl=TRUE), ".nii.gz")
    fnumbers <- as.numeric(sub("set feat_files\\((\\d+)\\) \"[^\"]+\"", "\\1", infile, perl=TRUE)) #always return in sorted order
    infile <- sub("set feat_files\\(\\d+\\) \"([^\"]+)\"", "\\1", infile, perl=TRUE)
    return(infile[fnumbers])
  })

  inputs_are_featdirs <- grepl(".g?feat$", inputs, perl=TRUE)
  if (isTRUE(recursive) && any(inputs_are_featdirs)) {
    inputs <- lapply(1:length(inputs), function(ii) {
      if (isTRUE(inputs_are_featdirs[ii])) {
        return(read_feat_inputs(inputs[ii], recursive=recursive))
      } else {
        return(list(inputs[ii])) #single-element list
      }
    })
  }

  #need some sort of wrap-up function here in the recursive case
  return(inputs)
}


#' helper function to look at whether feat ingredients exist for a .feat/.gfeat directory
#' 
#' @param feat_dir an expected .feat/.gfeat directory for an analysis
#' @param feat_fsf an expected .fsf file corresponding to the feat analysis
#' @param lg an optional lgr logger object to be used for logging. If not passed, the root
#'   logger will be used
#'
#' @return a data.frame containing information about whether the .feat analysis is complete
#'   and whether various ingredients are present
#'
#' @importFrom anytime anytime
#' @keywords internal
get_feat_status <- function(feat_dir, feat_fsf=NULL, lg=NULL, prefix=NULL) {
  checkmate::assert_string(feat_dir)
  checkmate::assert_string(feat_fsf, null.ok = TRUE)
  if (is.null(lg)) lg <- lgr::get_logger() #just use root logger

  # handle feat_dir with no extension
  if (!grepl(".feat$", feat_dir)) {
    feat_dir <- paste0(feat_dir, ".feat")
    if (!dir.exists(feat_dir) && dir.exists(sub(".feat$", ".gfeat", feat_dir))) {
      feat_dir <- sub(".feat$", ".gfeat", feat_dir)
    }
  }

  if (is.null(feat_fsf)) {
    # assume that fsf of same name as feat_dir exists at same level of filesystem hierarchy
    feat_fsf <- sub("\\.g?feat$", ".fsf", feat_dir)
  }

  if (dir.exists(feat_dir) && !file.exists(feat_fsf) && file.exists(file.path(feat_dir, "design.fsf"))) {
    # fall back to design.fsf inside the .feat folder if .fsf is not present in parent folder
    feat_fsf <- file.path(feat_dir, "design.fsf")
  }

  feat_checks <- list()
  feat_checks$feat_fsf <- feat_fsf
  feat_checks$feat_fsf_modified_date <- if (file.exists(feat_fsf)) file.info(feat_fsf)$mtime else as.POSIXct(NA)
  feat_checks$feat_fsf_exists <- file.exists(feat_fsf)
  feat_checks$feat_dir <- feat_dir
  feat_checks$feat_dir_exists <- dir.exists(feat_dir)
  feat_checks$feat_execution_start <- as.POSIXct(NA)
  feat_checks$feat_execution_end <- as.POSIXct(NA)
  feat_checks$feat_execution_min <- NA_real_
  feat_checks$feat_complete <- FALSE # default: FALSE for anything but a .feat_complete outcome
  feat_checks$feat_failed <- NA # default: NA if feat hasn't even been run

  if (dir.exists(feat_dir)) {
    if (file.exists(file.path(feat_dir, ".feat_complete"))) {
      lg$debug("Feat directory is complete: %s", feat_dir)
      timing_file <- file.path(feat_dir, ".feat_complete")
      if (file.exists(file.path(feat_dir, ".feat_fail"))) {
        lg$warn("Both .feat_complete and .feat_fail objects exist in %s", feat_dir)
        lg$warn("Assuming that .feat_complete reflects a successful completion of feat")
      }
      feat_checks$feat_complete <- TRUE
      feat_checks$feat_failed <- FALSE
    } else if (file.exists(file.path(feat_dir, ".feat_fail"))) {
      lg$debug("Detected feat failure in: %s", feat_dir)
      timing_file <- file.path(feat_dir, ".feat_fail")
      feat_checks$feat_failed <- TRUE
    } else {
      timing_file <- NULL
    }

    if (!is.null(timing_file)) {
      timing <- readLines(timing_file)
       if (length(timing) > 0L) {
         # convert to POSIXct object to allow for any date calculations
         timing <- anytime::anytime(timing)
         feat_checks$feat_execution_start <- timing[1L]
         if (length(timing) == 2L) {
           feat_checks$feat_execution_end <- timing[2L]
           feat_checks$feat_execution_min <- as.numeric(difftime(timing[2L], timing[1L], units = "mins"))
         } else {
           lg$warn("Did not find two timing entries in %s.", timing_file)
           lg$warn("File contents: %s", timing)
         }
       }
    }
  }

  df <- as.data.frame(feat_checks)
  if (!is.null(prefix)) { #support naming prefix like "l1_"
    names(df) <- paste0(prefix, names(df))
  }
  return(df)
}

#' Helper function to find/replace or insert additional arguments to include
#'   in FSF syntax for a FEAT model.
#'
#' @param fsf_syntax a character vector containing all fsf syntax
#' @param feat_args a named list 
#' @return a modified version of \code{fsf_syntax} that now includes all custom
#'   arguments from \code{feat_args}
#'
#' @details If a tag already exists in the fsf syntax, it is replaced by the custom argument
#'   in \code{feat_args}. If the tag does not exist, it is inserted.
#' @keywords internal
add_custom_feat_syntax <- function(fsf_syntax, feat_args, lg=NULL) {

  checkmate::assert_character(fsf_syntax)
  checkmate::assert_list(feat_args, null.ok = TRUE)

  if (is.null(feat_args)) {
    #no custom arguments, return fsf_syntax unchanged
    return(fsf_syntax)
  }

  if (is.null(lg)) {
    lg <- lgr::get_logger() #use root logger if not passed
  }

  # for now, this only supports single name/value pairs, even though the list may contain
  # things like list(arg=c(1,2,3,4)). The function skips these for now.
  # TODO: Provide better enforcement of expected structure in finalize_pipeline_configuration
  # or setup_glm_pipeline so that user rectifies things at setup.

  # handle additional custom feat level 1 fields in fsf syntax

  for (ii in seq_along(feat_args)) {
    this_name <- names(feat_args)[ii]
    this_value <- feat_args[[ii]]
    if (length(this_value) > 1L) {
      lg$warn("feat_args: Cannot handle settings with multiple values: %s. Skipping out.", this_name)
      next
    } else {
      lg$debug("Adding custom feat setting: %s = %s", this_name, this_value)
      if (any(grepl(paste0("set fmri(", this_name, ")"), fsf_syntax, fixed = TRUE))) {
        lg$debug("Substituting existing value of feat setting: %s", this_name)
        fsf_syntax <- gsub(paste0("(set fmri\\s*\\(", this_name, "\\))\\s*(.*)"), paste0("\\1 ", this_value),
          fsf_syntax,
          perl = TRUE
        )
      } else {
        fsf_syntax <- c(fsf_syntax, paste0("set fmri(", this_name, ") ", this_value))
      }
    }
  }

  return(fsf_syntax)

}

#' Helper function to look at feat outputs and files to determine if execution is complete
#'
#' @param gpa a \code{glm_pipeline_arguments object}
#' @param level the level of analysis to be refreshed (1, 2, or 3)
#' @param lg an optional lgr logger object used for logging
#' @return a modified copy of \code{gpa} with the feat columns of
#'   $l1_model_setup, $l2_model_setup, or $l3_model_setup refreshed
#'
#' @export
#' @importFrom dplyr select
#' @importFrom purrr pmap_dfr
refresh_feat_status <- function(gpa, level = 1L, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)
  setup_name <- paste0("l", level, "_model_setup")

  if ("fsl" %in% gpa$glm_software && setup_name %in% names(gpa)) {
    if (is.null(lg)) lg <- lgr::get_logger()
    orig <- gpa[[setup_name]]$fsl
    if (is.null(orig) || (is.data.frame(orig) && nrow(orig) == 0L)) {
      lg$warn("Could not find populated $fsl object in gpa$%s$fsl", setup_name)
    } else if (!"feat_fsf" %in% names(orig)) {
      lg$warn("No $feat_fsf field in gpa$%s$fsl", setup_name)
    } else {
      lg$info("Found existing %s field. Refreshing status of L%d feat execution and outputs.", setup_name, level)
      refresh <- orig %>%
        dplyr::select(feat_dir, feat_fsf) %>%
        purrr::pmap_dfr(get_feat_status, lg = lg)

      # copy back relevant columns into data structure
      gpa[[setup_name]]$fsl[, names(refresh)] <- refresh
    }
  }

  return(gpa)
}

#' helper function to look up core stats outputs from a .gfeat folder
#' @param gfeat_dir a .gfeat folder containing the outputs of an FSL analysis
#' @param what What to parse in each folder. Currently just passed through to read_feat_dir
#' @return a list containing sorted vectors of each stat output
#' @importFrom glue glue
#' @importFrom checkmate assert_directory_exists assert_file_exists test_directory_exists
#' @importFrom readr read_delim
#' @export
read_gfeat_dir <- function(gfeat_dir, what = "all") {
  gfeat_dir <- normalizePath(gfeat_dir) # convert to absolute path

  mask_file <- file.path(gfeat_dir, "mask.nii.gz")
  mask_file <- ifelse(checkmate::test_file_exists(mask_file), mask_file, NA_character_)

  cope_dirs <- grep("/cope[0-9]+\\.feat",
    list.dirs(path = gfeat_dir, full.names = TRUE, recursive = FALSE),
    value = TRUE, perl = TRUE
  )

  # sort copes in numeric order to match ascending expectation
  cope_nums <- as.numeric(sub(".*/cope(\\d+).feat", "\\1", cope_dirs, perl = TRUE))
  cope_dirs <- cope_dirs[order(cope_nums)]

  cope_list <- lapply(cope_dirs, read_feat_dir, what = what)
  names(cope_list) <- basename(cope_dirs)
  ret_list <- list(cope_dirs = cope_list, design_files = get_design_files(gfeat_dir), mask_file = mask_file)
  class(ret_list) <- c("list", "gfeat_info")
  return(ret_list)
}

# incomplete function stub for subselecting parts of the gfeat_info into a combined data.frame
# initially developed for getting specific cope dirs and zstats for ptfce setup. Low priority at present
get_feat_files <- function(gfeat_info, lower_cope_numbers = NULL, what = "") {
  checkmate::assert_class(gfeat_info)
  checkmate::assert_integerish(lower_cope_numbers)
  if (is.null(gfeat_info$cope_dir)) {
    stop("Cannot find cope_dirs in gfeat_info")
  }

  # look for specific lower-level cope .feat directories that are in the .gfeat folder
  if (!is.null(lower_cope_numbers)) {
    expect_cope_dirs <- paste0("cope", lower_cope_numbers, ".feat")
  } else {
    expect_cope_dirs <- names(gfeat_info$cope_dirs)
  }

  checkmate::assert_subset(expect_cope_dirs, names(gfeat_info$cope_dirs))
  #file_list <- c
}

get_design_files <- function(dir) {
  sapply(
    c(
      "design.con", "design_cov.png", "design_cov.ppm", "design.frf",
      "design.fsf", "design.mat", "design.min", "design.png", "design.ppm",
      "design.grp", "design.lcon", "design.lev"
    ), function(x) {
      ff <- file.path(dir, x)
      if (file.exists(ff)) {
        return(ff)
      } else {
        return(character(0))
      }
    },
    simplify = FALSE
  )
}

#' helper function to look up core stats outputs from a .feat folder
#' @param feat_dir a .feat folder containing the outputs of an FSL analysis
#' @param what which parts of the .feat folder should be parsed. If "all", extract everything.
#' @return a list containing sorted vectors of each stat output
#' @importFrom glue glue
#' @importFrom checkmate assert_directory_exists assert_file_exists test_directory_exists
#' @importFrom readr read_delim
#' @importFrom gsubfn strapply
#' @export
read_feat_dir <- function(feat_dir, what = "all") {
  if (!checkmate::test_directory_exists(feat_dir)) return(NULL)
  feat_dir <- normalizePath(feat_dir) # convert to absolute path
  stats_dir <- file.path(feat_dir, "stats")
  if (!checkmate::test_directory_exists(stats_dir)) return(NULL)

  checkmate::assert_character(what)
  if (what[1L] == "all") {
    what <- c("stat_files", "aux_files", "parsed_txt")
  }

  # inside the stats directory we will have pes, copes, varcopes, and zstats
  z_files <- list.files(path = stats_dir, pattern = "zstat[0-9]+\\.nii.*", full.names = TRUE)
  z_nums <- as.numeric(sub(".*zstat(\\d+)\\.nii.*$", "\\1", z_files, perl = TRUE))
  z_files <- z_files[order(z_nums)] # need to sort stat files numerically for labeling to work appropriately

  t_files <- file.path(dirname(z_files), sub("zstat", "tstat", basename(z_files), fixed = TRUE))
  cope_files <- file.path(dirname(z_files), sub("zstat", "cope", basename(z_files), fixed = TRUE))
  varcope_files <- file.path(dirname(z_files), sub("zstat", "varcope", basename(z_files), fixed = TRUE))
  zptfce_files <- file.path(dirname(z_files), sub("zstat(\\d+)", "zstat\\1_ptfce", basename(z_files), perl = TRUE))
  zptfce_exists <- sapply(zptfce_files, checkmate::test_file_exists)
  zptfce_files[!zptfce_exists] <- NA_character_ # NA out the missing ptfce files

  checkmate::assert_file_exists(t_files)
  checkmate::assert_file_exists(cope_files)
  checkmate::assert_file_exists(varcope_files)

  #thresh files live in parent folder
  zthresh_files <- file.path(feat_dir, sub("zstat", "thresh_zstat", basename(z_files), fixed = TRUE))
  rendered_zthresh_files <- file.path(feat_dir, sub("zstat", "rendered_thresh_zstat", basename(z_files), fixed = TRUE))

  # parameter estimate (PE) images
  pe_files <- list.files(path = stats_dir, pattern = "^pe[0-9]+\\.nii.*", full.names = TRUE)
  pe_nums <- as.numeric(sub(".*pe(\\d+)\\.nii.*$", "\\1", pe_files, perl=TRUE))
  pe_files <- pe_files[order(pe_nums)]

  get_file <- function(str) {
    ff <- file.path(feat_dir, str)
    ff <- ifelse(checkmate::test_file_exists(ff), ff, NA_character_)
  }

  mask_file <- get_file("mask.nii.gz")
  smoothness_file <- get_file("stats/smoothness")
  dof_file <- get_file("stats/dof")

  aux_files <- sapply(
    c(
      "sigmasquareds", "threshac1", "res4d",
      "mean_outlier_random_effects_var1", "mean_random_effects_var1", "prob_outlier1", "global_prob_outlier1"
    ),
    function(x) {
      ff <- list.files(path = stats_dir, pattern = glue("^{x}\\.nii.*"), full.names = TRUE)
      if (length(ff) > 0L) {
        return(ff)
      } else {
        return(character(0))
      }
    }
  )

  parsed_txt <- sapply(c("dof", "smoothness", "lmax_zstat", "cluster_zstat"), simplify = FALSE, FUN = function(x) {
    if (!"parsed_txt" %in% what) {
      return(list()) # return empty list if parsing not requested
    }

    if (x == "lmax_zstat") {
      pat <- "lmax_zstat\\d+_std\\.txt"
    } else if (x == "cluster_zstat") {
      pat <- "cluster_zstat\\d+_std\\.txt"
    } else {
      pat <- x # literal
    }
    files <- list.files(path = feat_dir, pattern = glue("^{pat}$"), full.names = TRUE, recursive=TRUE)
    if (length(files) > 0L) {
      if (x == "dof") {
        stopifnot(length(files) == 1L)
        return(as.numeric(readLines(files)))
      } else if (x == "smoothness") {
        stopifnot(length(files) == 1L)
        ss <- readLines(files)
        dlh <- NA_real_
        volume <- NA_integer_
        resels <- NA_real_
        fwhm_voxel <- rep(NA_real_, 3)
        fwhm_mm <- rep(NA_real_, 3)
        dlh_line <- grep("^\\s*DLH\\s+", ss, value = TRUE)
        if (length(dlh_line) == 1L && nchar(dlh_line[1L]) > 0L) {
          dlh <- as.numeric(sub("^\\s*DLH\\s+(-?[0-9.]+)\\s*$", "\\1", dlh_line[1L], perl = TRUE))
        }
        volume_line <- grep("^\\s*VOLUME\\s+", ss, value = TRUE)
        if (length(volume_line) == 1L && nchar(volume_line[1L]) > 0L) {
          volume <- as.integer(sub("^\\s*VOLUME\\s+([0-9]+)\\s*$", "\\1", volume_line[1L], perl = TRUE))
        }
        resels_line <- grep("^\\s*RESELS\\s+", ss, value = TRUE)
        if (length(volume_line) == 1L && nchar(volume_line[1L]) > 0L) {
          resels <- as.numeric(sub("^\\s*RESELS\\s+([0-9.]+)\\s*$", "\\1", resels_line[1L], perl = TRUE))
        }
        fwhm_voxel_line <- grep("^\\s*FWHMvoxel\\s+", ss, value = TRUE)
        if (length(fwhm_voxel_line) == 1L && nchar(fwhm_voxel_line[1L]) > 0L) {
          fwhm_voxel <- as.numeric(strapply(fwhm_voxel_line, "[0-9.]+", simplify = c))
          stopifnot(length(fwhm_voxel) == 3L)
        }
        fwhm_mm_line <- grep("^\\s*FWHMmm\\s+", ss, value = TRUE)
        if (length(fwhm_mm_line) == 1L && nchar(fwhm_mm_line[1L]) > 0L) {
          fwhm_mm <- as.numeric(strapply(fwhm_mm_line, "[0-9.]+", simplify = c))
          stopifnot(length(fwhm_mm) == 3L)
        }

        return(named_list(dlh, volume, resels, fwhm_voxel, fwhm_mm))
      } else if (x == "cluster_zstat") {
        f_out <- lapply(files, read.table, header = T, sep = "\t")
        names(f_out) <- basename(files)
        return(f_out)
      } else if (x == "lmax_zstat") {
        f_out <- lapply(files, function(f) {
          inp <- readLines(f)
          # these files have a trailing tab on the first line and also have a space in the header row
          inp[1] <- trimws(sub("Cluster Index", "Cluster_Index", inp[1], fixed = TRUE))
          res <- readr::read_delim(I(inp), delim = "\t", show_col_types = FALSE)
        })
        names(f_out) <- basename(files)
        return(f_out)
      } else {
        return(lapply(files, readLines))
      }
    } else {
      return(character(0))
    }
  })

  # get design files
  design_files <- get_design_files(feat_dir)

  # lookup contrast names
  contrast_file <- design_files$design.con
  if (checkmate::test_file_exists(contrast_file)) {
    # design.con contains names of contrasts
    dcon <- readLines(contrast_file)
    contrast_names <- sub("/ContrastName\\d+\\s+([\\w_.]+).*", "\\1", grep("/ContrastName", dcon, value = TRUE), perl = TRUE)
    contrast_names <- gsub("\\s", "_", contrast_names, perl = TRUE) # replace spaces with underscores to make labels accurate in AFNI
  } else {
    contrast_names <- NULL
  }

  ret_list <- named_list(mask_file, smoothness_file, dof_file,
    pe_files, cope_files, varcope_files, z_files, zthresh_files, zptfce_files,
    rendered_zthresh_files, t_files, contrast_names, aux_files, design_files, parsed_txt
  )

  # create a combined data.frame of cope-level statistics
  expect_length <- length(cope_files)

  # these elements should all have the same length (one file per contrast)
  to_stitch <- ret_list[c(
    "contrast_names", "cope_files", "varcope_files", "z_files", "zthresh_files",
    "rendered_zthresh_files", "t_files", "zptfce_files"
  )]

  to_stitch <- do.call(data.frame, lapply(to_stitch, function(x) {
    if (length(x) != expect_length) {
      return(rep(NA_character_, expect_length))
    } else {
      return(x)
    }
  })) %>% setNames(c("contrast_name", "cope", "varcope", "z", "zthresh", "rendered_zthresh", "t", "zptfce"))

  to_stitch <- to_stitch %>%
    dplyr::mutate(cope_number = 1:n()) %>%
    dplyr::select(cope_number, contrast_name, everything())

  ret_list$cope_df <- to_stitch

  return(ret_list)

}

gfeat_stats_to_brik <- function(gfeat_dir, cope_names = NULL, level=NULL) {
  checkmate::assert_directory_exists(gfeat_dir)
  gfeat_dir <- normalizePath(gfeat_dir) # convert to absolute path
  checkmate::assert_character(cope_names, null.ok=TRUE)
  checkmate::assert_integerish(level, lower=1L, upper=3L, null.ok=TRUE)

  cope_dirs <- copedirs <- grep("/cope[0-9]+\\.feat",
    list.dirs(path = gfeat_dir, full.names = TRUE, recursive = FALSE),
    value = TRUE, perl = TRUE
  )

  # sort copes in numeric order to match ascending expectation
  cope_nums <- as.numeric(sub(".*/cope(\\d+).feat", "\\1", cope_dirs, perl = TRUE))
  cope_dirs <- copedirs[order(cope_nums)]

  # use external names if provided
  if (!is.null(cope_names)) names(cope_dirs) <- cope_names

  lower_img_list <- list()
  for (d in seq_along(cope_dirs)) {
    ## run the L1 -> AFNI conversion for each separate cope
    cat("Processing: ", cope_dirs[d], "\n\n")

    lower_img <- feat_stats_to_brik(cope_dirs[d], what = c("cope_files", "z_files"))

    # if .feat folders are used as the input to a .gfeat analysis, then the lower level
    # contrast name will be placed in design.lev. If not, it can be passed in as cope_names
    if (isTRUE(checkmate::test_file_exists(lower_img$feat_info$design_files$design.lev))) {
      lower_name <- readLines(lower_img$feat_info$design_files$design.lev)
      if (nchar(lower_name) == 0L) lower_name <- NULL
      names(cope_dirs)[d] <- lower_name
    }

    lower_img$lower_brik_name <- names(cope_dirs)[d]
    if (!is.null(level)) lower_img$level <- level - 1 # for combining later

    gfeat_cope_inputs <- file.path(cope_dirs[d], "filtered_func_data.nii.gz")
    if (checkmate::test_file_exists(gfeat_cope_inputs)) lower_img$gfeat_cope_inputs <- gfeat_cope_inputs

    gfeat_varcope_inputs <- file.path(cope_dirs[d], "var_filtered_func_data.nii.gz")
    if (checkmate::test_file_exists(gfeat_varcope_inputs)) lower_img$gfeat_varcope_inputs <- gfeat_varcope_inputs


    lower_img_list[[d]] <- lower_img
  }

  names(lower_img_list) <- basename(cope_dirs)

  return(list(gfeat_dir = gfeat_dir, level = level, lower_img_list = lower_img_list))

}

#' Function to combine Feat L3 analyses into AFNI BRIK+HEAD files according to user-specified combinations
#' @param gpa a glm_pipeline_arguments object having a populated l3_model_setup field. The L3 analyses should also
#'   be complete (i.e., after run_glm_pipeline).
#' @param feat_l3_combined_filename a glue expression for the path and filename prefix
#' @param feat_l3_combined_briknames a glue expression for naming the subbriks in the AFNI output
#' @param template_brain an optional filename for the MNI template that should be used as an underlay in AFNI. This
#'   image will be symbolically linked into each directory created by \code{combine_feat_l3_to_afni}.
#' @details To specify the folder and filename structure for the combined feat analyses, use a \code{glue} expression
#'   that indicates how outputs should be structured. In particular, variables in gpa$l3_model_setup$fsl can be used for
#'   dynamically naming of afni outputs.
#'   
#' \describe{
#'   \item{l1_model}{the name of the level 1 model}
#'   \item{l2_model}{the name of the level 2 model}
#'   \item{l3_model}{the name of the level 3 model}
#'   \item{l1_cope_name}{the level 1 contrast name}
#'   \item{l2_cope_name}{the level 2 contrast name}
#'   \item{l3_cope_name}{the level 3 contrast name}
#' }
#' 
#' @examples
#' \dontrun{
#'   feat_l3_combined_filename <- "{gpa$output_directory}/afni_combined/L1m-{l1_model}/l1c-{l1_cope_name}/L3m-{l3_model}_stats"
#'   feat_l3_combined_briknames = "l2c-{l2_cope_name}_l3c-{l3_cope_name}"
#'   template_brain <- "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_brain.nii"
#' 
#'   combine_feat_l3_to_afni(gpa, feat_l3_combined_filename, feat_l3_combined_briknames, template_brain)
#' }
#'
#' @export
#' @importFrom tidyr unnest pivot_longer
#' @importFrom glue glue_data
combine_feat_l3_to_afni <- function(gpa, feat_l3_combined_filename=NULL, feat_l3_combined_briknames=NULL, template_brain=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (!is.null(template_brain)) {
    checkmate::assert_file_exists(template_brain)
    template_ext <- sub(".*(\\.nii(\\.gz)*)", "\\1", template_brain, perl=TRUE)
  }

  if (is.null(feat_l3_combined_filename)) {
    feat_l3_combined_filename <- gpa$output_locations$feat_l3_combined_filename
  }

  if (is.null(feat_l3_combined_briknames)) {
    feat_l3_combined_briknames <- gpa$output_locations$feat_l3_combined_briknames
  }

  checkmate::assert_string(feat_l3_combined_filename)
  checkmate::assert_string(feat_l3_combined_briknames)

  # run the gfeat -> afni conversion for each cope*.feat directory in the corresponding l3 .gfeat folder
  # l3_stats <- lapply(gpa$l3_model_setup$fsl$feat_dir, gfeat_stats_to_brik, level = 3)

  # Get a list of information about each cope1.feat directory in the .gfeat folder.
  # Note: in the current implementation of l3 analysis, the cope1.feat folder in the l3 .gfeat folder always corresponds
  # to the single lower-level cope that was fed forward for group analysis. This relates to the 'Inputs are lower-level
  # FEAT directories' (many .feat folders in .gfeat) versus 'Inputs are 3D cope images from FEAT directores' (one .feat folder).
  l3_stats <- lapply(file.path(gpa$l3_model_setup$fsl$feat_dir, "cope1.feat"), read_feat_dir, what = "stat_files") #don't parse text files

  if (isTRUE(gpa$multi_run)) {
    meta_df <- gpa$l3_model_setup$fsl %>% dplyr::select(l1_model, l1_cope_name, l2_model, l2_cope_name, l3_model, feat_dir)
  } else {
    meta_df <- gpa$l3_model_setup$fsl %>% dplyr::select(l1_model, l1_cope_name, l3_model, feat_dir)
  }

  # note that keep_empty = FALSE will drop any .feat folder that failed to be parsed (usually a failed run)
  meta_df$cope_df <- lapply(l3_stats, "[[", "cope_df") # add a list-column, then unnest to expand
  meta_df <- meta_df %>% unnest(cope_df, keep_empty = FALSE) %>% # expand so that multiple rows of copes in cope_df are added
    dplyr::rename(l3_cope_name=contrast_name, l3_cope_number=cope_number)

  # use glue_data to evaluate the glue file and brik expressions for every row of the data.frame
  meta_df <- meta_df %>%
    mutate(afni_out=glue_data(., !!feat_l3_combined_filename), afni_briks=glue_data(., !!feat_l3_combined_briknames))

  meta_split <- split(meta_df, meta_df$afni_out)

  lapply(meta_split, function(ss) {
    ss <- ss %>% pivot_longer(cols = c(cope, z), names_to = "image_type", values_to = "nii_file") # varcope, [leave out for now]
    afni_out <- sub("(\\+tlrc)*$", "+tlrc", ss$afni_out[1]) # force +tlrc extension
    afni_dir <- dirname(afni_out)
    if (!dir.exists(afni_dir)) {
      dir.create(afni_dir, recursive = TRUE)
      if (!is.null(template_brain)) file.symlink(template_brain, file.path(afni_dir, paste0("template_brain", template_ext)))
    }
    tcatcall <- paste("3dTcat -overwrite -prefix", afni_out, paste(ss$nii_file, collapse = " "))
    run_afni_command(tcatcall)

    z_briks <- which(ss$image_type %in% c("z", "zthresh")) - 1 # subtract 1 because AFNI uses 0-based indexing

    # tack on statistic type as suffix to sub-brik name (avoid ambiguity)
    ss$afni_briks <- paste(ss$afni_briks, ss$image_type, sep="_")

    # label images
    refitcall <- paste0(
      "3drefit -fbuc ", paste("-substatpar", z_briks, "fizt", collapse = " "),
      " -relabel_all_str '", paste(ss$afni_briks, collapse = " "), "' ", afni_out
    )
    run_afni_command(refitcall)

  })
  return(meta_df)

}

#' internal function to convert a feat level 1 directory to an AFNI brik/head file
#' @param feat_dir a .feat directory
#' @param out_filename the directory and filename prefix for the output brik/head file
#' @param what which elements of the feat structure to add to the brik.
#' @param label_prefix an optional character string to add as a prefix to the labels
#' @return a list containing the feat_info for the directory, the full path of the output image,
#'   a vector of the brik names, and a vector of the contrast names
#' @keywords internal
#' @importFrom checkmate assert_directory_exists assert_subset
feat_stats_to_brik <- function(feat_dir, out_filename = "feat_stats",
                               what = c("cope_files", "z_files", "varcope_files"), label_prefix = NULL) {
  if (!checkmate::test_directory_exists(feat_dir)) return(NULL)
  if (!checkmate::test_directory_exists(file.path(feat_dir, "stats"))) return(NULL)
  feat_dir <- normalizePath(feat_dir) # convert to absolute path

  feat_info <- read_feat_dir(feat_dir)
  checkmate::assert_subset(what, names(feat_info))
  checkmate::assert_string(label_prefix, null.ok = TRUE)

  # always insist that outputs are placed inside corresponding feat directory
  afni_out <- file.path(feat_dir, paste0(out_filename, "+tlrc"))

  # interleave images
  to_stitch <- feat_info[what]
  lens <- sapply(to_stitch, length)
  if (length(unique(lens)) != 1L) {
    print(to_stitch)
    stop("unable to sort out how to interleave feat statistics")
  }

  # fancy shortcut for interleaving: https://coolbutuseless.github.io/2019/01/29/interleaving-vectors-and-matrices-part-1/
  stat_order <- do.call(rbind, to_stitch)
  attributes(stat_order) <- NULL
  names(stat_order) <- sub(".*/(.*)\\.nii(\\.gz)*$", "\\1", stat_order, perl = TRUE)
  z_briks <- which(grepl("zstat", names(stat_order))) - 1 # subtract 1 because AFNI uses 0-based indexing

  # get contrast names
  contrast_names <- feat_info$contrast_names
  stat_numbers <- as.numeric(sub(".*(\\d+)$", "\\1", names(stat_order), perl = TRUE))
  stat_names <- names(stat_order)
  brik_names <- sapply(seq_along(stat_numbers), function(x) {
    # these four correspond to contrasts and should get labeled
    if (grepl("(zstat|tstat|cope|varcope)", stat_names[x])) {
      return(paste(contrast_names[stat_numbers[x]], stat_names[x], sep = "_"))
    } else {
      return(stat_names[x]) # don't label with contrast
    }
  })

  # concatenate stat images
  tcatcall <- paste("3dTcat -overwrite -prefix", afni_out, paste(stat_order, collapse = " "))
  run_afni_command(tcatcall)

  if (!is.null(label_prefix)) { # allow user to add label prefix to brik names
    brik_names <- paste0(label_prefix, brik_names)
  }

  # label images
  refitcall <- paste0(
    "3drefit -fbuc ", paste("-substatpar", z_briks, "fizt", collapse = " "),
    " -relabel_all_str '", paste(brik_names, collapse = " "), "' ", afni_out
  )
  run_afni_command(refitcall)

  return(list(feat_dir = feat_dir, afni_img = afni_out, feat_info = feat_info, brik_names = brik_names, contrast_names = contrast_names))
}
