#' Populates parcel-wise statistics from an atlas or set of clusters back into NIfTIs in the original space
#' 
#' @param atlas_nifti The filename of the NIfTI image containing atlas values to match against \code{stat_dt}
#' @param stat_dt a \code{data.table} object containing the statistics to write to parcels in the atlas
#' @param stat_cols The column names in \code{stat_dt} that should be written to parcels in the NIfTI. By default,
#'   each column will be written as a separate sub-brik (i.e., along the 4th dimension of the NIfTI), but this can be
#'   modified by the \code{stack_along} argument.
#' @param stat_labels A character vector of the same length as \code{stat_cols} giving the names for each statistic
#'   that will be used in the output NIfTI files. By default the labels are simply the \code{stat_cols}.
#' @param atlas_col The column name in \code{stat_dt} containing the integer values that correspond to parcel numbers
#'   in \code{atlas_nifti}. This is used to map the statistics back into the \code{atlas_nifti} space.
#' @param split_on An optional character vector that yields an output image for each unique combination of splits. For
#'   example, if \code{stat_dt} contains multiple GLM contrasts in \code{"contrast"} and multiple subjects in \code{"subject"},
#'   then \code{split_on = c("contrast", "subject")} would yield a separate NIfTI with the statistics for each contrast and
#'   subject. This allows the fill operation to be repeated over an arbitrarily stacked dataset.
#' @param stack_along By default, this is \code{NULL}, in which case each output NIfTI will have a sub-briks for each statistic
#'   in \code{stat_cols}. Alternatively, \code{stack_along} can be a character vector of other columns in \code{stat_dt}. In this case,
#' @param img_prefix A string prefixing the filenames for all NIfTIs
#' @param out_dir The directory in which to place the stat NIfTIs. If not present, it will be created
#' @param overwrite Whether to allow existing files to be overwritten. Default: FALSE
#' @param afni_dir The location of your AFNI installation, which is used to cobble together NIfTIs
#'
#' @importFrom RNifti readNifti
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom checkmate assert_file_exists assert_directory_exists assert_subset assert_flag
#' @importFrom data.table setDT
#' @importFrom glue glue
#' @export
fill_atlas_with_stats <- function(atlas_nifti, stat_dt, stat_cols = c("t", "p"), stat_labels = NULL,
                                 atlas_col = "atlas_value", split_on = NULL, stack_along = NULL,
                                 img_prefix = "atlasfill", out_dir = getwd(), overwrite = FALSE, afni_dir = "~/abin") {

  #input validation
  checkmate::assert_file_exists(atlas_nifti)

  if (!checkmate::test_data_table(stat_dt)) {
    message("Coercing stat_dt to data.table object")
    data.table::setDT(stat_dt)
  }

  checkmate::assert_subset(stat_cols, names(stat_dt))

  # default labels to match stat column names
  if (is.null(stat_labels)) {
    stat_labels <- stat_cols
  } else {
    checkmate::assert_character(stat_labels, len = length(stat_cols))
  }

  checkmate::assert_string(atlas_col)
  checkmate::assert_subset(atlas_col, names(stat_dt))

  checkmate::assert_character(split_on, null.ok = TRUE)
  checkmate::assert_character(stack_along, null.ok = TRUE)
  checkmate::assert_string(img_prefix)
  checkmate::assert_string(out_dir, null.ok = TRUE)

  if (is.null(out_dir)) out_dir <- getwd()
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  checkmate::assert_flag(overwrite)

  checkmate::assert_directory_exists(afni_dir)
  afni_dir <- normalizePath(afni_dir)

  # handle non-standard stacking
  if (!is.null(stack_along)) {
    # message("Reshaping input dataset to put stat_cols as long and stack_along as wide")
    stat_dt <- stat_dt %>%
      dplyr::select(all_of(stat_cols), all_of(stack_along), all_of(split_on), all_of(atlas_col)) %>%
      tidyr::pivot_longer(cols = all_of(stat_cols), names_to = ".statistic", values_to = ".value") %>%
      #tidyr::pivot_wider(names_glue = "stack_{stack_along}_{.value}", names_from = all_of(stack_along), values_from = ".value") %>%
      tidyr::pivot_wider(names_prefix = "stack_", names_from = all_of(stack_along), values_from=".value") %>%
      data.table::setDT()

    split_on <- c(split_on, ".statistic")
    stack_cols <- grep("^stack_", names(stat_dt), value = TRUE)
    stat_labels <- sub("stack_", "", stack_cols)
  } else {
    stack_cols <- stat_cols
  }

  # handle splits
  if (!is.null(split_on)) {
    stat_dt <- split(stat_dt, by = split_on)
    if (is.null(img_prefix)) {
      img_names <- paste0(make.names(names(stat_dt)), ".nii.gz")
    } else {
      img_names <- paste0(img_prefix, "_", make.names(names(stat_dt)), ".nii.gz")
    }
  } else {
    stat_dt <- list(stat_dt) # single element
    img_names <- paste0(img_prefix, ".nii.gz")
  }

  use_afni <- TRUE # locked to TRUE based on prior testing -- some vestigial code for FALSE preserved
  if (isTRUE(use_afni)) {
    # force atlas to short data type to ensure that 3dUndump is okay with it
    orig <- atlas_nifti # user specified atlas
    atlas_nifti <- tempfile(pattern="atlas", fileext=".nii.gz") # where the dtype = short file goes
    system2(
      command = glue("{afni_dir}/nifti_tool"), stdout = FALSE, stderr = FALSE,
      args = glue("-copy_image -infiles {orig} -prefix {atlas_nifti} -convert2dtype DT_SIGNED_SHORT -convert_fail_choice fail -convert_verify")
    )
  } else {
    vol <- RNifti::readNifti(atlas_nifti)
  }

  # helper subfunction for populating nifti  
  fill_nifti <- function(data, out_file) {
    to_combine <- sapply(seq_along(stack_cols), function(x) {
      tempfile(fileext = ".nii.gz")
    })

    # loop over subbriks to create (e.g., t, p, and logp)
    for (ff in seq_along(stack_cols)) {
      if (file.exists(to_combine[ff]) && isFALSE(overwrite)) {
        message(glue("Skipping existing file {to_combine[ff]}"))
        next
      }

      if (use_afni) {
        # use 3dUndump, which is about 50% faster than internal code below
        stat_val <- data[[stack_cols[ff]]]
        stat_tmp <- tempfile(pattern="stat", fileext=".txt")
        writeLines(as.character(stat_val), con = stat_tmp)
        system2(
          command = glue("{afni_dir}/3dUndump"),
          args = glue("-overwrite -ROImask {atlas_nifti} -datum float -prefix {to_combine[ff]} {stat_tmp}"),
          stdout = FALSE, stderr = FALSE
        )
      } else {
        new_stat <- vol # copy
        new_stat <- new_stat * 0 # start with all zeros

        # prior benchmarking indicated that the loop-based assignment was no slower than vectorization over atlas values
        for (val in unique(data[[atlas_col]])) {
          stat_val <- data[[stack_cols[ff]]][data[[atlas_col]] == val]
          if (length(stat_val) > 1L) {
            stop("more than one stat value found for a parcel. Problem with splits?")
          }
          new_stat[vol == val] <- stat_val # the RHS should be one value... probably validate this if I extend the function
        }

        writeNifti(new_stat, file = to_combine[ff])
      }
    }

    # AFNI does not tolerate spaces in output file, so need to use a tmpfile and then move it
    tmpout <- tempfile(fileext = ".nii.gz")

    system2(
      command = glue("{afni_dir}/3dTcat"), stdout = FALSE, stderr = FALSE,
      args = glue("-overwrite -prefix '{tmpout}' {paste(to_combine, collapse=' ')}")
    )
    system2(
      command = glue("{afni_dir}/3drefit"), stdout = FALSE, stderr = FALSE,
      args = glue("-relabel_all_str '{paste(stat_labels, collapse=' ')}' '{tmpout}'")
    )

    file.copy(tmpout, out_file, overwrite = TRUE)
    unlink(c(tmpout, to_combine))
  }

  
  # loop over splits (e.g., contrasts)
  for (dd in seq_along(stat_dt)) {
    this_out_nifti <- glue("{out_dir}/{img_names[dd]}")
    if (checkmate::test_file_exists(this_out_nifti) && isFALSE(overwrite)) {
      message(glue("Skipping existing file: {this_out_nifti}"))
      next
    }

    data <- stat_dt[[dd]]

    fill_nifti(data, this_out_nifti)
  }
}
