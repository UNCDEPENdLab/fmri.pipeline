#https://bids-specification.readthedocs.io/en/stable/appendices/entities.html

#' Function for extracting fields from BIDS filenames (mostly generated correctly by ChatGPT)
#' @param filenames a character vector of BIDS filenames
#' @return a `data.frame` containing the fields extraced from `filenames`
#' @importFrom checkmate assert_character
#' @keywords internal
#' @examples 
#'  filenames <- c(
#'     "/proj/fmap-phase/task-memory_sub-01_ses-02_run-1_space-MNI2009c_acq-highres_desc-preproc_bold.nii.gz",
#'     "acq-lowres_desc-smoothed_sub-02_task-attention_run-2_bold.nii.gz",
#'    "sub-03_space-MNI152NLin6Asym_task-motor_desc-raw_echo-2_dir-PA_bold.nii.gz",
#'    "hemi-L_desc-denoised_task-vision_rec-magnitude_fmap-phase_sub-04_bold.nii.gz"
#'  )
#' 
#' extract_bids_info(filenames)
extract_bids_info <- function(filenames, drop_unused=FALSE) {
  checkmate::assert_character(filenames)
  filenames <- basename(filenames) # avoid matching on path components

  # Define regex patterns for each BIDS entity
  patterns <- list(
    subject = "sub-(\\d+)",
    session = "ses-(\\d+)",
    task = "task-([a-zA-Z0-9]+)",
    run = "run-(\\d+)",
    space = "space-([a-zA-Z0-9]+)",
    acquisition = "acq-([a-zA-Z0-9]+)",
    description = "desc-([a-zA-Z0-9]+)",
    modality = "mod-([a-zA-Z0-9]+)",
    echo = "echo-(\\d+)",
    direction = "dir-([a-zA-Z0-9]+)",
    hemisphere = "hemi-([a-zA-Z0-9]+)",
    reconstruction = "rec-([a-zA-Z0-9]+)",
    fieldmap = "fmap-([a-zA-Z0-9]+)"
  )
  
  # Function to extract an entity from a filename
  extract_entity <- function(filename, pattern) {
    match <- regmatches(filename, regexpr(pattern, filename))
    if (length(match) > 0) {
      return(sub(".*-", "", match))  # Extract value after the last "-"
    } else {
      return(NA)
    }
  }
  
  # Process each filename
  extracted_info <- lapply(filenames, function(filename) {
    # Extract each entity independently
    info <- lapply(patterns, extract_entity, filename = filename)
    return(as.data.frame(info, stringsAsFactors = FALSE))
  })
  
  # Combine results into a single data frame
  df <- do.call(rbind, extracted_info)
  if (isTRUE(drop_unused)) {
    all_na <- sapply(df, function(i) all(is.na(i)))
    df <- df[!all_na]
  }
  
  return(df)
}

#' Function to generate a run_data object from a BIDS-compliant folder
#' @param bids_dir a directory containing BIDS-compliant processed data for analysis
#' @param modality the subfolder within \code{bids_dir} that contains data of a certain modality.
#'   Almost always 'func', which is the default.
#' @param task_name the name of the task, which is appended with "task-"
#' @param suffix an optional suffix in the expected filename (just before the file extension)
#' @return a data.frame containing all run_nifti and confound_input_file results for subjects in the folder
#' @details The files should generally have a name like
#'   sub-220256_task-ridl3_space-MNI152NLin2009cAsym_desc-preproc_bold_postprocessed.nii.gz
#'   and be located in a folder like: /proj/mnhallqlab/proc_data/sub-220256/func/
#'   where 'func' is the \code{modality}, 'ridl' is the \code{task_name}, and
#'   '_postprocessed' is the \code{suffix}.
#'   If the expected \code{modality} folder is missing, the function will also
#'   search session-level folders (e.g., \code{sub-<id>/ses-<id>/}) or the subject
#'   root for NIfTI files, and fall back to any \code{.nii/.nii.gz} files found
#'   directly in those directories.
#' @importFrom dplyr bind_rows
#' @export
#' @examples 
#' \dontrun{
#'   run_df <- generate_run_data_from_bids("/proj/mnhallqlab/no_backup/flanker-fmriprep", task_name = "flanker", desc = "postproc")
#' }
generate_run_data_from_bids <- function(bids_dir, modality="func", task_name="ridl", desc="postproc", suffix="bold", space=NULL, anat_root=NULL, fmap_root=NULL) {
  checkmate::assert_string(bids_dir)
  checkmate::assert_directory_exists(bids_dir)
  checkmate::assert_string(modality)
  checkmate::assert_string(task_name)
  checkmate::assert_string(desc, na.ok=TRUE, null.ok=TRUE)
  checkmate::assert_string(suffix, na.ok = TRUE, null.ok = TRUE)
  checkmate::assert_string(space, na.ok=TRUE, null.ok=TRUE)
  sub_dirs <- grep("^.*/?sub-", list.dirs(bids_dir, recursive = FALSE), value = TRUE)

  find_search_dirs <- function(ss, modality) {
    ses_dirs <- list.dirs(ss, recursive = FALSE, full.names = TRUE)
    ses_dirs <- ses_dirs[grepl("/ses-", ses_dirs)]
    if (length(ses_dirs) > 0L) {
      dirs <- vapply(ses_dirs, function(sd) {
        mod_dir <- file.path(sd, modality)
        if (dir.exists(mod_dir)) mod_dir else sd
      }, character(1))
      return(unique(dirs))
    }
    mod_dir <- file.path(ss, modality)
    if (dir.exists(mod_dir)) mod_dir else ss
  }

  slist <- lapply(sub_dirs, function(ss) {
    id <- sub("^sub-", "", basename(ss))
    mr_dir <- ss

    if (is.null(anat_root)) anat_dir <- file.path(ss, "anat")
    else anat_dir <- file.path(anat_root, basename(ss), "anat")

    if (dir.exists(anat_dir)) {
      t1w <- Sys.glob(glue("{anat_dir}/sub-{id}*_T1w.nii.gz"))
      if (length(t1w) > 1L) {
        warning(glue("Using first of multiple T1w files: {paste(t1w, collapse=', ')}"))
        t1w <- t1w[1L]
      } else if (length(t1w) == 0L) {
        t1w <- NA_character_
      }
    } else {
      t1w <- NA_character_
    }

    # this is weakly developed -- works for A>>P and P>>A spin echos only
    if (is.null(fmap_root)) fmap_dir <- file.path(ss, "fmap")
    else fmap_dir <- file.path(fmap_root, basename(ss), "fmap")

    if (dir.exists(fmap_dir)) {
      se_pos <- Sys.glob(glue("{fmap_dir}/sub-{id}*-{task_name}_dir-PA*.nii.gz"))
       if (length(se_pos) > 1L) {
        warning(glue("Using first of multiple SE P>>A files: {paste(se_pos, collapse=', ')}"))
        se_pos <- se_pos[1L]
      } else if (length(se_pos) == 0L) {
        se_pos <- NA_character_
      }

      se_neg <- Sys.glob(glue("{fmap_dir}/sub-{id}*-{task_name}_dir-AP*.nii.gz"))
      if (length(se_neg) > 1L) {
        warning(glue("Using first of multiple SE A>>P files: {paste(se_neg, collapse=', ')}"))
        se_neg <- se_neg[1L]
      } else if (length(se_neg) == 0L) {
        se_neg <- NA_character_
      }
    } else {
      se_pos <- NA_character_
      se_neg <- NA_character_
    }

    if (!is.null(desc) && !is.na(desc)) desc <- paste0("_desc-", desc) else desc <- ""
    if (!is.null(suffix) && !is.na(suffix)) {
      if (!substr(suffix, 1, 1) == "_") suffix <- paste0("_", suffix)
    } else {
      suffix <- ""
    }

    if (!is.null(space) && !is.na(space)) space <- paste0("_space-", space) else space <- ""
      

    search_dirs <- find_search_dirs(ss, modality)
    if (!checkmate::test_directory_exists(search_dirs[1L])) {
      warning(glue("Cannot find expected modality/session directory under: {ss}"))
      return(NULL)
    }

    all_rows <- list()
    for (expect_dir in search_dirs) {
      pattern <- glue("^sub-\\d+(_ses-\\d+)?_task-{task_name}(_run-\\d+)?.*{space}.*{desc}.*{suffix}\\.nii(\\.gz)?$")
      # Recursively search for matching NIfTI files in the BIDS directory
      nii_files <- list.files(path = expect_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
      bids_found <- length(nii_files) > 0L

      if (!bids_found) {
        # fallback for non-BIDS layouts: look for .nii/.nii.gz in the directory itself
        nii_files <- list.files(path = expect_dir, pattern = "\\.nii(\\.gz)?$", recursive = FALSE, full.names = TRUE)
      }

      if (length(nii_files) == 0L) {
        warning(glue("No NIfTI file matches in: {expect_dir}"))
        next
      }

    # sbref_files <- Sys.glob(glue("{expect_dir}/sub*_task-{task_name}*sbref.nii.gz"))
    # sbref_files <- list.files(path = expect_dir, pattern = sub(pattern, recursive = TRUE, full.names = TRUE)
    #  (glue("{expect_dir}/sub*_task-{task_name}*sbref.nii.gz"))

      if (bids_found) {
        sbref_files <- sub(suffix, "_sbref", nii_files)
        sbref_files <- sapply(sbref_files, function(x) if (file.exists(x)) x else NA_character_, USE.NAMES=FALSE)
      } else {
        sbref_files <- rep(NA_character_, length(nii_files))
      }
    
    # if (length(sbref_files) > 0L && length(nii_files) != length(sbref_files)) {
    #   warning(glue("Cannot align nifti and sbref files for {expect_dir}"))
    #   return(NULL)
    # } else if (length(sbref_files) == 0L) {
    #   sbref_files <- NA_character_
    # }

      if (bids_found) {
        confound_files <- list.files(path = expect_dir, pattern = glue("sub.*_task-{task_name}.*confounds.*\\.tsv$"), recursive = TRUE, full.names = TRUE)
        if (length(confound_files) > 0L && length(nii_files) != length(confound_files)) {
          warning(glue("Cannot align nifti and confound files for {expect_dir}"))
          return(NULL)
        } else if (length(confound_files) == 0L) {
          confound_files <- rep(NA_character_, length(nii_files))
        }

        events_files <- list.files(path = expect_dir, pattern = glue("sub.*_task-{task_name}.*_events.*\\.tsv$"), recursive = TRUE, full.names = TRUE)
        if (length(events_files) > 0L && length(nii_files) != length(events_files)) {
          warning(glue("Cannot align nifti and events files for {expect_dir}"))
          return(NULL)
        } else if (length(events_files) == 0L) {
          events_files <- rep(NA_character_, length(nii_files))
        }
      } else {
        confound_files <- rep(NA_character_, length(nii_files))
        events_files <- rep(NA_character_, length(nii_files))
      }

    # if only one matching NIfTI is found, we likely have a single-run scenario and may not expect a run-<x> syntax
      if (all(grepl(glue(".*sub-.*_task-{task_name}_run-\\d+.*"), nii_files, perl=TRUE))) {
        run_number <- as.integer(sub(glue(".*sub-.*_task-{task_name}_run-(\\d+).*"), "\\1", nii_files, perl = TRUE))
      } else if (length(nii_files) == 1L) {
        run_number <- 1
      } else {
        run_number <- seq_along(nii_files)
      }

      all_rows[[length(all_rows) + 1L]] <- data.frame(
        id, task_name, mr_dir, run_number,
        run_nifti = nii_files,
        sbref_nifti = sbref_files,
        confound_input_file = confound_files,
        events_file = events_files,
        t1w, se_pos, se_neg,
        stringsAsFactors = FALSE
      )
    }

    if (length(all_rows) == 0L) return(NULL)
    dplyr::bind_rows(all_rows)
  })

  dplyr::bind_rows(slist)
}


#' Function to generate a trial_data object from a BIDS-compliant folder
#' @param bids_dir a directory containing BIDS-compliant processed data for analysis
#' @param modality the subfolder within \code{bids_dir} that contains data of a certain modality.
#'   Almost always 'func', which is the default.
#' @param type at present, always 'task' to denote this part of the BIDS filename... Not totally sure what else it could be
#' @param task_name the name of the task, which is appended with \code{type}
#' @param suffix an optional suffix in the expected filename (just before the file extension)
#' @return a data.frame containing all run_nifti and confound_input_file results for subjects in the folder
#' @details The files should generally have a name like
#'   sub-220256_task-ridl3_space-MNI152NLin2009cAsym_desc-preproc_bold_postprocessed.nii.gz
#'   and be located in a folder like: /proj/mnhallqlab/proc_data/sub-220256/func/
#'   where 'func' is the \code{modality}, 'task' is the \code{type}, 'ridl' is the \code{task_name}, and
#'   '_postprocessed' is the \code{suffix}.
#' @importFrom dplyr bind_rows
#' @importFrom tidyr unnest
#' @export
#' @examples
#' \dontrun{
#'   df <- generate_trial_data_from_bids("/proj/mnhallqlab/no_backup/flanker-fmriprep", task_name = "flanker")
#' }
generate_trial_data_from_bids <- function(bids_dir, modality="func", task_name="ridl") {
  checkmate::assert_directory_exists(bids_dir)
  checkmate::assert_string(modality)
  checkmate::assert_string(task_name)
  sub_dirs <- grep("^.*/?sub-", list.dirs(bids_dir, recursive = FALSE), value = TRUE)

  slist <- lapply(sub_dirs, function(ss) {
    id <- sub("^sub-", "", basename(ss))
    mr_dir <- ss

    expect_dir <- file.path(ss, modality)
    if (!checkmate::test_directory_exists(expect_dir)) {
      warning(glue("Cannot find expected modality directory: {expect_dir}"))
      return(NULL)
    }

    events_files <- Sys.glob(glue("{expect_dir}/sub*_task-{task_name}*_events*.tsv"))

    if (length(events_files) == 0L) {
      events_files <- NA_character_
    } else {
      events_df <- extract_bids_info(events_files, drop_unused = TRUE)
    }
    
    events_df$events_file <- events_files
    events_df$data <- lapply(events_df$events_file, read.delim, header = TRUE, sep = "\t")
    
    return(events_df)

  })

  dplyr::bind_rows(slist) %>% tidyr::unnest(data) %>% dplyr::select(-events_file)
}

#' Function to generate a subject_data object from a BIDS-compliant folder
#' @param bids_dir a directory containing BIDS-compliant processed data for analysis
#' @details This function basically just reads participants.tsv from the root of the BIDS folder and renames
#'   `participant_id`` to `id` to match this pipeline's conventions
#' @return a data.frame containing subject-level data contained in participants.tsv
#' @importFrom checkmate assert_string assert_directory_exists
#' @export
generate_subject_data_from_bids <- function(bids_dir) {
  checkmate::assert_string(bids_dir)
  checkmate::assert_directory_exists(bids_dir)

  if (!file.exists(file.path(bids_dir, "participants.tsv"))) {
    warning("No participants.tsv file found in: ", bids_dir)
    return(NULL)
  }

  subject_df <- read.table(file.path(bids_dir, "participants.tsv"),
    header = TRUE, stringsAsFactors = FALSE,
    colClasses = c(participant_id = "character")
  )
  
  # generate expected location of subject files
  subject_df$mr_dir <- file.path(bids_dir, subject_df$participant_id)

  # strip sub- prefix
  subject_df$participant_id <- sub("sub-", "", subject_df$participant_id)

  names(subject_df)[which(names(subject_df) == "participant_id")] <- "id" # to match our pipeline's conventions
  return(subject_df)
}
