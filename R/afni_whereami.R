#' wrapper class for AFNI whereami
#' @importFrom tidyselect everything
#' @importFrom dplyr select
#' @export
afni_whereami <- R6::R6Class("afni_whereami",
  private = list(
    pvt_method = "coord_file", # coord_file or coord_vector
    pvt_coord_file = NULL,
    pvt_coord_file_columns = NULL,
    pvt_coord_vector = NULL,
    pvt_atlases = c("MNI_Glasser_HCP_v1.0", "Brainnetome_1.0", "CA_ML_18_MNI"),
    pvt_omask = NULL,
    pvt_orient = "LPI",
    pvt_space = "MNI",
    pvt_call = NULL,
    pvt_call_omask = NULL,
    pvt_output_file = NULL,
    pvt_omask_output_file = "whereami_omask.txt",
    pvt_afnidir = NULL,
    pvt_whereami_df = NULL,

    parse_whereami_output = function(txt) {
      if (!is.null(private$pvt_whereami_df) && nrow(private$pvt_whereami_df) > 0L) {
        # data have already been parsed... may need a mechanism for forced reparsing?
        return(private$pvt_whereami_df)
      }

      # This row demarcates a new region in the output. Split on this basis
      section_map <- grep("+++++++ nearby Atlas structures +++++++", txt, fixed = TRUE)

      # Trim off any lines that precede the first region, then adjust the section boundaries accordingly
      txt <- txt[min(section_map):length(txt)]
      section_map <- section_map - min(section_map) + 1

      # Create an integer marker for each section, used for splitting into a list
      section_split <- rep.int(seq_along(section_map), times = diff(c(section_map, length(txt) + 1)))
      txt_split <- split(txt, section_split)

      roi_list <- lapply(seq_along(txt_split), function(ii) {
        sec <- txt_split[[ii]]
        atlas_lines <- grep("^Atlas .*", sec, perl = TRUE)
        nomatch <- grep("***** Not near any region stored in databases *****", sec, fixed = TRUE)
        if (length(nomatch) > 0L) {
          # return("Unable to identify label")
          return(list(roi_num = ii)) # need at least one value for bind_rows to work as expected
        } else {
          atlas_names <- sub("^Atlas\\s+([^:]+).*", "\\1", sec[atlas_lines], perl = TRUE)
          best_guess <- lapply(atlas_lines, function(aa) {
            return(trimws(sec[aa + 1])) # first match after atlas for each cluster
          }) %>% setNames(atlas_names)
          best_guess[["roi_num"]] <- ii
          return(best_guess)
        }
      })

      # create a data.frame where each ROI has one row and the columns are the best labels from each atlas specified
      roi_df <- bind_rows(roi_list)

      coordlines <- grep("Focus point (LPI)=", txt, fixed = TRUE)
      if (length(coordlines) == 0L) {
        stop("Parser cannot find focus point lines in output. Cannot continue.")
      }

      coords <- txt[coordlines + 2] # first line after header is TLRC, second is MNI
      # coords <- sub("<a href=.*$", "", coords, perl=TRUE)
      coords <- sub("^\\s*(-?\\d+\\s*mm.*\\{MNI\\})\\s*<a href=.*$", "\\1", coords, perl = TRUE)

      coords_l <- as.numeric(sub("^\\s*(-*\\d+) mm.*", "\\1", coords, perl = TRUE))
      coords_p <- as.numeric(sub("^\\s*(-*\\d+) mm \\[(?:L|R)\\],\\s+(-*\\d+) mm.*", "\\2", coords, perl = TRUE))
      coords_i <- as.numeric(sub("^\\s*(-*\\d+) mm \\[(?:L|R)\\],\\s+(-*\\d+) mm \\[(?:A|P)\\],\\s+(-*\\d+) mm.*", "\\3", coords, perl = TRUE))

      stopifnot(length(coords) == nrow(roi_df)) # cannot continue if these don't match.

      roi_df <- roi_df %>%
        dplyr::mutate(x = coords_l, y = coords_p, z = coords_i) %>%
        dplyr::select(roi_num, x, y, z, everything())

      return(roi_df)
    },
    parse_whereami_bmask_output = function(txt) {

    },
    build_call = function() {
      orient_str <- ifelse(private$pvt_orient == "LPI", "-lpi", "-rai")
      atlas_str <- paste("-atlas", private$pvt_atlases, collapse=" ")

      if (private$pvt_method == "coord_file") {
        input_str <- glue("-coord_file {private$pvt_coord_file}'[{paste(private$pvt_coord_file_columns, collapse=',')}]'")
      } else if (private$pvt_method == "coord_vector") {
        input_str <- paste(private$pvt_coord_vector, collapse=" ")
      }

      str <- glue("whereami {input_str} {atlas_str} {orient_str} -space {private$pvt_space} > {private$pvt_output_file}")

      private$pvt_call <- str

      # also build omask call
      if (!is.null(private$pvt_omask)) {
        omask_str <- glue("whereami -omask {private$pvt_omask} {atlas_str} {orient_str} -space {private$pvt_space} > {private$pvt_omask_output_file}")
        private$pvt_call_omask <- omask_str
      }

    }
  ),
  public = list(
    initialize = function(afni_3dclusterize_obj = NULL, omask = NULL, coord_file = NULL, coord_file_columns = NULL,
                          coord_vector = NULL, atlases = NULL, coord_orientation = NULL, coord_space = NULL, 
                          output_file = NULL, omask_output_file = NULL, afnidir = NULL) {

      if (!is.null(afni_3dclusterize_obj)) {
        checkmate::assert_class(afni_3dclusterize_obj, "afni_3dclusterize")
        coord_orientation <- afni_3dclusterize_obj$get_orient()
        ofiles <- afni_3dclusterize_obj$get_output_files()
        omask <- ofiles["cluster_map"] # file containing integer-valued clusters
        coord_file <- ofiles["cluster_table"]
        coord_file_columns <- 1:3 # 0-based indexing in AFNI 1D parser, so this should be CM LR, CM PA, and CM IS
      }

      if (!is.null(coord_file)) {
        checkmate::assert_file_exists(coord_file)
        private$pvt_coord_file <- coord_file
        private$pvt_method <- "coord_file"

        if (is.null(coord_file_columns)) {
          stop("If a coord_file is provided, you must also provide the coord_file_columns corresponding to X, Y, Z coordinates.")
        } else {
          checkmate::assert_integerish(coord_file_columns, lower = 0, upper = 1e5, any.missing = FALSE, len = 3L)
          private$pvt_coord_file_columns <- as.integer(coord_file_columns)
        }
      }

      # use a triplet of X, Y, Z coordinates as input (rather than a -coord_file)
      if (!is.null(coord_vector)) {
        checkmate::assert_numeric(coord_vector, lower = -1e3, upper = 1e3, any.missing = FALSE, len = 3L)
        private$pvt_method <- "coord_vector"
        private$pvt_coord_vector <- coord_vector
      }

      if (!is.null(coord_orientation)) {
        checkmate::assert_string(coord_orientation)
        coord_orientation <- toupper(coord_orientation)
        checkmate::assert_subset(coord_orientation, c("LPI", "RAI"))
        private$pvt_orient <- coord_orientation
      }

      if (!is.null(atlases)) {
        checkmate::assert_character(atlases, all.missing = FALSE)
        # should probably validate them... but for now, let AFNI sort it
        private$pvt_atlases <- atlases
      }

      if (!is.null(omask) && !is.na(omask)) {
        checkmate::assert_file_exists(omask)
        private$pvt_omask <- omask

        # location of whereami output for omask
        if (!is.null(omask_output_file)) {
          checkmate::assert_string(omask_output_file)
          private$pvt_omask_output_file <- omask_output_file
        }

        # convert output file to absolute path, using location of omask if user doesn't provide absolute path.
        private$pvt_omask_output_file <- R.utils::getAbsolutePath(private$pvt_omask_output_file, workDirectory = dirname(private$pvt_omask))
      }

      # convert output file to absolute path, using location of coord_file if user doesn't provide absolute path.
      if (private$pvt_method == "coord_file") {
        wd <- dirname(private$pvt_coord_file)
      } else {
        wd <- getwd()
      }

      if (!is.null(output_file)) {
        checkmate::assert_string(output_file)
        private$pvt_output_file <- output_file
      }

      # default to naming clusters file according to the inset file
      if (is.null(private$pvt_output_file)) {
        if (private$pvt_method == "coord_file") {
          private$pvt_output_file <- paste0(basename(file_sans_ext(private$pvt_coord_file)), "_whereami.txt")
        } else {
          private$pvt_output_file <- "whereami.txt" # a bit dangerous if you have many stats in the same folder (overwriting)
        }
      }

      private$pvt_output_file <- R.utils::getAbsolutePath(private$pvt_output_file, workDirectory = wd)

      if (!is.null(coord_space)) {
        checkmate::assert_string(coord_space)
        coord_space <- toupper(coord_space)
        checkmate::assert_subset(coord_space, c("MNI", "MNI_ANAT", "TLRC"))
        private$pvt_space <- coord_space
      }

      if (!is.null(afnidir)) {
        checkmate::assert_directory_exists(afnidir)
        private$pvt_afnidir <- afnidir
      }

    },
    get_call = function() {
      private$build_call()
      private$pvt_call
    },
    get_omask_call = function() {
      private$build_call()
      private$pvt_call_omask
    },
    get_whereami_df = function() {
      if (!checkmate::test_file_exists(private$pvt_output_file)) {
        warning("Expected whereami output directory does not exist. Cannot return data! ", private$pvt_output_file)
      }

      txt <- readLines(private$pvt_output_file)
      return(private$parse_whereami_output(txt))
    },
    get_omask_df = function() {

    },
    get_output_files = function(exclude_missing = TRUE) {
      whereami_from_coords <- private$pvt_output_file
      whereami_omask_overlap <- private$pvt_omask_output_file

      if (isTRUE(exclude_missing)) {
        if (!checkmate::test_file_exists(whereami_from_coords)) whereami_from_coords <- NA_character_
        if (!checkmate::test_file_exists(whereami_omask_overlap)) whereami_omask_overlap <- NA_character_
      }

      named_vector(whereami_from_coords, whereami_omask_overlap)
    },
    run = function(force = FALSE) {
      private$build_call()
      outfile_exists <- checkmate::test_file_exists(private$pvt_output_file)
      # NULL means that omask is irrelevant, FALSE means we need to run it, TRUE means we already have the omask result
      omaskfile_exists <- if (is.null(private$pvt_omask)) NULL else checkmate::test_file_exists(private$pvt_omask_output_file)

      if (isFALSE(force) && isTRUE(outfile_exists) && isTRUE(omaskfile_exists)) {
        message("whereami output file already exists: ", private$pvt_output_file, ". Use $run(force=TRUE) if you want to regenerate this file.")
        return(invisible(NULL))
      }

      result <- 0

      if (isFALSE(outfile_exists)) {
        result <- run_afni_command(private$pvt_call, afnidir = private$pvt_afnidir, stderr = "/dev/null")

        if (result != 0) {
          warning("whereami returned an exit status of: ", result)
        }
      }

      if (isFALSE(omaskfile_exists)) {
        result <- run_afni_command(private$pvt_call_omask, afnidir = private$pvt_afnidir, stderr = "/dev/null")

        if (result != 0) {
          warning("whereami returned an exit status of: ", result)
        }
      }

      return(invisible(result))
    }
  ),
)
