#' wrapper class for AFNI whereami
#' @export
afni_whereami <- R6::R6Class("afni_whereami",
  private = list(
    pvt_atlases = c("MNI_Glasser_HCP_v1.0", "Brainnetome_1.0", "CA_ML_18_MNI"),
    pvt_whereami_call = NULL,
    parse_whereami_output = function(txt) {
      # This row demarcates a new region in the output. Split on this basis
      section_map <- grep("+++++++ nearby Atlas structures +++++++", lookup, fixed = TRUE)

      # Trim off any lines that precede the first region, then adjust the section boundaries accordingly
      lookup <- lookup[min(section_map):length(lookup)]
      section_map <- section_map - min(section_map) + 1

      # Create an integer marker for each section, used for splitting into a list
      section_split <- rep.int(seq_along(section_map), times = diff(c(section_map, length(lookup) + 1)))
      lookup_split <- split(lookup, section_split)

      roi_lookup <- sapply(seq_along(lookup_split), function(ii) {
        sec <- lookup_split[[ii]]
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
      lookup_df <- bind_rows(roi_lookup)

      coordlines <- grep("Focus point (LPI)=", lookup, fixed = TRUE)
      coords <- lookup[coordlines + 2] # first line after header is TLRC, second is MNI
      # coords <- sub("<a href=.*$", "", coords, perl=TRUE)
      coords <- sub("^\\s*(-?\\d+\\s*mm.*\\{MNI\\})\\s*<a href=.*$", "\\1", coords, perl = TRUE)

      coords_l <- as.numeric(sub("^\\s*(-*\\d+) mm.*", "\\1", coords, perl = TRUE))
      coords_p <- as.numeric(sub("^\\s*(-*\\d+) mm \\[(?:L|R)\\],\\s+(-*\\d+) mm.*", "\\2", coords, perl = TRUE))
      coords_i <- as.numeric(sub("^\\s*(-*\\d+) mm \\[(?:L|R)\\],\\s+(-*\\d+) mm \\[(?:A|P)\\],\\s+(-*\\d+) mm.*", "\\3", coords, perl = TRUE))
    },
    parse_whereami_bset_output = function(txt) {

    },
    build_call = function() {

    }
  ),
  public = list(
    initialize = function(afni_3dclusterize_obj = NULL, coord_file = NULL, coord_file_columns = NULL,
                          coords_vector = NULL, atlases = NULL) {

    },
    get_call = function() {

    },
    run = function() {

      # get coordinates and names of regions
      lookup <- runAFNICommand(paste0("whereami -coord_file ", clust_1d, "'[1,2,3]' -space MNI -lpi -atlas CA_ML_18_MNIA"),
        stderr = "/dev/null", intern = TRUE
      )

      lookup <- readLines("/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/mtest")

      exitstatus <- attr(lookup, "status")
      if (!is.null(exitstatus) && exitstatus != 0) next # whereami failed, which occurs when there are no clusters. Skip to next tbrik

      # get voxel sizes of clusters
      vsizes <- read.table(clust_1d)$V1
    }
  ),
)