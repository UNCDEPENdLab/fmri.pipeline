#' R6 class for pTFCE specification for one or more z-statistic images
#' 
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_directory_exists assert_subset assert_string
#' @export
ptfce_spec <- R6::R6Class("ptfce_spec",
  private = list(
    gfeat_info = NULL,
    z_files = NULL,
    mask_files = NULL,
    residuals_file = NULL,
    dof = NULL,
    fsl_smoothest_file = NULL,
    write_thresh_imgs = TRUE,
    smoothness_method = NULL,
    pvt_two_sided = TRUE,
    pvt_fwe_p = NULL,
    expect_list = NULL,
    scheduler = "slurm",
    time_per_zstat = "12:00", # 12 minutes for one run of pTFCE (at 2.3mm, takes about 2.5 mins for me)
    memgb_per_command = 8, # 8 GB requested for every ptfce run
    set_fwep = function(vec) {
      checkmate::assert_numeric(vec, lower = 1e-10, upper = .9999, any.missing = FALSE, unique = TRUE, null.ok = FALSE)
      private$pvt_fwe_p <- vec
    },
    set_twosided = function(val) {
      checkmate::assert_logical(val, len = 1L, null.ok = FALSE)
      private$pvt_two_sided <- val
    },
    ptfce_outfile_from_zstat = function(z_file, what="ptfce", fwe_p = NULL) {
      z_dir <- normalizePath(dirname(z_file))
      ext <- sub(".*(\\.nii(\\.gz)?)$", "\\1", z_file, perl = TRUE)
      base <- sub(ext, "", basename(z_file), fixed = TRUE)
      if (what == "ptfce") {
        if (is.null(fwe_p)) {
          return(file.path(z_dir, paste0(base, "_ptfce", ext)))
        } else if (!is.null(fwe_p)) {
          fwe_p <- round(fwe_p, 3) # for consistency in file-naming
          return(file.path(z_dir, paste0(base, "_ptfce_fwep_", fwe_p, ext)))
        }
      } else if (what == "csv") {
        return(file.path(z_dir, paste0(base, "_ptfce_zthresh.csv")))
      }
    },
    set_expected_files = function() {
      f_list <- list()
      for (ii in seq_along(private$z_files)) {
        this_z <- private$z_files[ii]
        f_list[[ii]] <- list(
          #z_input = this_z,
          #mask_input = private$mask_files[ii],
          z_ptfce = private$ptfce_outfile_from_zstat(this_z),
          z_csv = private$ptfce_outfile_from_zstat(this_z, what = "csv")
        )
        for (pp in private$pvt_fwe_p) {
          pp <- round(pp, 3) # for consistency in file-naming
          f_list[[ii]][[paste0("z_ptfce_fwep_", pp)]] <- private$ptfce_outfile_from_zstat(this_z, fwe_p = pp)
        }
      }
      private$expect_list <- f_list
    }
  ),
  active = list(
    #' @field fwe_p a vector of p-values used for familywise error (FWE) z-statistic threshold calculations in pTFCE.
    fwe_p = function(vec) {
      if (missing(vec)) {
        return(private$pvt_fwe_p)
      } else {
        private$set_fwep(vec)
      }
    },
    two_sided = function(val) {
      if (missing(val)) {
        return(private$pvt_two_sided)
      } else {
        private$set_twosided(val)
      }
    }
  ),
  public = list(
    #' @param gfeat_dir a .gfeat folder containing a higher-level FSL analysis. This will be used for zstat images
    #'   and mask files.
    #' @param zstat_numbers if a \code{gfeat_dir} is used, a vector of zstat numbers can also be provided to
    #'   subset the zstat images that are used in pTFCE correction. Ignored if gfeat_dir is not provided.
    #' @param z_files A vector of z-statistic filenames that should be corrected using pTFCE
    #' @param mask_files A vector of mask filenames that correspond to \code{z_files}. If this is of length 1,
    #'   then the mask file will be recycled for all zstat images.
    #' @param fwe_p A vector of p-values for which z-statistic thresholds will be calculated.
    #' @param two_sided If TRUE, p-values for \code{fwe_p} are treated as two-tailed (i.e., p-values are divided by 2 in
    #'   the TFCE z-threshold calculation.)
    #' @param write_thresh_imgs If TRUE, then pTFCE thresholds for each fwe_p will be applied to the TFCE image and saved
    #'   to the same folder as the z-statistic. These thresholded files let you look at the map at a given FWE threshold.
    #' @param scheduler Which scheduler to use for submitting jobs. Options are 'local', 'slurm', and 'torque'.
    #' @param time_per_zstat The amount of time to budget for each zstat to run through pTFCE in dd-hh:mm:ss format.
    #'   Default is 10:00 (10 minutes).
    initialize = function(gfeat_dir = NULL, zstat_numbers = NULL, fsl_smoothest_file = NULL, dof = NULL, residuals_file = NULL,
      z_files = NULL, mask_files = NULL, fwe_p = .05, two_sided = TRUE, write_thresh_imgs = TRUE, 
      scheduler = NULL, time_per_zstat = NULL, memgb_per_command = NULL) {

      if (!is.null(gfeat_dir)) {
        if (!is.null(z_files)) {
          stop("Cannot provide both z_files and gfeat_dir as inputs")
        }

        checkmate::assert_directory_exists(gfeat_dir)
        private$gfeat_info <- read_gfeat_dir(gfeat_dir)
        # TODO: need to support subsetting by cope number (e.g., cope3.feat) and zstat number (e.g., zstat1.nii.gz)
        # not currently supported (no high-priority use case)
        z_list <- lapply(private$gfeat_info$cope_dirs, function(x) {
          if (!is.null(zstat_numbers)) {
            checkmate::assert_integerish(zstat_numbers)
            x$cope_df %>%
              filter(cope_number %in% !!zstat_numbers) %>%
              pull(z)
          } else {
            x$cope_df %>% pull(z)
          }
        })

        # number of z stats in each .feat folder
        z_lens <- sapply(z_list, length)

        # build out mask files to match each cope directory zstats
        mask_files <- rep(sapply(private$gfeat_info$cope_dirs, "[[", "mask_file", USE.NAMES = FALSE), z_lens)
        fsl_smoothest_file <- rep(sapply(private$gfeat_info$cope_dirs, "[[", "smoothness_file", USE.NAMES = FALSE), z_lens)
        
        # for now, if gfeat_dir is provided, just use smoothness file, which is computed by smoothest on the res4D
        # this is more efficient than having ptfce handle it internally
        #dof <- rep(sapply(private$gfeat_info$cope_dirs, function(x) x$parsed_txt$dof), z_lens)
        
        # To support multiple L2 cope folders in an L3 gfeat folder, we need to allow for multiple stats/smoothness and dof files
        # these can vary by the .feat folder. Thus, build out logic that class has an fsl_smoothest element for every z_file
        # example:
        # [mnhallq@longleaf-login5 v_entropy-Intercept-Age.gfeat]$ head -n 1 ./cope3.feat/stats/smoothness
        # DLH 0.139392
        # [mnhallq@longleaf-login5 v_entropy-Intercept-Age.gfeat]$ head -n 1 ./cope5.feat/stats/smoothness
        # DLH 0.137736

        private$z_files <- unname(do.call(c, z_list))
      } else if (!is.null(z_files)) {        
        checkmate::assert_file_exists(z_files)
        private$z_files <- z_files
      } else {
        stop("Either z_files or gfeat_dir must be provided as an input to ptfce_spec")
      }

      if (is.null(mask_files)) {
        stop("You must supply one or more mask files to be used in pTFCE calculations.")
      } else if (length(mask_files) == 1L && length(private$z_files) > 1L) {
        # replicate mask file for all inputs
        private$mask_files <- rep(mask_files, length(private$z_files))
      } else {
        stopifnot(length(mask_files) == length(private$z_files))
        private$mask_files <- mask_files
      }

      # verify that all mask files exist
      checkmate::assert_file_exists(private$mask_files)

      if (!is.null(residuals_file)) {
        checkmate::assert_file_exists(residuals_file)
        private$residuals_file <- residuals_file

        if (is.null(dof)) {
          stop("If a residuals file is provided, you also need to provide dof (either as an integer or a file containing that integer).")
        } else if (checkmate::test_integerish(dof)) {
          private$dof <- as.integer(dof)
        } else if (checkmate::test_file_exists(dof)) {
          private$dof <- as.integer(readLines(dof))
        } else {
          stop("Can't sort out this dof input")
        }

        private$smoothness_method <- "residuals"
      } else if (!is.null(fsl_smoothest_file)) {
        checkmate::assert_file_exists(fsl_smoothest_file)
        if (length(fsl_smoothest_file) == 1L) {
          # replicate smoothest for all inputs
          private$fsl_smoothest_file <- rep(fsl_smoothest_file, length(private$z_files))
        } else {
          stopifnot(length(fsl_smoothest_file) == length(private$z_files))
          private$fsl_smoothest_file <- fsl_smoothest_file
        }
        private$smoothness_method <- "smoothest"
      } else {
        message("Using smoothness from z-statistic image itself. Not recommended!")
        private$smoothness_method <- "zstat"
      }

      # populate fwep and twosided
      private$set_fwep(fwe_p)
      private$set_twosided(two_sided)
      private$set_expected_files()

      if (!is.null(write_thresh_imgs)) {
        checkmate::assert_logical(write_thresh_imgs, len = 1L)
        private$write_thresh_imgs <- write_thresh_imgs
      }

      if (!is.null(scheduler)) {
        checkmate::assert_string(scheduler)
        checkmate::assert_subset(scheduler, c("torque", "slurm", "local"))
        private$scheduler <- scheduler
      }

      if (!is.null(time_per_zstat)) {
        private$time_per_zstat <- validate_dhms(time_per_zstat)
      }

      if (!is.null(memgb_per_command)) {
        checkmate::assert_number(memgb_per_command)
        private$memgb_per_command <- memgb_per_command
      }

    },

    #' method to return calls to external ptfce_zstat.R script for each zstat
    #' @param include_complete if TRUE, return calls for zstats that already appear to have
    #'   pTFCE-corrected images in place. Default: FALSE.
    get_ptfce_calls = function(include_complete = FALSE) {
      if (private$smoothness_method == "smoothest") {
        method_string <- glue("--fsl_smoothest {private$fsl_smoothest_file}")
      } else if (private$smoothness_method == "residuals") {
        method_string <- glue("--residuals {private$residuals_file} --dof {private$dof}")
      } else {
        method_string <- ""
      }

      if (isTRUE(private$write_thresh_imgs)) {
        write_string <- "--write_thresh_imgs"
      } else {
        write_string <- ""
      }

      if (isTRUE(private$pvt_two_sided)) {
        side_string <- "--twosided"
      } else {
        side_string <- "--onesided"
      }

      calls <- rep(NA_character_, length(private$z_files))
      script_loc <- system.file("bin/ptfce_zstat.R", package = "fmri.pipeline")
      if (script_loc == "") {
        stop("Cannot find ptfce_zstat.R inside fmri.pipeline R installation folder.")
      }

      for (pp in seq_along(private$z_files)) {
        if (isTRUE(include_complete)) {
          run_me <- TRUE
        } else {
          #check whether all expected files for a given zstat input exist
          run_me <- !checkmate::test_file_exists(unlist(private$expect_list[[pp]]))
        }

        if (run_me == FALSE) next

        calls[pp] <- glue(
          "{script_loc} --zstat {private$z_files[pp]} --mask {private$mask_files[pp]}",
          "--fwep {paste(self$fwe_p, collapse=' ')} {method_string} {write_string} {side_string}", .sep=" "
        )
      }

      return(na.omit(calls))
    },
    get_expected_files = function() {
      private$expect_list
    },
    run = function(force = FALSE) {
      # run pTFCE right here
      require(pTFCE)
      stop("At present, we only work through the external ptfce_zstat.R script. Use $submit()")
    },

    #' method to submit all ptfce_zstat.R calls to a cluster based on the scheduler specified
    #' @param force if TRUE, re-run pTFCE for zstat images that already appear to have pTFCE-corrected outputs in place
    submit = function(force = FALSE) {
      ptfce_calls <- self$get_ptfce_calls(include_complete = force)
      child_job_ids <- fmri.pipeline::cluster_submit_shell_jobs(
        ptfce_calls,
        time_per_command = private$time_per_zstat,
        memgb_per_command = private$memgb_per_command, fork_jobs = TRUE,
        scheduler = private$scheduler
      )
      return(child_job_ids)
    },
    all_expected_exist = function() {
      checkmate::test_file_exists(unlist(private$expect_list))
    }
  )

)

# for testing
# zf <- list.files(
#   pattern = "zstat[0-9]+\\.nii\\.gz", path = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-age_sex/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats",
#   full.names = TRUE
# )

# x <- ptfce_spec$new(
#   #z_files = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/zstat1.nii.gz",
#   z_files = zf,
#   mask_files = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/mask.nii.gz",
#   fsl_smoothest_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/smoothness",
#   fwe_p = c(.05, .01)
# )

# y <- ptfce_spec$new(
#   gfeat_dir = "/proj/mnhallqlab/studies/MMClock/group_analyses/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy-Intercept-Age.gfeat"
#   fwe_p = c(.05, .01)
# )
# y$submit()

