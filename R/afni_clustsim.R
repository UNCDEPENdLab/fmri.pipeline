#' R6 class for 3dClustSim automation
#' @keywords internal
#' @importFrom tibble tibble
clustsim_spec <- R6::R6Class("clustsim_spec",
  private = list(
    out_dir = getwd(),
    clustsim_df = tibble::tibble(),
    out_files = NULL,
    prefix = "clustsim_",
    pvt_use_fwhmx_acf = FALSE,
    acf_params = setNames(rep(NA_real_, 3), c("a", "b", "c")),
    clustsim_batch = NULL,
    clustsim_call = NULL,
    clustsim_mask = NULL,
    clustsim_complete = FALSE,
    method = NULL, # or 'inset' for 3dttest++ approach or 'xyz' for voxel + matrix approach
    nopad = "",
    pthr = ".01 .005 .002 .001 .0005 .0002 .0001",
    athr = ".05 .02 .01 .005 .002 .001 .0005 .0002 .0001",
    iter = 30000L,
    nodec = "",
    seed = 0,
    sim_string = "",
    scheduler = "slurm",
    ncpus = 8L,
    residuals_file = NULL,
    residuals_mask_file = NULL,
    residuals_njobs = 32, # split iter into this many jobs for permutations
    which_sim = function() {
      if (private$method == "inset") {
        private$sim_string <- glue::glue("-inset {paste(self$inset_files, collapse=' ')}")
        if (!is.null(private$clustsim_mask)) {
          private$sim_string <- paste(private$sim_string, glue::glue("-mask {private$clustsim_mask}")) # -mask compatible with -inset
        }
      } else if (private$method == "insdat") {
        private$sim_string <- glue::glue("-insdat {self$insdat_mask_file} {self$insdat_file}")

        # note that the clustsim mask controls the bounds of the simulation for FWE, 
        # while the insdat_mask is for identifying values in the insdat data file
        if (!is.null(private$clustsim_mask)) {
          private$sim_string <- paste(private$sim_string, glue::glue("-mask {private$clustsim_mask}")) # -mask compatible with -inset
        }
      } else if (!is.null(private$clustsim_mask)) {
        private$sim_string <- glue::glue("-mask {private$clustsim_mask}")
      } else {
        # use voxel size and matrix size
        private$sim_string <- glue::glue("-nxyz {paste(private$nxyz, collapse=' ')} -dxyz {paste(private$dxyz, collapse=' ')}")
      }
    },
    build_call = function() {
      if (private$method %in% c("inset", "insdat")) {
        acf_string <- "" # irrelevant
      } else if (isTRUE(private$pvt_use_fwhmx_acf)) {
        if (isFALSE(self$fwhmx_set$is_fwhmx_complete())) {
          warning("Cannot build 3dClustSim call because some 3dFWHMx runs are incomplete.")
          return(invisible(NULL))
        }
        acf_string <- glue::glue("-acf {paste(self$fwhmx_set$get_acf_average(), collapse=' ')}")
      } else {
        acf_string <- glue::glue("-acf {paste(private$acf_params, collapse=' ')}")
      }

      private$which_sim()
      private$clustsim_call <- glue::glue(
        "3dClustSim {private$sim_string} -iter {private$iter} {acf_string} -prefix {private$prefix}",
        " -pthr {private$pthr} -athr {private$athr} -seed {private$seed}{private$nopad}{private$nodec}"
      )
    },
    build_df = function() {
      if (isTRUE(private$clustsim_complete) && nrow(private$clustsim_df) > 0L) {
        return(invisible(NULL)) # data frame has already been compiled
      }

      clustsim_files <- self$get_clustsim_output_files()
      if (length(clustsim_files) == 0L) {
        # most likely, the $submit method has not completed yet.
        return(invisible(NULL))
      } else if (length(clustsim_files) == 9L) {
        private$clustsim_complete <- TRUE
      } else {
        warning("Fewer than 9 clustsim files output files detected: ", paste(clustsim_files, collapse=", "))
      }

      read_clustsim_1d <- function(file) {
        checkmate::assert_file_exists(file)
        nn <- as.integer(sub("^.*NN([1-3])_.*", "\\1", file, perl = TRUE))
        sided <- sub("^.*NN[1-3]_(1|2|bi)sided.*", "\\1", file, perl = TRUE)
        data <- read.table(file, header = FALSE)
        header <- grep(pattern = "^\\s*#.*", readLines(file), value = TRUE)
        p_line <- gsub("[#|]", "", grep(pattern = "^\\s*#\\s*pthr.*", header, value = TRUE)) # strip # and | from header line
        c_names <- strsplit(trimws(p_line), "\\s+")[[1]]
        stopifnot(length(c_names) == ncol(data))
        names(data) <- c_names
        data_long <- data %>%
          tidyr::pivot_longer(cols = !pthr, names_to = "athr", values_to = "nvoxels") %>%
          dplyr::mutate(athr = as.numeric(athr), nn = !!nn, sided = !!sided) %>%
          dplyr::select(nn, sided, pthr, athr, nvoxels)
        return(data_long)
      }

      # read each clustsim output and compile into single indexed data.frame
      private$clustsim_df <- dplyr::bind_rows(lapply(clustsim_files, read_clustsim_1d))
    }
  ),
  public = list(
    #' @field fwhmx_set A fwhmx_set_spec object containing 3dFWHMx information for all fwhmx_input_files
    fwhmx_set = NULL,

    #' @field inset_files A character vector of files to use directly as volumes to threshold and clusterize
    inset_files = NULL,

    #' @field insdat_file An sdat file containing permutations to be passed to 3dClustSim through -insdat
    insdat_file = NULL,

    #' @field insdat_file A mask file corresponding to the insdat_file data that indicates where each value is in space
    insdat_mask_file = NULL,

    #' @field only used if residuals_file is passed in, this contains the object for running the permutations
    null_3dttest_obj = NULL,

    initialize = function(out_dir = NULL, prefix=NULL,
                          fwhmx_input_files = NULL, fwhmx_mask_files = NULL, 
                          residuals_file = NULL, residuals_mask_file = NULL, residuals_njobs = NULL,
                          inset_files = NULL, insdat_file = NULL, insdat_mask_file = NULL,
                          dxyz = NULL, nxyz = NULL,
                          clustsim_mask = NULL, acf_params = NULL, nopad = NULL, pthr = NULL, athr = NULL, iter = NULL, nodec = NULL, 
                          seed = NULL, scheduler = NULL, ncpus = NULL) {

      # only one method can be used for the volume over which 3dClustSim simulates. (the first two are parts of the same method)
      n_passed <- sapply(list(dxyz, nxyz, clustsim_mask, residuals_file, inset_files, insdat_file, insdat_mask_file), function(x) !is.null(x))
      methods_passed <- 0
      if (!is.null(dxyz) || !is.null(nxyz)) methods_passed <- methods_passed + 1
      if (!is.null(inset_files)) methods_passed <- methods_passed + 1
      if (!is.null(insdat_file) || !is.null(insdat_mask_file)) methods_passed <- methods_passed + 1
      if (!is.null(clustsim_mask) && is.null(insdat_file) && is.null(inset_files)) methods_passed <- methods_passed + 1

      if (methods_passed != 1L) {
        stop("Only one method can be passed for the 3dClustSim volume settings. Use dxyz + nxyz, clustsim_mask, residuals_file, insdat_file, or inset_files")
      }

      if (!is.null(scheduler)) {
        checkmate::assert_string(scheduler)
        checkmate::assert_subset(scheduler, c("torque", "slurm", "local"))
        private$scheduler <- scheduler
      }

      if (!is.null(out_dir)) {
        checkmate::assert_directory_exists(out_dir)
        private$out_dir <- out_dir
      }

      if (!is.null(prefix)) {
        checkmate::assert_string(prefix)
        private$prefix <- prefix
      }

      if (!is.null(fwhmx_input_files)) {
        private$pvt_use_fwhmx_acf <- TRUE
        self$fwhmx_set <- fwhmx_set_spec$new(input_files = fwhmx_input_files, mask_files = fwhmx_mask_files, scheduler=scheduler)
        if (!is.null(acf_params)) {
          stop("Cannot pass fwhmx_input_files and acf_params since the ACF parameters are calculated by 3dFWHMx!")
        }
      }

      if (!is.null(inset_files)) {
        if (!is.null(fwhmx_input_files)) {
          stop("Cannot specify both inset_files and fwhmx_input_files. See 3dClustSim documentation!")
        }

        private$method <- "inset"
        checkmate::assert_file_exists(inset_files)
        self$inset_files <- inset_files
      }

      if (!is.null(insdat_file)) {
        private$method <- "insdat"
        checkmate::assert_string(insdat_file)
        checkmate::assert_file_exists(insdat_file)

        checkmate::assert_string(insdat_mask_file)
        checkmate::assert_file_exists(insdat_mask_file)
        self$insdat_file <- insdat_file
        self$insdat_mask_file <- insdat_mask_file
      }

      # if a residuals file is provided, we pass it to 3dttest++ for permutation testing
      if (!is.null(residuals_file)) {
        checkmate::assert_file_exists(residuals_file)
        private$residuals_file <- residuals_file
        private$method <- "residuals"
        if (is.null(residuals_mask_file)) {
          stop("At present, a residuals_mask_file must be passed with the residuals_file.")
        } else {
          checkmate::assert_file_exists(residuals_mask_file)
          private$residuals_mask_file <- residuals_mask_file
        }

        if (!is.null(residuals_njobs)) {
          checkmate::assert_integerish(residuals_njobs, lower = 1, upper = 1e5)
          private$residuals_njobs <- as.integer(residuals_njobs)
        }

        # note that the number of permutations should always equal the number of iterations of 3dClustSim (AFNI resets this internally anyhow)
        self$null_3dttest_obj <- simulate_null_3dttest$new(
          residuals_file = residuals_file, mask_file = residuals_mask_file,
          n_permutations = private$iter, njobs = private$residuals_njobs
        )
      }

      if (!is.null(dxyz) || !is.null(nxyz)) {
        private$method <- "xyz"
        if (!is.null(dxyz)) {
          private$dxyz <- check_nums(dxyz, lower = .01, upper = 100) # millimeters
        } else {
          stop("Spatial domain 3dClustSim specification requires both dxyz and nxyz. I only received nxyz.")
        }

        if (!is.null(nxyz)) {
          private$nxyz <- check_nums(nxyz, lower = 1, upper = 10000) # voxels
        } else {
          stop("Spatial domain 3dClustSim specification requires both dxyz and nxyz. I only received dxyz.")
        }
      }

      if (!is.null(clustsim_mask)) {
        checkmate::assert_string(clustsim_mask)
        checkmate::assert_file_exists(clustsim_mask)
        private$clustsim_mask <- clustsim_mask
      } else {
        message("No clustsim_mask argument received. All voxels in volume will be used for cluster simulations!")
      }

      if (!is.null(acf_params)) {
        checkmate::assert_numeric(acf_params, len=3)
      }

      if (!is.null(nopad)) {
        checkmate::assert_logical(nopad, len = 1L)
        if (isTRUE(nopad)) private$nopad <- " -nopad"
      }

      if (!is.null(pthr)) {
        private$pthr <- check_nums(pthr, lower = 1e-10, upper = .999)
      }

      if (!is.null(athr)) {
        private$athr <- check_nums(athr, lower = 1e-10, upper = .999)
      }

      if (!is.null(iter)) {
        checkmate::assert_integerish(iter, lower = 10, upper = 1e8)
        private$iter <- as.integer(iter)
      }

      if (!is.null(nodec)) {
        checkmate::assert_logical(nodec, len = 1L)
        if (isTRUE(nodec)) private$nodec <- " -nodec"
      }

      if (!is.null(seed)) {
        checkmate::assert_integerish(seed, len = 1L, lower=0)
        private$seed <- as.integer(seed)
      }

      if (!is.null(ncpus)) {
        checkmate::assert_integerish(ncpus, lower = 1, upper = 1e3, len = 1)
        private$ncpus <- ncpus
      }

      # always put the clustsim method in the filename for clarity
      private$prefix <- paste(private$prefix, private$method, sep = "_")

      # populate 3dClustSim object at initialization (if output files are cached already)
      private$build_df()

      # always refresh the object at initialization so that if the underlying 3dFWHMx or 3dttest++ files are complete,
      # these are current in the created clustsim object.
      self$refresh()
    },
    submit = function(force=FALSE) {
      if (isTRUE(private$clustsim_complete) && isFALSE(force)) {
        message("3dClustSim already finished for these inputs. Use $submit(force=TRUE) to re-run or $get_clustsim_df() to retrieve results.")
        return(invisible(NULL))
      }

      # create batch job
      private$clustsim_batch <- R_batch_job$new(
        job_name = "run_3dclustsim", n_cpus = private$ncpus,
        wall_time = private$wall_time, scheduler = private$scheduler,
        input_objects = list(csim_obj = self), # pass the current object as input to the batch
        r_packages = "fmri.pipeline",
        r_code = c(
          "# refresh completeness of 3dFWHMx runs in case these were run by a preceding batch job",
          "csim_obj$refresh()", # any 3dFWHMx calculations or 3dttest++ permutations are now populated
          "setwd(csim_obj$get_out_dir())",
          "run_afni_command(csim_obj$get_call(), omp_num_threads = csim_obj$get_ncpus())"
        ),
        batch_directory = private$out_dir
      )

      # need to run 3dFWHMx before 3dClustSim can run
      if (isTRUE(private$pvt_use_fwhmx_acf) && isFALSE(self$fwhmx_set$is_fwhmx_complete())) {
        fwhmx_batch <- self$fwhmx_set$get_batch()
        batch_obj <- private$clustsim_batch
        batch_obj$depends_on_parents <- "run_3dfwhmx"
        clust_seq <- R_batch_sequence$new(fwhmx_batch, batch_obj)
        clust_seq$submit()
      } else if (isTRUE(private$method == "residuals") && isFALSE(self$null_3dttest_obj$is_complete())) {
        null_batch <- self$null_3dttest_obj$get_batch()
        batch_obj <- private$clustsim_batch
        batch_obj$depends_on_parents <- "run_3dttest"
        clust_seq <- R_batch_sequence$new(null_batch, batch_obj)
        clust_seq$submit()
      } else {
        private$clustsim_batch$submit()
      }

      return(invisible(self))
    },
    get_clustsim_df = function() {
      private$build_df()
      private$clustsim_df
    },
    get_clustsim_output_files = function() {
      list.files(path = private$out_dir, pattern = glue::glue("{private$prefix}.*sided\\.1D"), full.names = TRUE)
    },
    get_call = function() {
      private$build_call()
      private$clustsim_call
    },
    get_out_dir = function() {
      private$out_dir
    },
    get_ncpus = function() {
      private$ncpus
    },
    #' return TRUE/FALSE for whether this clustsim relies on the ACF estimates from 3dFWHMx
    use_fwhmx_acf = function() {
      private$pvt_use_fwhmx_acf
    },
    is_complete = function() {
      self$refresh() # always refresh object status in case it completed through a scheduled job
      private$clustsim_complete
    },
    #' Simple method to refresh the clustsim_df, fwhmx files, and 3dttest++ permutation files
    #' This is useful if the object needs to be updated just in time to determine whether permutations or ACF params
    #' are available
    refresh = function() {
      # read clustsim output files, if available
      private$build_df()
      if (isTRUE(self$use_fwhmx_acf())) {
        self$fwhmx_set$refresh()
      }

      # if 3dttest++ permutations completed, switch method from residuals to insdat and populate files
      if (isTRUE(self$null_3dttest_obj$is_complete())) {
        if (isTRUE(self$null_3dttest_obj$get_use_sdat())) {
          private$method <- "insdat"
          ofiles <- self$null_3dttest_obj$get_permutation_files()
          self$insdat_file <- ofiles["permutation_file"]
          self$insdat_mask_file <- ofiles["mask_file"]
        } else {
          stop("Not supported yet...")
          # need to switch to -inset files
        }
      }

    }

  )
)

# for testing
# x <- simulate_null_3dttest$new(
#   residuals_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/res4d.nii.gz",
#   mask_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/mask.nii.gz",
#   n_permutations = 100000,
#   njobs = 32
# )
# x$get_3dttest_calls(include_complete = TRUE)
# x$submit()
# x$get_permutation_files()


# mytest <- clustsim_spec$new(
#   insdat_file = x$get_permutation_files()["permutation_file"], insdat_mask_file = x$get_permutation_files()["mask_file"],
#   scheduler = "slurm", prefix = "test_sdat", out_dir = "/proj/mnhallqlab/users/michael/fmri.pipeline/local",
#   clustsim_mask = "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii", ncpus = 8
# )
# mytest$submit()


# mytest <- clustsim_spec$new(
#   residuals_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/stats/res4d.nii.gz",
#   residuals_mask_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-abspe/L2m-l2_l2c-overall/L3m-int_only/FEAT_l1c-EV_abspe.gfeat/cope1.feat/mask.nii.gz",
#   residuals_njobs = 32,
#   scheduler = "slurm", prefix = "test_sdat", out_dir = "/proj/mnhallqlab/users/michael/fmri.pipeline/local",
#   clustsim_mask = "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii", ncpus = 8
# )
# mytest$submit()

#' R6 class for a list of 3dClustSim runs
#' @keywords internal
#' @importFrom tibble tibble
clustsim_list_spec <- R6::R6Class("clustsim_list_spec",
  private = list(
    clustsim_objs = NULL
  ),
  public = list(
    initialize = function(obj_list=NULL, ...) {
      if (is.null(obj_list)) {
        # assume the ... contains a set of clustsim_spec objects
        obj_list <- list(...)
      }

      class_match <- sapply(obj_list, function(x) {
        checkmate::test_class(x, "clustsim_spec")
      })

      if (!all(class_match == TRUE)) {
        stop("At least one input is not a clustsim_spec object.")
      }

      private$clustsim_objs <- obj_list
    },
    #' submit all jobs in this list
    submit = function(force = FALSE) {
      lapply(private$clustsim_objs, function(x) { x$submit(force = force) })
      return(invisible(self))
    },
    get_objs = function() {
      private$clustsim_objs
    }
  )
)
