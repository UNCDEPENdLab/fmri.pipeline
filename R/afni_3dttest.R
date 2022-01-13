simulate_null_3dttest <- R6::R6Class("simulate_null_3dttest",
  private = list(
    residuals_file = NULL,
    mask_file = NULL,
    n_permutations = 20000,
    njobs = 16,
    job_output_files = NULL, # full path to permutation outputs for each job
    combined_output_file = NULL, # output file for combined permutations across all jobs
    o_prefix = NULL, # overall filename/prefix for job outputs
    split_iter = NULL, # vector of how many sign flips to perform for each job
    use_sdat = TRUE,
    scheduler = "slurm",
    method = "residuals", # in future, potentially support full 3dttest++ inputs like -setA and so on
    perm_batch = NULL,
    wall_time = "4:00:00", # default of 4 hours to complete all permutations once job begins
    memgb_per_3dttest = "8", # 8 GB RAM per cpu (job)
    memgb_combine = "32g", # 32 GB RAM for the combination of permutations step (the parent job)

    # handy one-liner from: https://stackoverflow.com/questions/64828789/creating-a-function-to-split-single-number-in-approximately-equal-groups-r
    int_split = function(n, p) n %/% p + (sequence(p) - 1 < n %% p),
    generate_batch = function(include_complete = TRUE) {
      # this will create an R_batch_job object for running all 3dttest++ calls

      perm_calls <- self$get_3dttest_calls(include_complete = include_complete)
      if (length(perm_calls) == 0L) {
        message("No 3dttest++ calls are returned")
        private$perm_batch <- NULL # reset to NULL for consistency
        return(invisible(NULL))
      }

      # TODO: need to handle situation where the output file exists -- 3dtoXdataset will just concatenate on top of it
      if (isTRUE(private$use_sdat)) {
        if (is.null(private$mask_file)) {
          stop("A mask file must be provided for sdat outputs.")
        }
        combine_call <- glue("3dtoXdataset -prefix {private$combined_output_file} {private$mask_file}")
      } else {
        combine_call <- glue("3dTcat -overwrite -prefix {private$combined_output_file}")
      }

      # create parent batch job to run all 3dFWHMx scripts
      private$perm_batch <- R_batch_job$new(
        job_name = "run_3dttest", n_cpus = 1, mem_total = private$memgb_combine,
        wall_time = private$wall_time, scheduler = private$scheduler,
        input_objects = list(perm_calls = perm_calls, job_output_files = private$job_output_files), # export this object to the job
        wait_for_children = TRUE, r_packages = "fmri.pipeline", repolling_interval = 60, # check completion on child jobs every minute
        r_code = c(
          glue(
            "child_job_ids <- cluster_submit_shell_jobs(perm_calls, memgb_per_command={private$memgb_per_3dttest}, ",
            "fork_jobs=TRUE, time_per_job='{private$wall_time}', scheduler='{private$scheduler}', job_script_prefix='job_3dttest')"
          )
        ),

        # code to combine permutations after the jobs complete
        post_children_r_code = c(
          glue("run_afni_command(paste(\"{combine_call}\", paste(job_output_files, collapse=' ')))"),
          "unlink(job_output_files, sub('\\\\.nii\\\\.gz$', '.minmax.1D', job_output_files))" # job_output_files shuttled to batch job environment
        )
      )

      return(private$perm_batch)
    }
  ),
  #' @param residuals_file the 4D file containing voxelwise residuals for all subjects (e.g., res4d.nii.gz in FEAT)
  #' @param mask_file A mask file used to specify which voxels should be analyzed/permuted
  #' @param njobs The number of independent jobs across which permutations are distributed
  #' @param n_permutations The total number of null datasets to be computed by sign-flipping
  #' @param use_sdat a logical indicating whether to output null datasets in sdat format (single-precision, serialized, I think)
  public = list(
    initialize = function(residuals_file = NULL, mask_file = NULL, njobs = NULL, n_permutations = NULL, use_sdat = NULL,
                          wall_time = NULL, memgb_per_3dttest = NULL, memgb_combine = NULL) {

      checkmate::assert_string(residuals_file) # for now, must be single input
      checkmate::assert_file_exists(residuals_file)
      private$residuals_file <- residuals_file

      # for now, we just support this method
      private$method <- "residuals"

      if (!is.null(use_sdat)) {
        checkmate::assert_logical(use_sdat, len = 1L)
        private$use_sdat <- use_sdat
      }

      if (isTRUE(private$use_sdat) && is.null(mask_file)) {
        stop("A mask_file must be provided for permutations to be output in an sdat format!")
      }

      if (!is.null(mask_file)) {
        checkmate::assert_string(mask_file)
        checkmate::assert_file_exists(mask_file)
        private$mask_file <- mask_file
      }

      if (!is.null(njobs)) {
        checkmate::assert_integerish(njobs, lower = 1, upper = 1e4)
        private$njobs <- as.integer(njobs)
      }

      if (!is.null(n_permutations)) {
        checkmate::assert_integerish(n_permutations, lower = 1, upper = 1e8)
        if (n_permutations < 100) {
          warning("Number of permutations chosen is less than 100! It's hard to know why this is a good idea!")
        }
        private$n_permutations <- as.integer(n_permutations)
      }

      if (!is.null(wall_time)) {
        private$wall_time <- validate_dhms(wall_time)
      }

      if (!is.null(memgb_per_3dttest)) {
        if (checkmate::test_string(memgb_per_3dttest)) {
          memgb_per_3dttest <- as.integer(memgb_per_3dttest)
        }
        checkmate::assert_integerish(memgb_per_3dttest, lower = 0.1, upper = 128, len = 1L)
        private$memgb_per_3dttest <- memgb_per_3dttest
      }

      if (!is.null(memgb_combine)) {
        if (checkmate::test_string(memgb_combine)) {
          memgb_combine <- as.integer(memgb_combine)
        }
        checkmate::assert_integerish(memgb_combine, lower = 0.1, upper = 128, len = 1L)
        private$memgb_combine <- paste0(memgb_combine, "g")
      }

      # populate expected files and settings
      private$o_prefix <- paste(file_sans_ext(private$residuals_file), "randomsign", sep = "_")
      private$combined_output_file <- paste0(private$o_prefix, ifelse(isTRUE(private$use_sdat), ".sdat", ".nii.gz"))

      # generate how many permutations per job
      private$split_iter <- private$int_split(private$n_permutations, private$njobs)

      # generate expected output file name for each job
      private$job_output_files <- sapply(seq_along(private$split_iter), function(i) {
        paste0(private$o_prefix, "_", sprintf("%04d", i), ".nii.gz") # always force nifti output
      })
    },

    is_complete = function() {
      checkmate::test_file_exists(private$combined_output_file)
    },
    submit = function(force = FALSE) {
      if (isTRUE(self$is_complete()) && isFALSE(force)) {
        message(glue("The permutations output file already exists: {private$combined_output_file}."))
        return(invisible(self))
      }
      private$generate_batch()
      # private$perm_batch$generate()
      private$perm_batch$submit()
    },
    get_3dttest_calls = function(include_complete = FALSE) {
      if (private$method == "residuals") {
        # 3dttest++ uses all caps for SDAT output and lower case for nii/BRIK output
        # at present, there is not an AFNI program to concatenate SDAT files, so we can't use this
        # if (isTRUE(private$use_sdat)) {
        #   fmt_str <- "-RANDOMSIGN"
        # } else {
        #   fmt_str <- "-randomsign"
        # }

        fmt_str <- "-randomsign" # nifti output

        mask_str <- ifelse(is.null(private$mask_file), "", glue("-mask {private$mask_file}"))

        # based on afni src: https://github.com/afni/afni/blob/9b6398061ab25afc95a14eed7f86eb9d1b1deb37/src/3dttest%2B%2B.c#L5335
        base <- glue(
          "3dttest++ -DAFNI_AUTOMATIC_FDR=NO -DAFNI_DONT_LOGFILE=YES ",
          "-nomeans -toz -setA {private$residuals_file} {mask_str}"
        )

        all_calls <- sapply(seq_along(private$split_iter), function(i) {
          glue(base, " {fmt_str} {private$split_iter[i]} -prefix {private$job_output_files[i]}")
        })

        # filter out complete intermediate files if requested
        if (isFALSE(include_complete)) {
          f_exist <- sapply(private$job_output_files, file.exists)
          all_calls <- all_calls[!f_exist]
        }

        return(all_calls)
      } else {
        stop("Not implemented!")
      }
    },
    get_permutation_files = function() {
      if (!self$is_complete()) {
        stop("Cannot provide the permutation file because it doesn't yet exist! Try $submit()")
      }

      c(mask_file = private$mask_file, permutation_file = private$combined_output_file)
    },
    #' method to return R_batch_job to run 3dttest permutations
    get_batch = function() {
      private$generate_batch()
      return(private$perm_batch)
    }
  )
)
