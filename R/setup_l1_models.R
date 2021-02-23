#' Function for setting up level 1 model specifications and corresponding files
#'
#' @param gpa A glm_model_arguments function setup by the \code{setup_glm_pipeline} function
#' @param to_setup A character vector of model names within \code{gpa$l1_models} whose l1
#'    inputs should be generated. If omitted, all models within \code{gpa$l1_models} will be generated.
#'
#' @details The \code{to_setup} argument allows the creation of multiple l1 models to be parallelized at
#'   a superordinate level or, even better, to be spread across independent jobs on a cluster. This function
#'   already provides the option to parallelize over subjects for a single model if \code{gpa$l1_setup_cpus}
#'   is greater than 1.
#'
#' @importFrom checkmate assert_class assert_character
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom parallel makeCluster stopCluster
#' @importFrom RNifti niftiHeader
setup_l1_models <- function(gpa, to_setup=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(to_setup, null.ok=TRUE)
  checkmate::assert_data_frame(gpa$subject_data)
  
  #if no model subset is requested, output all models
  if (is.null(to_setup)) { to_setup <- names(gpa$l1_models$models) }
 
  lg <- lgr::get_logger("glm_pipeline/l1_setup")
  if (isTRUE(gpa$log_json)) { lg$add_appender(lgr::AppenderJson$new(paste0(gpa$l1_setup_log, ".json"), name="json")) }
  if (isTRUE(gpa$log_txt)) { lg$add_appender(lgr::AppenderFile$new(paste0(gpa$l1_setup_log, ".txt"), name="txt")) }
  #if (isTRUE(gpa$log_console)) { lg$add_appender(lgr::AppenderConsole$new()) }
  
## model_clock_fmri_lvl1 <- function(trial_statistics, id_col=NULL, subject_data=NULL, drop_volumes=6, ncpus=1,
##                                   expectdir="mni_5mm_aroma", expectfile = "nfaswuktm_clock[0-9]_5.nii.gz",
##                                   sceptic_run_signals=c("v_chosen", "v_entropy", "d_auc", "pe_max"), #which signals to model jointly in LVL1
##                                   l1_contrasts=NULL,
##                                   outdir=NULL, glm_software="fsl", ...) {

  
  
  #setup parallel worker pool, if requested
  if (gpa$parallel$l1_setup_cores > 1L) {
    lg$info("Initializing l1 setup cluster with %d cores", gpa$parallel$l1_setup_cores)
    cl <- parallel::makeCluster(gpa$parallel$l1_setup_cores)
    doParallel::registerDoParallel(cl)    
    on.exit(try(stopCluster(cl))) #cleanup pool upon exit of this function
  } else {
    lg$info("Initializing l1 setup with serial execution")
    foreach::registerDoSEQ() #formally register a sequential 'pool' so that dopar is okay
  }

  idvec <- gpa$subject_data[[ gpa$vm["id"] ]]
    
  # loop over each subject, identify relevant fMRI data, and setup FSL level 1 files
  #by_subj <- by_subj[4]
  ll <- foreach(subid = iter(idvec), .inorder=FALSE, .packages=c("dependlab"),
    .export=c("lg", "gpa", "truncateRuns", "fsl_l1_model", "spm_l1_model", "runFSLCommand") ) %dopar% {

      browser()
      #use specific run NIfTIs included in run_data, rather finding these by regex
      if (gpa$vm["run_nifti"] %in% names(gpa$run_data)) {
        lg$info("Using run_data to identify NIfTI files for analysis")
        subj_runs <- gpa$run_data %>% filter(!!sym(gpa$vm["id"]) == !!subid)
        if (nrow(subj_runs) != 1L) {
          print(subj_runs);
          log$fatal("Found an unexpected number of runs (%d) for subject %s", nrow(subj_runs), subid)
          stop("Unexpected number of runs in subj_runs")
        }
        
        mrfiles <- subj_runs[[ gpa$vm["run_nifti"] ]]

        mr_found <- file.exists(mrfiles)
        if (any(mr_found != TRUE)) {
          lg$debug("Cannot find any files: %s", paste(mrfiles, collapse=", "))
          mrdirs <- subj_runs[[ gpa$vm["mr_dir"] ]]
          mrfiles[!mr_found] <- file.path(mrdirs[!mr_found], mrfiles[!mr_found]) #try prepending the mr_dir path if run_nifti is relative
        }

        mr_found <- file.exists(mrfiles)
        if (any(mr_found != TRUE)) {
          lg$warn("Could not find the following run files: %s", paste(mrfiles[!mr_found], collapse=", "))
        }

        lg$debug(paste("MR file to analyze:", mrfiles[mr_found])) #log files that were found
        
      } else {
        lg$info("Using regex-based find approach to identify run NIfTIs")
        
        #find run nifti files based on directory and regular expression settings
        mrmatch <- gpa$subject_data %>% filter(!!sym(gpa$vm["id"]) == !!subid) %>% pull(!!gpa$vm["mr_dir"])

        if (length(mrmatch) != 1L) {
          lg$warn("Unable to find fMRI directory record in subject_data for subid: %s", subid)
          return(NULL)
        }
        
        if (!dir.exists(file.path(mrmatch))) {
          lg$warn("Unable to find subject data directory: %s for subid: %s", mrmatch, subid)
          return(NULL)
        }
        
        ## Find processed fMRI run-level data for this subject
        #mrfiles <- list.files(mrmatch, pattern=gpa$fmri_file_regex, full.names=TRUE, recursive=TRUE)
        ##cat(paste0("command: find ", mrmatch, " -iname '", expectfile, "' -ipath '*", expectdir, "*' -type f\n"))

        # -ipath '*", expectdir, "*' -type f | sort -n"), intern=TRUE)
        if (!is.null(gpa$fmri_path_regex)) { addon <- paste0(" -ipath '*/", gpa$fmri_path_regex, "/*'") } else { addon <- "" }
        find_string <- paste0("find ", mrmatch, " -regextype posix-egrep -iregex '.*", gpa$fmri_file_regex, "'", addon, " -type f | sort -n")
        lg$debug("mrfiles find syntax: %s", find_string)
        mrfiles <- system(find_string, intern=TRUE)

        ##mr_run_nums <- as.integer(sub(paste0(".*", expectfile, "$"), "\\1", mrfiles, perl=TRUE))
        mr_run_nums <- as.integer(sub(paste0(gpa$run_number_regex), "\\1", mrfiles, perl=TRUE)) #extract run number from file name

        ##NB. If we reorder the mrfiles, then the run numbers diverge unless we sort(mr_run_nums). Remove for now for testing
        ##mrfiles <- mrfiles[order(mr_run_nums)] #make absolutely sure that runs are ordered ascending

      }

      #process files
      if (length(mrfiles) == 0L) {
        lg$warn("Unable to find any preprocessed fMRI files in dir: %s", mrmatch)
        return(NULL)
      }

      ##read number of volumes from NIfTI header
      #RNifti is unexpectedly slow
      #run_lengths <- unname(sapply(mrfiles, function(x) { RNifti::niftiHeader(x)$dim[5L] }))
      run_lengths <- unname(sapply(mrfiles, function(x) { oro.nifti::readNIfTI(x, read_data=FALSE)@dim_[5L] }))
      lg$debug("Run lengths of mrfiles: %s", paste(run_lengths, collapse=", "))
      
      ## create truncated run files to end analysis 12s after last ITI (or big head movement)
      ## also handle removal of N volumes from the beginning of each run due to steady state magnetization
      
      #mrdf <- truncate_runs(b, mrfiles, mr_run_nums, run_lengths, drop_volumes=drop_volumes)
      mrdf <- data.frame(mrfile_to_analyze=mrfiles, run=mr_run_nums, last_vol_analysis=run_lengths, drop_volumes=gpa$drop_volumes)
      
      mrfiles <- mrdf$mrfile_to_analyze
      run_lengths <- mrdf$last_vol_analysis

      #loop over models to output
      for (this_model in to_setup) {
        #setup design matrix for any given software package
        m_signals <- gpa$l1_models$signals[gpa$l1_models$models[[this_model]]$model_signals]

        t_out <- gpa$glm_software
        if (isTRUE(gpa$use_preconvolve)) { t_out <- c("convolved", t_out) } #compute preconvolved regressors
        ba <- gpa$additional$ba
        ba$events <- gpa$l1_models$events
        ba$signals <- m_signals
        ba$tr <- gpa$tr
        ba$t_out <- t_out
        ba$drop_volumes <- gpa$drop_volumes
        ba$run_volumes <- run_lengths
        ba$runs_to_output <- mr_run_nums
        d <- do.call(build_design_matrix, ba)
        
      }
      
      if ("fsl" %in% gpa$glm_software) {
        #Setup FSL run-level models for each combination of signals
        tryCatch(fsl_l1_model(b, sceptic_run_signals, l1_contrasts, mrfiles, run_lengths, mrrunnums, drop_volumes=drop_volumes, outdir=outdir, ...),
          error=function(e) {
            cat("Subject: ", b[[id_col]][1], ", run variant: ", paste(sceptic_run_signals, collapse="-"), " failed with mrfiles: \n",
              paste(mrfiles, collapse="\n"), "\n", "error: ", as.character(e), "\n\n", file="lvl1_crashlog.txt", append=TRUE)
          })
      }

      if ("spm" %in% gpa$glm_software) {
        #Setup spm run-level models for each combination of signals
        tryCatch(spm_l1_model(b, sceptic_run_signals, mrfiles, run_lengths, mrrunnums, drop_volumes=drop_volumes, outdir=outdir, ...),
          error=function(e) {
            cat("Subject: ", b[[id_col]][1], ", run variant: ", paste(sceptic_run_signals, collapse="-"), " failed with mrfiles: \n",
              paste(mrfiles, collapse="\n"), "\n", "error: ", as.character(e), "\n\n", file="lvl1_crashlog.txt", append=TRUE)
          })

      }

      message("completed processing of subject: ", subid)
      cat("\n\n\n")
    }
  
}
