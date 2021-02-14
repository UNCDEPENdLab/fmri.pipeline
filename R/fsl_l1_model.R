#note: this is a small adapation from the original fslSCEPTICModel to avoid use of the clockfit objects and to move to the
#simpler build_design_matrix approach and the use of the trial_statistics csv files from vba_fmri
fsl_l1_model <- function(subj_data, sceptic_signals, l1_contrasts=NULL, mrfiles, runlengths, mrrunnums, execute_feat=FALSE, force=FALSE,
                              drop_volumes=0, outdir=NULL, usepreconvolve=FALSE, spikeregressors=FALSE, model_suffix="", ...) {

  # subj_data is the trial-level data for one subject, as produced by parse_sceptic_outputs
  # sceptic_signals is a character vector of column names in subj_data used for parametric modulator regressors
  # mrfiles is a character vector of processed data to analyze
  # runlengths is the number of volumes in each run
  # mrrunnums is the numeric vector of the run numbers corresponding to mrfiles
  # execute_feat specifies whether to execute feat for the runs after setting up the fsf files
  # drop_volumes specifies how many volumes were dropped from the beginning of the run for elements of mrfiles.
  #   This is used by build_design_matrix to ensure that the convolved regressors line up properly with the fMRI data
  # outdir is the base name of folder for the specified model. The resulting directories are nested inside the data folder for the subject
  # usepreconvolve is a TRUE/FALSE denoting whether to use convolved regressors from build_design_matrix (TRUE) or let FSL handle 3-column format convolution (FALSE)
  # any additional arguments trapped by ... are passed forward to build_design_matrix  
  
  require(Rniftilib)
  require(dplyr)
  require(tidyr)
  require(dependlab)
  
  if (is.null(outdir)) {
    outdir=paste0("sceptic-", paste(sceptic_signals, collapse="-")) #define output directory based on combination of signals requested
    if (usepreconvolve) { outdir=paste(outdir, "preconvolve", sep="-") }

    outdir <- paste0(outdir, model_suffix) #add any model suffix, if requested
  }

  #determine which feat template is relevant
  use_new <- TRUE #generate EV and contrast syntax dynamically
  if (use_new) {
    fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl1_clock_sceptic_nparam_template.fsf"))
  } else {    
    if (length(sceptic_signals) == 1L) {
      ##single model-based regressor
      if (usepreconvolve) {
        fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl1_clock_sceptic_univariate_preconvolve_template.fsf"))
      } else {
        fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl1_clock_sceptic_univariate_template.fsf"))
      }
    } else if (length(sceptic_signals) == 4L) { #pemax, dauc, vchosen, ventropy
      if (usepreconvolve) {
        fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl1_clock_sceptic_4param_preconvolve_template.fsf"))
      } else {
        stop("not implemented yet")
      }
    } else if (length(sceptic_signals) == 5L) { #pemax, dauc, vchosen, ventropy, vtime
      if (usepreconvolve) {
        fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl1_clock_sceptic_5param_preconvolve_template.fsf"))
      } else { stop("not implemented yet") }      
    } else { stop("not implemented yet") }
  }
  
  #note: normalizePath will fail to evaluate properly if directory does not exist
  fsl_run_output_dir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), outdir)

  if (file.exists(fsl_run_output_dir) && force==FALSE) { message(fsl_run_output_dir, " exists. Skipping."); return(0) }
  cat("fsl_run_output_dir create: ", fsl_run_output_dir, "\n")
  dir.create(fsl_run_output_dir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(fsl_run_output_dir, "run_timing_sceptic")

  #create the events and signals structures for the build_design_matrix call
  #this assumes that we have received a data.frame with the structure from the _trial_statistics.csv.gz generated
  #thus, we have $clock_onset and $feedback_onset, and $iti_onset available
  events <- subj_data %>% dplyr::select(id, run, trial, clock_onset, feedback_onset, iti_onset, rt_csv) %>%
    dplyr::mutate(clock_duration=rt_csv/1000, feedback_duration=iti_onset - feedback_onset) %>%
    dplyr::select(-iti_onset, -rt_csv) %>% tidyr::gather(key="key", value="value", -id, -run, -trial) %>%
    tidyr::separate(col = key, into = c("event", "onset_duration")) %>%
    tidyr::spread(key=onset_duration, value=value) %>% dplyr::select(event, run, trial, onset, duration)

  signals <- populate_sceptic_signals(sceptic_signals, subj_data)
  
  #not currently handling vtime
  # else if (thisName == "vtime") {
  #      #vtime is a runs x trials list with a data.frame per trial containing onsets, durations, and values within trial
  #      onsets[v] <- NA #not relevant
  #      durations[v] <- NA #not relevant
  #      normalizations[v] <- "none" #should not try to normalize the within-trial regressor since this starts to confused within/between trial variation

  #save(file=file.path(fsl_run_output_dir, "bdm_call.RData"), events, signals, timingdir, drop_volumes, mrfiles, mrrunnums)
  
  #NB. The tr argument should be passed in as part of ...
  d <- build_design_matrix(events=events, signals=signals, baseline_coef_order=2, write_timing_files = c("convolved"), #, "FSL"),
    center_values=TRUE, plot=FALSE, convolve_wi_run=TRUE, output_directory=timingdir, drop_volumes=drop_volumes,
    run_volumes=mrfiles, runs_to_output=mrrunnums, ...)

  save(d, subj_data, events, signals, timingdir, runlengths, mrrunnums, file=file.path(fsl_run_output_dir, "designmatrix.RData"))
  
  allFeatFiles <- list()
  
  #FSL computes first-level models on individual runs
  for (r in 1:length(mrfiles)) {
    stopifnot(file.exists(file.path(dirname(mrfiles[r]), "motion.par"))) #can't find motion parameters
    
    runnum <- sub("^.*/clock(\\d+)$", "\\1", dirname(mrfiles[r]), perl=TRUE)
    nvol <- nifti.image.read(mrfiles[r], read_data=0)$dim[4L]

    ##just PCA motion on the current run
    ##mregressors <- pca_motion(mrfiles[r], runlengths[r], motion_parfile="motion.par", numpcs=3, drop_volumes=drop_volumes)$motion_pcs_concat

    ##Add volumes to censor here. Use censor_intersection.mat, which flags fd > 0.9 and DVARS > 20
    ##15Jun2016: Switch to FD > 0.9mm censoring in general (moving away from wavelet)
    ##If fd_0.9.mat doesn't exist, it means no spike regressors were generated at this threshold
    ##Thus, do not include in the nuisance set. Also do not include PCA motion regressors
    ##censorfile <- file.path(dirname(mrfiles[r]), "motion_info", "censor_intersection.mat")
    ##if (file.exists(censorfile) && file.info(censorfile)$size > 0) {
    ##  censor <- read.table(censorfile, header=FALSE)$V1
    ##  censor <- censor[(1+drop_volumes):runlengths[r]]
    ##  mregressors <- cbind(mregressors, censor)
    ##}

    mregressors <- NULL #start with NULL

    if (spikeregressors) { #incorporate spike regressors if requested (not used in conventional AROMA)
      censorfile <- file.path(dirname(mrfiles[r]), "motion_info", "fd_0.9.mat")
      if (file.exists(censorfile) && file.info(censorfile)$size > 0) {
        censor <- read.table(censorfile, header=FALSE)
        censor <- censor[(1+drop_volumes):runlengths[r],,drop=FALSE] #need no drop here in case there is just a single volume to censor
        #if the spikes fall outside of the rows selected above, we will obtain an all-zero column. remove these
        censor <- censor[,sapply(censor, sum) > 0,drop=FALSE]
        if (ncol(censor) == 0L) { censor <- NULL } #no volumes to censor within valid timepoints
        mregressors <- censor
      }
    }
    
    ##add CSF and WM regressors (with their derivatives)
    nuisancefile <- file.path(dirname(mrfiles[r]), "nuisance_regressors.txt")
    if (file.exists(nuisancefile)) {
      nuisance <- read.table(nuisancefile, header=FALSE)
      nuisance <- nuisance[(1+drop_volumes):runlengths[r],,drop=FALSE]
      nuisance <- as.data.frame(lapply(nuisance, function(col) { col - mean(col) })) #demean
      #cat("about to cbind with nuisance\n")
      #print(str(mregressors))
      #print(str(nuisance))
      if (!is.null(mregressors)) { mregressors <- cbind(mregressors, nuisance) #note that in R 3.3.0, cbind with NULL or c() is no problem...
      } else { mregressors <- nuisance }
    }
    
    motfile <- file.path(fsl_run_output_dir, paste0("run", runnum, "_confounds.txt"))
    write.table(mregressors, file=motfile, col.names=FALSE, row.names=FALSE)

    #search and replace within fsf file for appropriate sections
    ##.OUTPUTDIR. is the feat output location
    ##.NVOL. is the number of volumes in the run
    ##.FUNCTIONAL. is the fmri data to process (sans extension)
    ##.CONFOUNDS. is the confounds file for GLM
    ##.CLOCK_TIMES. is the three-column file for clock onset
    ##.FEEDBACK_TIMES. is the three-column file for feedback onset
    ##.VNAME. is the signal name in a univariate model
    ##.V_TIMES. is the three-column file for the signal
    ##.V_CON. is the contrast name for the signal
    
    thisTemplate <- fsfTemplate
    thisTemplate <- gsub(".OUTPUTDIR.", file.path(fsl_run_output_dir, paste0("FEAT_LVL1_run", runnum)), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".NVOL.", nvol, thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".FUNCTIONAL.", gsub(".nii(.gz)*$", "", mrfiles[r]), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".CONFOUNDS.", motfile, thisTemplate, fixed=TRUE)

    if (use_new) {
      #generate ev syntax
      dmat <- d$design_convolved[[paste0("run", runnum)]] %>% select(-matches("base\\d+")) #drop baseline columns
      regressors <- as.list(names(dmat))

      #add common ingredients for preconvolved regressors
      regressors <- lapply(regressors, function(x) { list(name=x, waveform="custom_1", convolution="none", tempfilt=1, timing_file=file.path(timingdir, paste0("run", runnum, "_", x, ".1D"))) })

      ev_syn <- dependlab::generate_fsf_lvl1_ev_syntax(regressors)

      #creation of l1 contrast matrices, including the diagonal contrasts, now abstracted to finalize_pipeline_configuration.R
      #thus, l1_contrasts is already a contrast matrix ready to be passed to the generate_fsf_contrast_syntax function
      cmat_syn <- dependlab::generate_fsf_contrast_syntax(l1_contrasts)
      
      thisTemplate <- c(thisTemplate, ev_syn, cmat_syn)      
    } else {
      if (usepreconvolve) {
        thisTemplate <- gsub(".CLOCK_TIMES.", file.path(timingdir, paste0("run", runnum, "_clock.1D")), thisTemplate, fixed=TRUE)
        thisTemplate <- gsub(".FEEDBACK_TIMES.", file.path(timingdir, paste0("run", runnum, "_feedback.1D")), thisTemplate, fixed=TRUE)
      } else {
        thisTemplate <- gsub(".CLOCK_TIMES.", file.path(timingdir, paste0("run", runnum, "_clock_FSL3col.txt")), thisTemplate, fixed=TRUE)
        thisTemplate <- gsub(".FEEDBACK_TIMES.", file.path(timingdir, paste0("run", runnum, "_feedback_FSL3col.txt")), thisTemplate, fixed=TRUE)
      }

      for (s in 1:length(sceptic_signals)) {
        if (usepreconvolve) {
          thisTemplate <- gsub(paste0(".V", s, "_TIMES."), file.path(timingdir, paste0("run", runnum, "_", sceptic_signals[s], ".1D")), thisTemplate, fixed=TRUE)
        } else {
          thisTemplate <- gsub(paste0(".V", s, "_TIMES."), file.path(timingdir, paste0("run", runnum, "_", sceptic_signals[s], "_FSL3col.txt")), thisTemplate, fixed=TRUE)
        }
        thisTemplate <- gsub(paste0(".V", s, "NAME."), sceptic_signals[s], thisTemplate) #define EV name
        thisTemplate <- gsub(paste0(".V", s, "_CON."), sceptic_signals[s], thisTemplate) #define contrast name
      }

    }
    
    featFile <- file.path(fsl_run_output_dir, paste0("FEAT_LVL1_run", runnum, ".fsf"))
    if (file.exists(featFile) && force==FALSE) { next } #skip re-creation of FSF and do not run below unless force==TRUE 
    cat(thisTemplate, file=featFile, sep="\n")
    
    allFeatFiles[[r]] <- featFile
  }

  #if execute_feat is TRUE, execute feat on each fsf files at this stage, using an 8-node socket cluster (since we have 8 runs)
  #if execute_feat is FALSE, just create the fsf files but don't execute the analysis
  if (execute_feat == TRUE) {
    require(parallel)
    cl_fork <- makeForkCluster(nnodes=8)
    runfeat <- function(fsf) {
      runname <- basename(fsf)
      runFSLCommand(paste("feat", fsf), stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
    }
    clusterApply(cl_fork, allFeatFiles, runfeat)
    stopCluster(cl_fork)
  }
  
}
