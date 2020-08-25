#cutting down on redundancy across glm setup scripts

#wrapper for running an AFNI command safely within R
#if AFNI does not have its environment setup properly, commands may not work
runAFNICommand <- function(args, afnidir=NULL, stdout=NULL, stderr=NULL, ...) {
  #look for AFNIDIR in system environment if not passed in
  if (is.null(afnidir)) {
    env <- system("env", intern=TRUE)
    if (length(afnidir <- grep("^AFNIDIR=", env, value=TRUE)) > 0L) {
      afnidir <- sub("^AFNIDIR=", "", afnidir)
    } else {
      warning("AFNIDIR not found in environment. Defaulting to ", paste0(normalizePath("~/"), "/afni"))
      afnidir <- paste0(normalizePath("~/"), "/afni")
    }
  }
  
  Sys.setenv(AFNIDIR=afnidir) #export to R environment
  afnisetup=paste0("AFNIDIR=", afnidir, "; PATH=${AFNIDIR}:${PATH}; DYLD_FALLBACK_LIBRARY_PATH=${AFNIDIR}; ${AFNIDIR}/")
  afnicmd=paste0(afnisetup, args)
  if (!is.null(stdout)) { afnicmd=paste(afnicmd, ">", stdout) }
  if (!is.null(stderr)) { afnicmd=paste(afnicmd, "2>", stderr) }
  cat("AFNI command: ", afnicmd, "\n")
  retcode <- system(afnicmd, ...)
  return(retcode)
}


#wrapper for running an fsl command safely within R
#if FSL does not have its configuration setup properly, commands such as feat don't work, or hang strangely
runFSLCommand <- function(args, fsldir=NULL, stdout=NULL, stderr=NULL) {
  #look for FSLDIR in system environment if not passed in
  if (is.null(fsldir)) {
    #check for FSLDIR in sourced .bashrc
    bashrc_fsldir <- character(0)
    if (file.exists("~/.profile")) {
      bashrc_fsldir <- system("source ~/.profile && echo $FSLDIR", intern=TRUE)
    }
    
    #check for FSLDIR in current environment
    env <- system("env", intern=TRUE)
    if (length(fsldir <- grep("^FSLDIR=", env, value=TRUE)) > 0L) {
      fsldir <- sub("^FSLDIR=", "", fsldir)
    } else if (!identical(bashrc_fsldir, character(0))) {
      fsldir <- bashrc_fsldir      
    } else {
      warning("FSLDIR not found in environment. Defaulting to /usr/local/fsl.")
      fsldir <- "/usr/local/fsl"
    }
  }
  
  #Sys.setenv(LD_LIBRARY_PATH="/gpfs/group/mnh5174/default/sw/openblas/lib")
  Sys.setenv(FSLDIR=fsldir) #export to R environment
  fslsetup=paste0("FSLDIR=", fsldir, "; PATH=${FSLDIR}/bin:${PATH}; . ${FSLDIR}/etc/fslconf/fsl.sh; LD_LIBRARY_PATH=/gpfs/group/mnh5174/default/sw/openblas/lib ${FSLDIR}/bin/")
  fslcmd=paste0(fslsetup, args)
  if (!is.null(stdout)) { fslcmd=paste(fslcmd, ">", stdout) }
  if (!is.null(stderr)) { fslcmd=paste(fslcmd, "2>", stderr) }
  cat("FSL command: ", fslcmd, "\n")
  retcode <- system(fslcmd)
  return(retcode)
}

gen_emo_interaction_regressors <- function(examplefile, regressors, emotions=c("fear","scram","happy"), timingdir, mrrunnums, runlengths, dropVolumes) {
    ##make between-session regressors based on emotion condition.
    ##fmriDir <- "/Volumes/Serena/MMClock/MR_Raw"
    fmriDir <- "/Volumes/Serena/MMClock/MR_Proc"
    ##fmriDir <- "/Volumes/Serena/SPECC/MR_Proc"
    fitDir <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits")
    ##/Volumes/Serena/MMClock/MR_Raw/10997_20140308/MBclock_recon/clock1/nfswudktm_clock1_5_trunc282.nii.gz
    subid <- factor(sub(paste0(fmriDir, "/([0-9]{5})_\\d+/(?:mni_5mm_wavelet|native_nosmooth)/.*$"), "\\1", examplefile, perl=TRUE))
    ##subid <- as.integer(sub(paste0(fmriDir, "/([0-9]{3})[A-z]{2}_.*/mni_5mm_wavelet/.*$"), "\\1", examplefile, perl=TRUE)) #for SPECC
    
    loc <- local({load(file.path(fitDir, paste0(as.character(subid), "_fitinfo.RData"))); environment()})$f #time-clock fit object (load as local var)
    emocon <- data.frame(emotion=loc$run_condition[mrrunnums], contingency=loc$rew_function[mrrunnums]) #vector of emotion and contingency

    ##generate interactions for ev, rpe_neg, and rpe_pos with emotion
    csum <- cumsum(runlengths - dropVolumes)
    
    for (reg in regressors) {
        for (emo in emotions) {
            vec <- rep(0, max(csum))
            rmatch <- which(emocon$emotion == emo)
            for (r in rmatch) {
                f <- read.table(file.path(timingdir, paste0("run", r, "_", reg, ".1D")))$V1
                start <- ifelse(r > 1, csum[r-1] + 1, 1)
                end <- csum[r] #if (r < nrow(emocon)) csum[r] else 
                vec[start:end] <- f
            }
            write.table(vec, file=file.path(timingdir, paste0(reg, "_", emo, "_concat.1D")), row.names=FALSE, col.names=FALSE)
        }
    }
}


truncateRuns <- function(s, mrfiles, mrrunnums, niftivols, drop_volumes=0) {
  ##Identify the last valid volume acquired in a given run.
  ##Subjects often exhibit head movement after run ends (MATLAB closes), but scan hasn't stopped
  ##This occurs because the MB raw transfer of the prior run is occurring, but does not finish before the current run
  ##Thus, truncate mr files to be 12 seconds after final feedback presentation, which is how the paradigm timing files are setup
  ##note that all of this would need to be reworked if TR were not 1.0 (i.e., 1 second = 1 volume)

  require(dplyr)
  
  mrdf <- do.call(rbind, lapply(1:length(mrfiles), function(r) {
    #iti_durations <- s$runs[[ mrrunnums[r] ]]$orig_data_frame$iti_ideal #for clockfit objects
    #last_iti <- s$runs[[ mrrunnums[r] ]]$iti_onset[length(s$runs[[ mrrunnums[r] ]]$iti_onset)]
    
    iti_durations <- s %>% dplyr::filter(run == mrrunnums[r]) %>% pull(iti_ideal) #for outputs from parse_sceptic_outputs (2018+)
    last_iti <- s %>% dplyr::filter(run == mrrunnums[r]) %>% pull(iti_onset) %>% tail(n=1)

    last_vol_behavior <- floor(last_iti + iti_durations[length(iti_durations)]) #use floor to select last vol in the iti window
    first_vol <- drop_volumes #first volume to use for analysis 
    
    if (last_vol_behavior < niftivols[r]) {
      ##more vols were acquired than presented in paradigm. Thus, truncation may be needed
      ##check framewise displacement and truncate earlier than 12 second ITI if a big movement occurred
      fd <- read.table(file.path(dirname(mrfiles[r]), "motion_info", "fd.txt"))$V1
      badfd <- do.call(c, sapply(1:length(fd), function(x) { if (x >= last_iti && fd[x] > 0.9) x else NULL })) #flag volumes after last_iti with high FD
      if (length(badfd) == 0L) {
        ##no frames flagged in last volumes
        last_vol_analysis <- last_vol_behavior
      } else {
        ##use either the last volume of the task or the volume before the earliest bad movement 
        last_vol_analysis <- min(last_vol_behavior, (min(badfd) - 1))
      }
      
      #length of truncated file
      truncLength <- last_vol_analysis - first_vol
      
      #generate filename for truncated volume
      if (first_vol > 0) {
        truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$", paste0("\\1_drop", drop_volumes, "_trunc", last_vol_analysis, ".nii.gz"), mrfiles[r], perl=TRUE)  
      } else {
        truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$", paste0("\\1_trunc", last_vol_analysis, ".nii.gz"), mrfiles[r], perl=TRUE)
      }
      
      if (!file.exists(truncfile)) { runFSLCommand(paste("fslroi", mrfiles[r], truncfile, first_vol, truncLength)) } #create truncated volume
      mrfile_to_analyze <- truncfile
    } else {
      last_vol_analysis <- niftivols[r] 
      if (drop_volumes > 0) {
        truncLength <- niftivols[r] - drop_volumes
        truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$", paste0("\\1_drop", drop_volumes, ".nii.gz"), mrfiles[r], perl=TRUE)
        if (!file.exists(truncfile)) { runFSLCommand(paste("fslroi", mrfiles[r], truncfile, first_vol, truncLength)) } #create truncated volume
        mrfile_to_analyze <- truncfile
      } else {
        mrfile_to_analyze <- mrfiles[r] #just use original file  
      }
      
    }
    #cat(paste0(paste(mrfiles[r], niftivols[r], floor(last_iti), truncLength, sep="\t"), "\n"), file="trunclog", append=TRUE)
    return(data.frame(last_vol_analysis, mrfile_to_analyze, stringsAsFactors=FALSE))
  }))
  
  mrdf
  
}

pca_motion <- function(mrfiles, runlengths, motion_parfile="motion.par", verbose=FALSE,  numpcs=3, drop_volumes=0) {
  #based on a vector of mr files to be analyzed, compute the PCA decomposition of motion parameters and their derivatives
  #do this for each run separately (e.g., for FSL or R glm), as well as concatenated files
  motion_runs <- lapply(1:length(mrfiles), function(i)  {
        mot <- read.table(file.path(dirname(mrfiles[i]), motion_parfile), col.names=c("r.x", "r.y", "r.z", "t.x", "t.y", "t.z"))
        mot <- mot[(1+drop_volumes):runlengths[i],]
        motderiv <- as.data.frame(lapply(mot, function(col) { c(0, diff(col)) }))
        names(motderiv) <- paste0("d.", names(mot)) #add delta to names
        cbind(mot, motderiv)
      })
  
  motion_pcs_runs <- lapply(1:length(motion_runs), function(r) {
        pc <- prcomp(motion_runs[[r]], retx=TRUE)
        cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))
        if (verbose) message("first", numpcs, "motion principal components account for: ", round(cumvar[numpcs], 3))
        mregressors <- pc$x[,1:numpcs] #cf Churchill et al. 2012 PLoS ONE
      })
  
  motion_concat <- do.call(rbind, motion_runs)
  pc <- prcomp(motion_concat, retx=TRUE)
  cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))
  
  if (verbose) message("first", numpcs, "motion principal components account for: ", round(cumvar[numpcs], 3))
  motion_pcs_concat <- pc$x[,1:numpcs] #cf Churchill et al. 2012 PLoS ONE
  
  if (verbose) {
    cat("correlation of motion parameters:\n\n")
    print(round(cor(motion_concat), 2))
  }
  
  return(list(motion_pcs_runs=motion_pcs_runs, motion_pcs_concat=motion_pcs_concat))
  
}

generateRunMask <- function(mrfiles, outdir=getwd(), outfile="runmask") {
  if (file.exists(file.path(outdir, paste0(outfile, ".nii.gz")))) { return(invisible(NULL)) }
  ##generate mask of mrfiles where temporal min is > 0 for all runs
  for (f in 1:length(mrfiles)) {
    runFSLCommand(paste0("fslmaths ", mrfiles[f], " -Tmin -bin ", outdir, "/tmin", f))#, fsldir="/usr/local/ni_tools/fsl")
  }
  
  ##sum mins together over runs and threshold at number of runs
  runFSLCommand(paste0("fslmaths ", paste(paste0(outdir, "/tmin", 1:length(mrfiles)), collapse=" -add "), " ", outdir, "/tminsum"))#, fsldir="/usr/local/ni_tools/fsl")
  runFSLCommand(paste0("fslmaths ", outdir, "/tminsum -thr ", length(mrfiles), " -bin ", outdir, "/", outfile))#, fsldir="/usr/local/ni_tools/fsl")
  runFSLCommand(paste0("imrm ", outdir, "/tmin*"))#, fsldir="/usr/local/ni_tools/fsl") #cleanup 
  
}

visualizeDesignMatrix <- function(d, outfile=NULL, runboundaries=NULL, events=NULL, includeBaseline=TRUE) {
  require(ggplot2)
  require(reshape2)
  
  if (!includeBaseline) {
    d <- d[,!grepl("run[0-9]+base", colnames(d))]
  }
  
  print(round(cor(d), 3))
  d <- as.data.frame(d)
  d$volume <- 1:nrow(d)
  d.m <- melt(d, id.vars="volume")
  g <- ggplot(d.m, aes(x=volume, y=value)) + geom_line(size=1.2) + theme_bw(base_size=15) + facet_grid(variable ~ ., scales="free_y")
  
  colors <- c("black", "blue", "red", "orange") #just a hack for color scheme right now
  
  if (!is.null(runboundaries)) {
    g <- g + geom_vline(xintercept=runboundaries, color=colors[1L])
  }
  
  
  if (!is.null(events)) {
    for (i in 1:length(events)) {
      g <- g + geom_vline(xintercept=events[[i]], color=colors[i+1])
    }
  }
  
  if (!is.null(outfile)) {
    ggsave(filename=outfile, plot=g, width=21, height=9)
  }
  return(invisible(g))
}

#compute the mean of each cluster in roimask
#within a 4d array: x, y, z, subbrik
#
#roimask is an x,y,z array with integer mask values
#ni4d is a concatenated set of cope values, either across subjects
# or, in the case of beta series, within a subject (over trials)
get_cluster_means <- function(roimask, ni4d) {
  #ni4d should be a 4s array in which the fourth dimension is either subject (single beta) or a beta series (single subject)
  maskvals <- sort(unique(as.vector(roimask)))
  maskvals <- maskvals[!maskvals == 0]

  roimats <- sapply(maskvals, function(v) {
    mi <- which(roimask==v, arr.ind=TRUE)
    nsubbriks <- dim(ni4d)[4]
    nvox <- nrow(mi)
    mi4d <- cbind(pracma::repmat(mi, nsubbriks, 1), rep(1:nsubbriks, each=nvox))
    
    mat <- matrix(ni4d[mi4d], nrow=nvox, ncol=nsubbriks) #need to manually reshape into matrix from vector
    #for each subject, compute huber m-estimator of location/center Winsorizing at 2SD across voxels (similar to voxel mean)
    #clusavg <- apply(mat, 2, function(x) { MASS::huber(x, k=2)$mu })
    clusavg <- apply(mat, 2, mean)
    
    return(clusavg)
  })          
}

#function to extract mean beta series
get_beta_series <- function(inputs, roimask, n_bs=50) {
  #inputs <- inputs[1:5] #speed up testing

  beta_res <- foreach(i=iter(1:length(inputs)), .packages=c("reshape2", "oro.nifti", "dplyr", "abind"), .export="get_cluster_means") %dopar% {
  #beta_res <- lapply(1:length(inputs), function(i) {
    run_dirs <- list.files(path=inputs[i], pattern="FEAT_LVL1_run\\d+\\.feat", recursive=FALSE, full.names=TRUE)
    run_betas <- list()
    
    for (r in 1:length(run_dirs)) {
      runnum <- as.numeric(sub(".*/FEAT_LVL1_run(\\d+)\\.feat", "\\1", run_dirs[r], perl=TRUE))
      copes <- file.path(run_dirs[r], "stats", paste0("cope", 1:n_bs, ".nii.gz"))

      if (!all(file.exists(copes))) {
        cat("File listing in dir:", run_dirs[r], list.files(run_dirs[r], recursive=FALSE), sep="\n", append=TRUE)
        cat("Cannot find all ", length(copes), " beta series copes in: ", run_dirs[r], "\n", file="beta_series_errors.txt", append=TRUE)
        next #don't error for now, just omit
      }
      
      #in testing, it is faster to use readNIfTI to read each cope file individually (7.5s) than to use fslmerge + readNIfTI on the 4d (13s)
      #concat_file <- tempfile()
      #system.time(runFSLCommand(paste("fslmerge -t", concat_file, paste(copes, collapse=" "))))
      #system.time(cout <- readNIfTI(concat_file, reorient=FALSE)@.Data)
      
      cout <- do.call(abind, list(along=4, lapply(copes, function(x) { readNIfTI(x, reorient=FALSE)@.Data })))
      beta_series_cluster_means <- get_cluster_means(roimask, cout)
      beta_melt <- reshape2::melt(beta_series_cluster_means, value.name="bs_value", varnames=c("trial", "cluster_number")) %>%
        mutate(feat_input_id=i, run=runnum)

      run_betas[[r]] <- beta_melt
    }

    return(do.call(rbind, run_betas))
    
  }

  return(do.call(rbind, beta_res))
}


populate_sceptic_signals <- function(sceptic_signals, subj_data) {
  require(dplyr)
  
  signals <- list()
  if ("clock_bs" %in% sceptic_signals) {
    #beta series variant of clock onset
    #NB. Using RT convolution with a normalization of "none" yields a peculiar beta distribution where
    #longer responses have smaller betas because the convolved signal goes higher (scaling problem)
    signals[["clock_bs"]] <- list(event="clock", normalization="evtmax_1", value=1, beta_series=TRUE)
  }

  if ("clock" %in% sceptic_signals) {
    #clock event task regressor
    signals[["clock"]] <- list(event="clock", normalization="none", value=1)
  }

  if ("feedback" %in% sceptic_signals) {
    #feedback event task regressor
    signals[["feedback"]] <- list(event="feedback", normalization="none", value=1)
  }

  if ("feedback_bs" %in% sceptic_signals) {
    #beta series variant of feedback onset
    #NB. Using RT convolution with a normalization of "none" yields a peculiar beta distribution where
    #longer responses have smaller betas because the convolved signal goes higher (scaling problem)
    signals[["feedback_bs"]] <- list(event="feedback", normalization="evtmax_1", value=1, beta_series=TRUE)
  }

  if ("rt_swing" %in% sceptic_signals) {
    signals[["rt_swing"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, rt_swing) %>% rename(value=rt_swing))
  }

  #sqrt transform of rt swing
  if ("rt_swing_sqrt" %in% sceptic_signals) {
    signals[["rt_swing_sqrt"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, rt_swing_sqrt) %>% rename(value=rt_swing_sqrt))
  }
  
  if ("v_chosen" %in% sceptic_signals) {
    #value of chosen action, aligned with choice
    signals[["v_chosen"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_chosen) %>% rename(value=v_chosen))
  }

  #uncertainty of chosen action
  if ("u_chosen" %in% sceptic_signals) {
    #uncertainty of chosen action, aligned with choice
    signals[["u_chosen"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, u_chosen) %>% rename(value=u_chosen))
  }

  #uncertainty of chosen action with sqrt transformation
  if ("u_chosen_sqrt" %in% sceptic_signals) {
    #uncertainty of chosen action, aligned with choice
    signals[["u_chosen_sqrt"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, u_chosen_sqrt) %>% rename(value=u_chosen_sqrt))
  }

  #uncertainty of chosen action with run-level z-scoring
  if ("u_chosen_z" %in% sceptic_signals) {
    #uncertainty of chosen action, aligned with choice, z-scored
    signals[["u_chosen_z"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, u_chosen_z) %>% rename(value=u_chosen_z))
  }

  #quantile for uncertainty of chosen action
  if ("u_chosen_quantile" %in% sceptic_signals) {
    #uncertainty of chosen action, aligned with choice
    signals[["u_chosen_quantile"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, u_chosen_quantile) %>% rename(value=u_chosen_quantile))
  }
  
  #change in uncertainty of chosen action
  if ("u_chosen_change" %in% sceptic_signals) {
    #uncertainty of chosen action, aligned with choice
    signals[["u_chosen_change"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, u_chosen_change) %>% rename(value=u_chosen_change))
  }
  
  if ("v_trial_fixed" %in% sceptic_signals) {
    #value of trial according to a fixed learning rate of 0.1 (just trial - outcome)
    signals[["v_trial_fixed"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_trial_fixed) %>% rename(value=v_trial_fixed))
  }
    
  if ("v_auc" %in% sceptic_signals) {
    #value of chosen action, aligned with choice
    signals[["v_auc"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_auc) %>% rename(value=v_auc))
  }
  
  if ("v_max" %in% sceptic_signals) {
    #value of best action, aligned with choice
    signals[["v_max"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_max) %>% rename(value=v_max))
  }

  if ("rt_vmax_change" %in% sceptic_signals) {
    #absolute value of change in best RT (vmax)
    signals[["rt_vmax_change"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, rt_vmax_change) %>% rename(value=rt_vmax_change))
  }

  if ("rt_vmax_change_dir" %in% sceptic_signals) {
    #directed change in value of best RT (vmax)
    signals[["rt_vmax_change_dir"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, rt_vmax_change_dir) %>% rename(value=rt_vmax_change_dir))
  }
  
  if ("v_entropy" %in% sceptic_signals) {
    #entropy of values, computed on normalized basis weights, aligned with choice
    signals[["v_entropy"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy) %>% rename(value=v_entropy))
  }

  if ("v_entropy_1h" %in% sceptic_signals) {
    #first half entropy of values, computed on normalized basis weights, aligned with choice
    signals[["v_entropy_1h"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy_1h) %>% rename(value=v_entropy_1h))
  }

  if ("v_entropy_2h" %in% sceptic_signals) {
    #second half entropy of values, computed on normalized basis weights, aligned with choice
    signals[["v_entropy_2h"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy_2h) %>% rename(value=v_entropy_2h))
  }
  
  #drop first 5 trials
  if ("v_entropy_no5" %in% sceptic_signals) {
    #entropy of values, computed on normalized basis weights, aligned with choice
    signals[["v_entropy_no5"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy_no5) %>% rename(value=v_entropy_no5))
  }
  
  if ("v_entropy_func" %in% sceptic_signals) {
    #entropy of values, computed on discretized evaluated value function, aligned with choice
    signals[["v_entropy_func"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy_func) %>% rename(value=v_entropy_func))
  }

  if ("v_entropy_change" %in% sceptic_signals) {
    #directed change in entropy on trial t (current) versus t-1
    signals[["v_entropy_change"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy_change) %>% rename(value=v_entropy_change))
  }

  if ("v_entropy_change_pos" %in% sceptic_signals) {
    #only the positive change in entropy
    signals[["v_entropy_change_pos"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy_change_pos) %>% rename(value=v_entropy_change_pos))
  }

  if ("v_entropy_change_neg" %in% sceptic_signals) {
    #only the negative change in entropy
    signals[["v_entropy_change_neg"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy_change_neg) %>% rename(value=v_entropy_change_neg))
  }

  #ascending trial within run
  if ("run_trial" %in% sceptic_signals) {
    signals[["run_trial"]] <- list(event="clock", normalization="none",
      value=subj_data %>% select(run, trial, run_trial) %>% rename(value=run_trial))
  }
  
  if ("pe_max" %in% sceptic_signals) {
    # PE, aligned with outcome
    signals[["pe_max"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, pe_max) %>% rename(value=pe_max))
  }

  # first half PE, aligned with outcome
  if ("pe_1h" %in% sceptic_signals) {
    signals[["pe_1h"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, pe_1h) %>% rename(value=pe_1h))
  }

  # second half PE, aligned with outcome
  if ("pe_2h" %in% sceptic_signals) {
    signals[["pe_2h"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, pe_2h) %>% rename(value=pe_2h))
  }
  
  #PEs from fixed learning rate analysis
  if ("pe_trial_fixed" %in% sceptic_signals) {
    # PE from fixed 0.1 learning rate, non-SCEPTIC representation (just trial - outcome)
    signals[["pe_trial_fixed"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, pe_trial_fixed) %>% rename(value=pe_trial_fixed))
  }

  #clunky copy-paste for different learning rates -- don't feel like refactoring
  if ("pe_trial_fixed_p05" %in% sceptic_signals) {
    signals[["pe_trial_fixed_p05"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, pe_trial_fixed_p05) %>% rename(value=pe_trial_fixed_p05))
  }

  if ("pe_trial_fixed_p10" %in% sceptic_signals) {
    signals[["pe_trial_fixed_p10"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, pe_trial_fixed_p10) %>% rename(value=pe_trial_fixed_p10))
  }

  if ("pe_trial_fixed_p15" %in% sceptic_signals) {
    signals[["pe_trial_fixed_p15"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, pe_trial_fixed_p15) %>% rename(value=pe_trial_fixed_p15))
  }

  if ("pe_trial_fixed_p20" %in% sceptic_signals) {
    signals[["pe_trial_fixed_p20"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, pe_trial_fixed_p20) %>% rename(value=pe_trial_fixed_p20))
  }
  
  if ("rew_om" %in% sceptic_signals) {
    signals[["rew_om"]] <- list(event="feedback", normalization="none",
      value=subj_data %>% select(run, trial, rew_om) %>% rename(value=rew_om))
  }
  
  if ("d_auc" %in% sceptic_signals) {
    # decay AUC, aligned with outcome
    signals[["d_auc"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, d_auc) %>% rename(value=d_auc))
  }

  if ("d_auc_clock" %in% sceptic_signals) {
    # decay AUC, aligned with choice
    signals[["d_auc_clock"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, d_auc) %>% rename(value=d_auc))
  }
  
  if ("d_auc_sqrt" %in% sceptic_signals) {
    # decay AUC, aligned with outcome
    signals[["d_auc_sqrt"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, d_auc_sqrt) %>% rename(value=d_auc_sqrt))
  }

  #align entropy at feedback, not clock
  if ("v_entropy_feedback" %in% sceptic_signals) {
    #entropy of values, computed on normalized basis weights, aligned with choice
    signals[["v_entropy_feedback"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, v_entropy) %>% rename(value=v_entropy))
  }

  #K-L distance measures
  if ("intrinsic_discrepancy" %in% sceptic_signals) {
    #intrinsic discrepancy measure of transition in value vector (basis weights) on t-1 vs. t
    signals[["intrinsic_discrepancy"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, intrinsic_discrepancy) %>% rename(value=intrinsic_discrepancy))
  }

  if ("intrinsic_discrepancy_feedback" %in% sceptic_signals) {
    #feedback-aligned intrinsic discrepancy measure of transition in value vector (basis weights) on t-1 vs. t
    signals[["intrinsic_discrepancy_feedback"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, intrinsic_discrepancy) %>% rename(value=intrinsic_discrepancy))
  }

  if ("mean_kld" %in% sceptic_signals) {
    #intrinsic discrepancy measure of transition in value vector (basis weights) on t-1 vs. t
    signals[["mean_kld"]] <- list(event="clock", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, mean_kld) %>% rename(value=mean_kld))
  }

  if ("mean_kld_feedback" %in% sceptic_signals) {
    #feedback-aligned intrinsic discrepancy measure of transition in value vector (basis weights) on t-1 vs. t
    signals[["mean_kld_feedback"]] <- list(event="feedback", normalization="evtmax_1",
      value=subj_data %>% select(run, trial, mean_kld) %>% rename(value=mean_kld))
  }

  return(signals)
}
