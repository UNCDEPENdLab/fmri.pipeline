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


truncate_runs <- function(s, mrfiles, mrrunnums, niftivols, drop_volumes=0) {
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

#' helper function for generating motion regressors from raw 6-parameter motion coregistration
#'
#' @param motion_params file containing 6 motion parameters
#' @param col.names names of columns in \code{motion_params}, in order from left to right
#' @param regressors names of regressors to generate and return to caller
#' @param drop_volumes number of volumes to drop from beginning of motion params
#'
#' @importFrom data.table fread
#' @importFrom checkmate assert_file_exists
generate_motion_regressors <- function(motion_params = "motion.par",
                                       col.names = c("rx", "ry", "rz", "tx", "ty", "tz"),
                                       regressors = c("rx", "ry", "rz", "tx", "ty", "tz"),
                                       demean = TRUE, drop_volumes = 0, last_volume=NULL) {

  checkmate::assert_file_exists(motion_params)
  checkmate::assert_subset(col.names, c("rx", "ry", "rz", "tx", "ty", "tz"))
  checkmate::assert_subset(regressors,
    c("rx", "ry", "rz", "tx", "ty", "tz",
    "drx", "dry", "drz", "dtx", "dty", "dtz",
    "qdrx", "qdry", "qdrz", "qdtx", "qdty", "qdtz")
  )
  
  mot <- data.table::fread(motion_params, col.names=col.names)

  if (is.null(last_volume)) { last_volume <- nrow(mot) }
  checkmate::assert_integerish(last_volume, upper=nrow(mot))
  
  mot <- mot[(1+drop_volumes):last_volume,]

  if ("fd" %in% regressors || any(derivcols <- grepl("^q?d{1}.*", regressors, perl=TRUE))) {
    motderiv <- mot[, lapply(.SD, function(x) { c(0, diff(x)) })]
    setnames(motderiv, paste0("d", names(mot))) #add delta to names
    mot <- cbind(mot, motderiv)
  }

  #quadratics always computed after derivative calculation
  if (any(quadcols <- grepl("^q{1}.*", regressors, perl=TRUE))) {
    motquad <- mot[, lapply(.SD, function(x) { x^2 })]
    setnames(motquad, paste0("q", names(mot))) #add delta to names
    mot <- cbind(mot, motquad)
  }

  if ("fd" %in% regressors) {
    #need to adapt in case of degrees (this is based on radians)
    #https://wiki.cam.ac.uk/bmuwiki/FMRI
    fd <- apply(mot.deriv[,c("drx", "dry", "drz")], 1, function(x) sum(2*pi*50*(abs(x)/360))) +
      apply(mot.deriv[,c("dtx", "dty", "dtz")], 1, function(x) sum(abs(x)))
    mot <- cbind(mot, fd=fd)
  }
  
  mot <- mot[, ..regressors] #keep regressors of interest
  if (isTRUE(demean)) {
    mot <- mot[, lapply(.SD, function(x) { x - mean(x, na.rm=TRUE) }) ] #demean all columns
  }
  
  ##just PCA motion on the current run
  ##mregressors <- pca_motion(mr_files[r], runlengths[r], motion_parfile="motion.par", numpcs=3, drop_volumes=drop_volumes)$motion_pcs_concat
  
  return(mot)
}

#' small helper function for compressing motion parameters using PCA
#'
#' @param motion_df volumes x motion parameters data frame for PCA compression
#' @param num_pcs number of principal components to extract
#' @param zscore whether to standardize motion parameters prior to PCA (this is a good idea)
#' @param verbose whether to print the variance explained by PCs and the multiple correlations of
#'   among motion parameters
#'
#' @return A matrix of PCA-compressed motion regressors
#' 
pca_motion <- function(motion_df, num_pcs=3L, zscore=TRUE, verbose=FALSE) {
  checkmate::assert_data_frame(motion_df)
  checkmate::assert_integerish(num_pcs, lower=1, upper=50)
  checkmate::assert_logical(zscore)
  checkmate::assert_logical(verbose)
  
  #compute the PCA decomposition of motion parameters and their derivatives
  pc <- prcomp(motion_df, retx=TRUE, scale.=zscore)
  cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))
  
  if (verbose) message("first", num_pcs, "motion principal components account for: ", round(cumvar[num_pcs], 3))
  mregressors <- pc$x[,1:num_pcs] #cf Churchill et al. 2012 PLoS ONE
  attr(mregressors, "variance.explained") <- cumvar[num_pcs]
  
  if (verbose) {
    cat("Multiple correlation of motion parameters:\n\n")
    print(round(sqrt(psych::smc(motion_df)), 2))
  }
  
  return(as.data.frame(mregressors))
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
