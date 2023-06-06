#!/usr/bin/env Rscript

# run like
# ./ROI_TempCorr.R -ts /Volumes/Phillips/WorkingMemory/11155/20130418/SpatialWM_v1_run1_768x720.5/nfswkmtd_functional.nii.gz -brainmask  /Volumes/Phillips/gm_50mask.nii 

printHelp <- function() {
  cat("ROI_TempCorr is a script that computes temporal correlations among ROIs defined by an integer-valued mask (e.g., a set of numbered spheres).",
      "",
      "Required inputs are:",
      "  -ts       <4D file>: The file containing time series of interest. Usually preprocessed resting-state data. Can be NIfTI or AFNI BRIK/HEAD",
      "  -rois <3D file>: The integer-valued mask file defining the set of ROIs. Computed correlation matrices will be nROIs x nROIs in size based on this mask.",
      "  -out_file <filename for output>: The file to be output containing correlations among ROIs.",
      "",
      "Optional arguments are:",
      "  -pcorr_method <pearson|spearman|kendall|adalasso.net|lasso.net|pls.net|ridge.net|pcor.shrink>: Method to compute partial correlation",
      "      pearson, spearman, and kendall compute the corresponding partial correlations (pcor function in ppcor package). BEWARE: Uses pseudoinverse for rank-deficient matrices.",
      "      adalasso.net uses adaptive LASSO regression with cross-validation to compute partial correlations (adalasso.net function in parcor package)",
      "      lasso.net uses LASSO regression (no adaptation) with cross-validation to compute partial correlations (adalasso.net function in parcor package)",
      "      pls.net uses partial least squares regression models to estimate partial correlations (pls.net function in parcor package)",
      "      ridge.net uses optimal ridge regression models to estimate partial correlations (ridge.net function in parcor package",
      "  -pcorr_cvsplit <10>: For adalasso.net, lasso.net, pls.net, and ridge.net methods, the number of splits for k-fold cross-validation.",      
      "  -corr_method <pearson|spearman|kendall|robust|mcd|weighted|donostah|M|pairwiseQC|pairwiseGK|cor.shrink>: Method to compute correlations among time series.",
      "      Default: pearson is the standard Pearson correlation",
      "      spearman is Spearman correlation based on ranks",
      "      kendall is Kendall's tau correlation (also based on ranks, but with somewhat better statistical properties than Spearman)",
      "      robust uses the covRob function from the robust package with estim=\"auto\" to obtain robust estimate of correlation (reduce sensitivity to outliers)",
      "      mcd, weighted, donostah, M, pairwiseQC, pairwiseGK are different robust estimators of correlation. See ?covRob in the robust package for details.",
      "      cor.shrink computes a shrinkage estimate of the correlation matrix by shrinking correlations towards the identity matrix. See ?cor.shrink in the corpcor package.",
      "  -roi_reduce <pca|median|mean|huber>: Method to obtain a single time series for voxels within an ROI. Default: pca",
      "      pca takes the first eigenvector within the ROI, representing maximal shared variance",
      "      median uses the median within each ROI",
      "      mean uses the mean within each ROI",
      "      huber uses the huber robust estimate of the center of the distribution (robust to outliers)",
      "  -brainmask <3D file>: A 0/1 mask file defining voxels in the brain. This will be applied to the ROI mask before computing correlations.",
      "  -dropvols <integer=0>: The number of initial volumes to drop from all voxel time series prior to computing correlations (e.g., for steady state problems)",
      "  -censor <1D file>: An AFNI-style 1D censor file containing a single column of 0/1 values where 0 represents volumes to be censored (e.g., for motion scrubbing)",
      "  -fisherz: Apply Fisher's z transformation (arctanh) to normalize correlation coefficients. Not applied by default.",
      "  -njobs <n>: Number of parallel jobs to run when computing correlations. Default: 4.",
      "  -na_string: Character string indicating how to represent missing correlations in output file. Default NA.",
      "  -ts_out_file <filename for time series output>: Output a file containing the average time series for each region before computing correlations.",
      "  -ts_only: stop before running correlations. useful with -ts_out_file",
      "  -prewhiten_arima <p> <d> <q>: prewhiten all time series by fitting an ARIMA model of order p (autoregressive), d (integrated), q (moving average) before computing correlation.",
      "  -detrend_ts <0,1,2>: Demean (0), linear detrending (1) or quadratic detrending (2) is applied to ROI aggregated time series before correlations are computed. Default: 2",
      "  -nuisance_regressors <matrix.txt>: If passed in, the nuisance regressors (columns) of this white space-separated file will be regressed out of ROI aggregated time series prior to correlations.",
      "  -bandpass_filter dt freq_low freq_high: apply a FIR1 bandpass filter to ROI averaged data and nuisance regressors (if specified).",
      "      dt is the repetition time in seconds, freq_low is the lowest frequency to pass in Hz, freq_high is the highest frequency to pass in Hz",
      "  -roi_diagnostics <file>: Create a file with information about the number of voxels per ROI in the mask versus the data, including the proportion missing.",
      "  -write_header 1/0: Whether to include a header row with variable names in the output. Variables are named according to their corresponding mask values. Default: 0",
      "",
      "If the -ts file does not match the -rois file, the -ts file will be resampled to match the -rois file using 3dresample. This requires that the images be coregistered,",
      "  in the same stereotactic space, and have the same grid size.",
      "",
      "Multiple correlation methods can be computed at once by passing a comma-separated list. Both partial and full correlations computation can also be combined.",
      "  Example: -corr_method pearson,cor.shrink -pcor_method pcor.shrink,adalasso.net",
      "",
      "The script depends on the following R libraries: foreach, doParallel, MASS, oro.nifti, parcor, ppcor, pracma, and robust. These can be installed using:",
      "  install.packages(c(\"foreach\", \"doParallel\", \"MASS\", \"oro.nifti\", \"parcor\", \"ppcor\", \"pracma\", \"ppcor\", \"robust\", \"forecast\", \"lmtest\"))",
      sep="\n")
}

#read in command line arguments.
args <- commandArgs(trailingOnly = FALSE)

scriptpath <- dirname(sub("--file=", "", grep("--file=", args, fixed=TRUE, value=TRUE), fixed=TRUE))
argpos <- grep("--args", args, fixed=TRUE)
if (length(argpos) > 0L) {
  args <- args[(argpos+1):length(args)]
} else {
  args <- c()
}

#contains run_afni_command
source(file.path(scriptpath, "R_helper_functions.R"))

if (is.null(args) || length(args) == 0L) {
  message("ROI_TempCorr expects at least -ts <4D file> -rois <3D file> -out_file <filename for output>.\n")
  printHelp()
  quit(save="no", 1, FALSE)
}

#defaults
njobs <- 4
out_file <- "corr_rois.txt"
ts_out_file <- ""
fname_censor1D <- NULL
fname_brainmask <- NULL
cm <- "pearson"
cm_partial <- FALSE #TRUE for partial correlation methods
pm <- c()
pm_partial <- c()
roi_reduce <- "mean" #see Varoquax and Craddock 2013 for why mean is a better choice than first eigenvariate
roi_diagnostics_fname <- NULL
fisherz <- FALSE
drop_vols <- 0
pcorr_cvsplit <- 10 #default 10-fold cross-validation for parcor functions
ts_only <- FALSE # only compute timeseries, not correlation matrix
write_header <- FALSE #whether to include mask value as a header row in outputs
na_string <- "NA"
clustersocketport <- 11290 #default port for parallel cluster
fit_p <- NULL #autoregressive order for ARIMA
fit_d <- NULL #differencing order for ARIMA
fit_q <- NULL #moving average order from ARIMA
delta_t <- NULL #time in seconds, used for bandpass filtering
freq_low <- NULL #filter low frequency cutoff in Hz
freq_high <- NULL #filter high frequency cutoff in Hz
nuisance_regressors <- NULL #filename of nuisance regressors to be removed from ROI time series
detrend_order <- 2 #deterministic detrending order (0, 1, or 2) to be projected from ROI time series
roi_arima_fits <- NULL #contains the fits for ARIMA models to each time series, if relevant
white_lags <- 6 #hard code for now: the number of lags to test with Breusch-Godfrey whiteness test in ARIMA

#for testing
#fname_rsproc <- "/gpfs/group/mnh5174/default/MMClock/MR_Proc/11336_20141204/mni_nosmooth_aroma_hp/rest1/Abrnawuktm_rest1.nii.gz" #name of preprocessed fMRI data
#fname_rsproc <- "/gpfs/group/mnh5174/default/MMClock/MR_Proc/11336_20141204/mni_5mm_3ddespike/rest1/brnswudktm_rest1_5.nii.gz" #name of preprocessed fMRI data
#fname_rsproc <- "nawuktm_rest1.nii.gz" #name of preprocessed fMRI data
#fname_roimask <- "/gpfs/group/mnh5174/default/lab_resources/parcellation/final_combined/7_Networks/Schaefer421_full_final_revised2018_groupmask95.nii.gz"
#fname_roimask <- "Schaefer_422_final_jul2018.nii.gz"
#fname_brainmask <- "/gpfs/group/mnh5174/default/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii"
#fname_brainmask <- "MNI152_T1_2.3mm_brain_mask.nii"
#nuisance_regressors <- "/gpfs/group/mnh5174/default/MMClock/MR_Proc/11336_20141204/mni_5mm_3ddespike/rest1/nuisance_regressors.txt"
#nuisance_regressors <- "unfiltered_nuisance_regressors.txt" #name of preprocessed fMRI data
#nuisance_regressors <- NULL
#fname_censor1D <- "/gpfs/group/mnh5174/default/MMClock/MR_Proc/10637_20140304/mni_5mm_wavelet/clock1/motion_info/censor_union.1D"
#fit_p <- 3
#fit_q <- 2
#fit_d <- 0
#delta_t = 1
#freq_low = .009
#freq_high = 1
#detrend_order = 2
#roi_reduce = "mean"

argpos <- 1
while (argpos <= length(args)) {
  if (args[argpos] == "-ts") {
    fname_rsproc <- args[argpos + 1] #name of preprocessed fMRI data
    stopifnot(file.exists(fname_rsproc))
    argpos <- argpos + 2
  } else if (args[argpos] == "-rois") {
    fname_roimask <- args[argpos + 1] #name of integer-valued ROI mask file
    argpos <- argpos + 2
    stopifnot(file.exists(fname_roimask))
  } else if (args[argpos] == "-out_file") {
    out_file <- args[argpos + 1] #name of file to be written
    argpos <- argpos + 2
  } else if (args[argpos] == "-ts_out_file") {
    ts_out_file <- args[argpos + 1] #name of file to be written
    argpos <- argpos + 2
  } else if (args[argpos] == "-censor") {
    fname_censor1D <- args[argpos + 1] #name of censor file
    argpos <- argpos + 2
    stopifnot(file.exists(fname_censor1D))
  } else if (args[argpos] == "-brainmask") {
    fname_brainmask <- args[argpos + 1] #mask file for brain voxels
    argpos <- argpos + 2
    stopifnot(file.exists(fname_brainmask))    
  } else if (args[argpos] == "-njobs") {
    njobs <- as.integer(args[argpos + 1])
    argpos <- argpos + 2
    if (is.na(njobs)) { stop("-njobs must be an integer") }
  } else if (args[argpos] == "-roi_diagnostics") {
    roi_diagnostics_fname <- args[argpos + 1]
    argpos <- argpos + 2
  } else if (args[argpos] == "-roi_reduce") {
    roi_reduce <- args[argpos + 1]
    argpos <- argpos + 2
    stopifnot(roi_reduce %in% c("pca", "mean", "median", "huber"))
  } else if (args[argpos] == "-corr_method") {
    cm <- strsplit(args[argpos + 1], ",")[[1]]
    cm_partial <- rep(FALSE, length(cm))
    argpos <- argpos + 2
    stopifnot(all(cm %in% c("pearson", "spearman", "robust", "kendall", "mcd", "weighted", "donostah", "M", "pairwiseQC", "pairwiseGK", "cor.shrink")))
  } else if (args[argpos] == "-pcorr_method") {
    pm <- strsplit(args[argpos + 1], ",")[[1]]
    pm_partial <- rep(TRUE, length(pm))
    argpos <- argpos + 2
    stopifnot(all(pm %in% c("pearson", "spearman", "kendall", "adalasso.net", "lasso.net", "pls.net", "ridge.net", "pcor.shrink")))
  } else if (args[argpos] == "-pcorr_cvsplit") {
    pcorr_cvsplit <- as.numeric(args[argpos + 1])
    if (is.na(pcorr_cvsplit)) { stop("Could not understand argument ", args[argpos+1], "to -pcorr_cvsplit") }
    argpos <- argpos + 2
  } else if (args[argpos] == "-fisherz") {
    fisherz <- TRUE
    argpos <- argpos + 1
  } else if (args[argpos] == "-ts_only") {
    ts_only <- TRUE
    argpos <- argpos + 1
  } else if (args[argpos] == "-na_string") {
    na_string <- args[argpos + 1]
    argpos <- argpos + 2
  } else if (args[argpos] == "-dropvols") {
    drop_vols <- as.numeric(args[argpos + 1]) #number of vols to drop
    if (is.na(drop_vols)) { stop("Could not understand argument ", args[argpos+1], "to -dropvols") }
    argpos <- argpos + 2
  } else if (args[argpos] == "-port") {
    clustersocketport <- as.integer(args[argpos + 1])
    argpos <- argpos + 2
    if (is.na(njobs)) { stop("-port must be an integer") }
  } else if (args[argpos] == "-prewhiten_arima") {
    fit_p <- as.integer(args[argpos + 1])
    fit_d <- as.integer(args[argpos + 2])
    fit_q <- as.integer(args[argpos + 3])
    if (any(is.na(c(fit_p, fit_d, fit_q)))) { stop ("You must pass three integer arguments (p, d, q) to -prewhiten_arima") }
    argpos <- argpos + 4
  } else if (args[argpos] == "-nuisance_regressors") {
    nuisance_regressors <- args[argpos + 1]
    stopifnot(file.exists(nuisance_regressors))
    argpos <- argpos + 2
  } else if (args[argpos] == "-bandpass_filter") {
    delta_t <- as.numeric(args[argpos + 1]) #repetition time (sampling rate) in seconds
    f_low <- as.numeric(args[argpos + 2]) #low frequency cutoff in Hz
    f_high <- as.integer(args[argpos + 3]) #high frequency cutoff in Hz
    
    if (any(is.na(c(delta_t, f_low, f_high)))) { stop ("You must pass three numeric arguments (dt, freq_low, freq_high) to -bandpass_filter") }
    argpos <- argpos + 4
  } else if (args[argpos] == "-detrend_ts") {
    detrend_order <- as.integer(args[argpos + 1])
    
    if (is.na(detrend_order) || ! detrend_order %in% c(1,2)) { stop ("-detrend_ts order must be an integer (1 or 2)") }
    argpos <- argpos + 2
  } else if (args[argpos] == "-write_header") {
    write_header <- as.integer(args[argpos + 1])
    argpos <- argpos + 2
    if (is.na(write_header) || (!write_header %in% c(0,1))) { stop("-write_header must be 1 or 0")
    } else { write_header <- as.logical(write_header) }
  } else {
    stop("Not sure what to do with argument: ", args[argpos])
  }
}

corr_method <- c(cm, pm) #put together full and partial methods
partial <- c(cm_partial, pm_partial)

#robust package uses "auto" to choose best robust estimator given problem complexity (matrix size)
corr_method <- ifelse(corr_method == "robust", "auto", corr_method)

# check for input sanity before spending time loading things in
if (ts_only && nchar(ts_out_file) == 0L) {
  stop('-ts_only is set without -ts_out_file! No point in running.')
}

# handle package dependencies
for (pkg in c("methods", "foreach", "doParallel", "oro.nifti", "MASS", "corpcor", "parcor", "tictoc", "pracma", "forecast", "lmtest")) {
  if (!suppressMessages(require(pkg, character.only=TRUE))) {
    message("Installing missing package dependency: ", pkg)
    install.packages(pkg)
    suppressMessages(require(pkg, character.only=TRUE))
  }
}

if (!is.null(fname_censor1D)) {
  stopifnot(file.exists(fname_censor1D))
  censor1D <- read.table(fname_censor1D, header=FALSE)$V1
  censor_vols <- which(censor1D == 0.0)      
} else {
  censor_vols <- c()
}

if (!is.null(nuisance_regressors)) {
  message("Reading white space-separated nuisance regressors from: ", nuisance_regressors)
  nuisance_df <- as.matrix(read.table(nuisance_regressors, header=FALSE))
  message("Nuisance file has: ", nrow(nuisance_df), " timepoints (rows) and ", ncol(nuisance_df), " regressors (columns)")
}

#generate (robust or partial) correlation matrix given a set of time series.
genCorrMat <- function(roits, method="auto", fisherz=FALSE, partial=FALSE, roi_arima_fits=NULL) {
  #roits should be an time x roi data.frame
  
  suppressMessages(require(robust))
  
  #assume that parallel has been setup upstream
  njobs <- getDoParWorkers()
  
  #sapply only works for data.frame
  if (!inherits(roits, "data.frame")) stop("genCorrMat only works properly with data.frame objects.")
  
  #remove missing ROI columns for estimating correlation
  nacols <- which(sapply(roits, function(col) all(is.na(col))))
  if (length(nacols) > 0) nona <- roits[,nacols*-1]
  else nona <- roits
  
  #All partial correlation methods are necessarily full matrix methods (since partial correlation is conditioned on remaining variables)
  #Standard correlations are also computed on full matrices. Only the robust estimators of correlation blow up on singular matrices
  if (!is.null(roi_arima_fits)) {
    pairwise <- TRUE
  } else if (method %in% c("pearson", "spearman", "kendall", "cor.shrink", "pcor.shrink", "adalasso.net", "lasso.net", "pls.net", "ridge.net")) {
    pairwise <- FALSE
  } else {
    #estimators using covRob from robust package tend to expect positive definite matrix
    pairwise <- TRUE
  }
  
  if (!pairwise) {
    if (method=="cor.shrink") {
      message("Estimating shrinkage estimates (always positive definite) of correlation matrix using cor.shrink")
      rcorMat <- corpcor::cor.shrink(as.matrix(nona)) #omit lambda to estimate shrinkage using analytic formula
    } else if (method=="pcor.shrink") {
      message("Estimating shrinkage estimates (always positive definite) of partial correlation matrix using pcor.shrink")
      rcorMat <- corpcor::pcor.shrink(as.matrix(nona)) #omit lambda to estimate shrinkage using analytic formula
    } else if (method=="adalasso.net") {
      rcorMat <- parcor::adalasso.net(as.matrix(nona), k=pcorr_cvsplit, both=TRUE)$pcor.adalasso
    } else if (method=="lasso.net") {
      rcorMat <- parcor::adalasso.net(as.matrix(nona), k=pcorr_cvsplit, both=FALSE)$pcor.lasso
    } else if (method=="pls.net") {
      rcorMat <- parcor::pls.net(as.matrix(nona), k=pcorr_cvsplit)$pcor
    } else if (method=="ridge.net") {
      rcorMat <- parcor::ridge.net(as.matrix(nona), k=pcorr_cvsplit)$pcor
    } else if (partial) {
      stopifnot(require(ppcor))
      message("ppcor pcor func")
      rcorMat <- ppcor::pcor(as.matrix(nona), method=method)
    } else {
      #conventional cor using specified method
      rcorMat <- cor(as.matrix(nona), method=method)
    }
    
    #adalasso.net(X, k = 10,use.Gram=FALSE,both=TRUE,verbose=FALSE,intercept=TRUE)
    #ridge.net(X, lambda, plot.it = FALSE, scale = TRUE, k = 10,verbose=FALSE)
    #pls.net(X, scale = TRUE, k = 10, ncomp = 15,verbose=FALSE)
    
    #come back to this: diagnostics on partial correlation?
    #if (partial) {
    #  require(parcor)
    #  #at the moment, GeneNet does not expose ggm.test.edges, so set fdr=FALSE to avoid crash
    #  perf <- parcor::performance.pcor(rcorMat, true.pcor=NULL, fdr=FALSE, cutoff.ggm=0.8, verbose=FALSE, plot.it=FALSE)
    #  print(perf)
    #}
  } else {
    
    #Due to rank degeneracy of many RS-fcMRI roi x time matrices, correlations are sometimes filled in pairwise.
    #This is slow, of course, but necessary.
    rcorMat <- matrix(NA, nrow=ncol(nona), ncol=ncol(nona))
    diag(rcorMat) <- 1
    
    #indices of lower triangle
    lo.tri <- which(lower.tri(rcorMat), arr.ind=TRUE)
    
    # how many chucks per core
    # for small datasets, need to make sure we didn't pick too high a number
    chunksPerProcessor <- 8
    
    # if we only have one job, we want all the chunks together
    # that is we want chunksize == nrow(lo.tri)
    if(njobs == 1) { chunksPerProcessor <- 1 }
    
    repeat {
      chunksize <- floor(nrow(lo.tri)/njobs/chunksPerProcessor)
      if(chunksize >= 1) break
      # decrease chunk size and check we can still go lower
      chunksPerProcessor <- chunksPerProcessor - 1
      cat('WARNING: job is small, decreasing chunks to',chunksPerProcessor,'\n')
      if (chunksPerProcessor < 1) { stop('too many jobs for too little work, lower -n') }
    }
    
    # let us know when we are asking for a lot of work
    # the threshold is less for non-full cor types
    chunksizewarn <- 2000
    if (chunksize > chunksizewarn) message(sprintf('Lots of datapoints (%d) going to each processor. Consider increasing -njobs if this is slow',chunksize))
    
    corrpair <- function(ts1, ts2, method="pearson") {
      if (method %in% c("robust", "mcd", "weighted", "donostah", "M", "pairwiseQC", "pairwiseGK")) {
        robust::covRob(na.omit(cbind(ts1, ts2)), estim=method, corr=TRUE)$cov[1,2]
      } else {
        cor(na.omit(cbind(ts1, ts2)), method=method)[1,2]
      }
    }
    
    #do manual chunking: divide correlations across processors, where each processor handles 10 chunks in total (~350 corrs per chunk)
    if (method %in% c("robust", "mcd", "weighted", "donostah", "M", "pairwiseQC", "pairwiseGK")) {
      tictoc::tic("Estimating pairwise robust correlation estimates")
      corrvec <- foreach(pair=iter(lo.tri, by="row", chunksize=chunksize), .inorder=TRUE, .combine=c, .multicombine=TRUE, .packages="robust") %dopar% {
        #iter will pass entire chunk of lo.tri, use apply to compute row-wise corrs
        apply(pair, 1, function(x) { corrpair(nona[,x[1]], nona[,x[2]]) })
      }
      tictoc::toc()
    } else if (!is.null(roi_arima_fits)) {
      #prewhiten ROI A based on ROI B
      tictoc::tic("Estimating pairwise correlation estimates after prewhitening using ARIMA")
      corrvec <- foreach(pair=iter(lo.tri, by="row", chunksize=chunksize), .inorder=TRUE, .combine=cbind, .multicombine=TRUE) %dopar% {
        #iter will pass entire chunk of lo.tri, use apply to compute row-wise corrs
        apply(pair, 1, function(x) {
          cvec <- rep(0, 2) #filter with target as x, y, or just regular old correlation
          for (filter_direction in 1:2) {
            response <- ifelse(filter_direction==1, 1, 2)
            input <- ifelse(filter_direction==1, 2, 1)
            filter_model <- roi_arima_fits[[ x[response] ]]
            ts_response <- residuals(filter_model)
            
            mcoefs <- coef(filter_model)
            armanames <- grep("(ar\\d+|ma\\d+)", names(mcoefs), value=TRUE)
            newcoefs <- setNames(rep(NA_real_, length(mcoefs)), names(mcoefs))
            newcoefs[armanames] <- mcoefs[armanames]
            
            #fix the ARMA coefficients, leave the rest of the model freely estimated.
            #note that include.drift the first time around will add a column called 'drift' to $xreg.
            xreg <- filter_model$xreg
            if ("drift" %in% colnames(xreg)) { 
              xreg <- xreg[,grep("^drift$", colnames(xreg), value=TRUE, invert=TRUE)]
              if (ncol(xreg) == 0) { xreg = NULL } #if drift was the only column (no exogenous covariates)
            }

            input_filtered <- tryCatch(forecast::Arima(nona[, x[input]], order=forecast::arimaorder(filter_model), fixed=newcoefs, xreg=xreg, include.mean=TRUE, include.drift=TRUE, method="CSS-ML"),
              error=function(e) { print(e); forecast::Arima(nona[, x[input]], order=forecast::arimaorder(filter_model), fixed=newcoefs, xreg=xreg, include.mean=TRUE, include.drift=TRUE, method="ML") })
              
            ts_input <- residuals(input_filtered)
            cvec[filter_direction] <- corrpair(residuals(input_filtered), ts_response, method=method)
            #cvec[filter_direction] <- cor(residuals(input_filtered), ts_response)
          }
          return(cvec)
        })
      }
      tictoc::toc()
      
      #for now, average the estimates of correlations from the two filter directions
      message("For ARIMA filtering, we prewhiten in each direction (i.e., alternating input and response series), then average the lag-0 cross-correlation")
      message("Similarity of cross-correlations alternating the input and response ARIMA filters: r = ", round(cor(corrvec[1,], corrvec[2,]),3))
      corrvec <- apply(corrvec, 2, mean)
    }
    
    stopifnot(length(corrvec)==nrow(lo.tri)) #be afraid is we haven't estimated every cell in the correlation matrix
    rcorMat[lo.tri] <- corrvec
    #duplicate the lower triangle to upper
    rcorMat[upper.tri(rcorMat)] <- t(rcorMat)[upper.tri(rcorMat)] #transpose flips filled lower triangle to upper
  }
  
  #populate lower triangle of correlation matrix
  if (fisherz == TRUE) { 
    message("Applying the Fisher z transformation to correlation coefficients.")
    rcorMat <- atanh(rcorMat)
    diag(rcorMat) <- 15 #avoid Inf in output by a large value. tanh(15) is ~1
  }
  
  #add back in NA cols
  if (length(nacols) > 0) {
    processedCorrs <- matrix(NA, nrow=ncol(roits), ncol=ncol(roits))
    
    #fill in all non-NA cells row-wise
    processedCorrs[!1:ncol(roits) %in% nacols, !1:ncol(roits) %in% nacols] <- rcorMat
  } else { processedCorrs <- rcorMat }
  
  processedCorrs #return processed correlations
  
}

### BEGIN DATA PROCESSING
message("Reading roi mask: ", fname_roimask)
if (grepl("^.*\\.(HEAD|BRIK|BRIK.gz)$", fname_roimask, perl=TRUE)) {
  roimask <- readAFNI(fname_roimask, vol=1)
  #afni masks tend to read in as 4D matrix with singleton 4th dimension. Fix this
  if (length(dim(roimask)) == 4L) {
    roimask@.Data <- roimask[,,,,drop=T]    
  }
} else {
  roimask <- readNIfTI(fname_roimask, reorient=FALSE)
}

#optional: apply brain mask
if (!is.null(fname_brainmask)) {
  message("Applying brain mask to ROIs: ", fname_brainmask)
  stopifnot(file.exists(fname_brainmask))
  if (grepl("^.*\\.(HEAD|BRIK|BRIK.gz)$", fname_brainmask, perl=TRUE)) {
    brainmask <- readAFNI(fname_brainmask)
  } else {
    brainmask <- readNIfTI(fname_brainmask, reorient=FALSE)
  }
  
  #brain mask and roi mask must be of same dimension
  stopifnot(identical(dim(brainmask)[1:3], dim(roimask)[1:3]))
  
  roimask_prebrainmask <- roimask #for tracking filtering
  message("  ROI voxels before applying brainmask: ", sum(roimask > 0, na.rm=TRUE))
  
  #remove non-brain voxels
  roimask[which(brainmask == 0.0)] <- NA_real_
  
  message("  ROI voxels after applying brainmask:  ", sum(roimask > 0, na.rm=TRUE))
  
}

message("Reading in 4D file: ", fname_rsproc)
#read in processed resting-state data
if (grepl("^.*\\.(HEAD|BRIK|BRIK.gz)$", fname_rsproc, perl=TRUE)) {
  rsproc <- readAFNI(fname_rsproc)
} else {
  rsproc <- readNIfTI(fname_rsproc, reorient=FALSE)
}

if (!identical(dim(rsproc)[1:3], dim(roimask)[1:3])) {
  message("Resampling rs proc file from: ", paste(dim(rsproc)[1:3], collapse="x"), " to: ", paste(dim(roimask)[1:3], collapse="x"), " using nearest neighbor")
  message("This assumes that the files are in the same space and have the same grid size. Make sure this is what you want!!")
  
  run_afni_command(paste0("3dresample -overwrite -inset ", fname_rsproc, " -rmode NN -master ", fname_roimask, " -prefix tmpResamp.nii.gz"))
  stopifnot(file.exists("tmpResamp.nii.gz"))
  
  rsproc <- readNIfTI("tmpResamp.nii.gz", reorient=FALSE)
  unlink("tmpResamp.nii.gz")
}

#obtain vector of mask values 
maskvals <- sort(unique(as.vector(roimask)))
maskvals <- maskvals[which(maskvals != 0)] #omit zero

if (length(maskvals) > 1000) {
  warning("More than 1000 putative ROIs identified in mask file: ", fname_roimask)
}

if (njobs > 1) {
  clusterobj <- makePSOCKcluster(njobs, master="localhost", port=clustersocketport)
  registerDoParallel(clusterobj)
} else {
  registerDoSEQ()
}

#to reduce RAM overhead of having to copy rsproc_censor to each worker, obtain list of vox x time mats for rois

#even though this seems more elegant, it is much slower (400x!) than the use of 4d lookup and reshape below
#system.time(roimats <- lapply(maskvals, function(v) {
#      apply(rsproc, 4, '[', which(roimask==v, arr.ind=TRUE))
#    }))

#generate a 4d mat of indices
roimats <- lapply(maskvals, function(v) {
  mi <- which(roimask==v, arr.ind=TRUE)
  nvol <- dim(rsproc)[4]
  nvox <- nrow(mi)
  mi4d <- cbind(pracma::repmat(mi, nvol, 1), rep(1:nvol, each=nvox))
  mat <- matrix(rsproc[mi4d], nrow=nvox, ncol=nvol) #need to manually reshape into matrix from vector
  attr(mat, "maskval") <- v #add mask value as attribute so that information about bad ROIs can be printed below
  t(mat) #transpose matrix so that it is time x voxels
})

#add numeric index of mask (1:max) for keeping track of number of good voxels inside roimats loop
#this is necessary if there is discontinuity in the mask values (e.g., jumps between 10 and 12)
roimats <- lapply(1:length(roimats), function(i) {
  attr(roimats[[i]], "maskindex") <- i
  return(roimats[[i]])
})

rm(rsproc) #clear imaging file from memory now that we have obtained the roi time series 

if (!is.null(roi_diagnostics_fname)) {
  nvox_observed_per_roi <- sapply(roimats, ncol)
  nvox_good_per_roi <- rep(NA, length(nvox_observed_per_roi))
}

message("Obtaining a single time series within each ROI using: ", roi_reduce)
roiavgmat <- foreach(roivox=iter(roimats), .packages=c("MASS"), .combine=cbind, .noexport=c("rsproc")) %do% { #minimal time savings from dopar here, and it prevents message output
  ##roivox is a time x voxels matrix for a single ROI
  ##data cleaning steps: remove voxels that are 1) partially or completely missing; 2) all 0; 3) variance = 0 (constant)
  ##leave out variance > mean check because bandpass-filtered data are demeaned
  badvox <- apply(roivox, 2, function(voxts) {
    if (any(is.na(voxts))) TRUE #any missing values
    else if (all(voxts == 0.0)) TRUE #all zeros
    else if (var(voxts) == 0.0) TRUE #constant time series
    ##else if (var(voxts) > mean(voxts)) TRUE #variance exceeds mean (very unstable heuristic)
    else FALSE #good voxel
  })
  
  if (!is.null(roi_diagnostics_fname)) { nvox_good_per_roi[attr(roivox, "maskindex")] <- sum(!badvox) }
  
  if (sum(!badvox) < 5) {
    ##only reduce if there are at least 5 voxels to average over after reduction above
    ##otherwise return NA time series
    
    ##cat("  ROI ", attr(roivox, "maskval"), ": fewer than 5 voxels had acceptable time series. Removing this ROI from correlations.\n", file=".roilog", append=TRUE)
    message("  ROI ", attr(roivox, "maskval"), ": fewer than 5 voxels had acceptable time series. Removing this ROI from correlations.")
    ts <- rep(NA_real_, nrow(roivox))
  } else {
    if (sum(badvox) > 0) {
      ##cat("  ROI ", attr(roivox, "maskval"), ": ", sum(badvox), " voxels had bad time series (e.g., constant) and were removed prior to ROI averaging.\n", file=".roilog", append=TRUE)
      message("  ROI ", attr(roivox, "maskval"), ": ", sum(badvox), " voxels had bad time series (e.g., constant) and were removed prior to ROI averaging.")
      roivox <- roivox[,!badvox] #remove bad voxels (columns)
    }
    
    if (roi_reduce == "pca") {
      ts <- prcomp(roivox, scale.=TRUE)$x[,1] #first eigenvector
      tsmean <- apply(roivox, 1, mean, na.rm=TRUE)
      #flip sign of component to match observed data (positive correlation)
      if (cor(ts, tsmean) < 0) { ts <- -1*ts }
    } else if (roi_reduce == "mean") {
      ts <- apply(roivox, 1, mean, na.rm=TRUE) #mean time series across voxels
    } else if (roi_reduce == "median") {
      ts <- apply(roivox, 1, median, na.rm=TRUE)
    } else if (roi_reduce == "huber") {
      ts <- apply(roivox, 1, getRobLocation)
    }
  }
  
  return(ts)
}


##drop initial volumes if requested
##this should be implemented before censoring so that ARIMA models and bandpass filtering are applied on the truncated data
##ARIMA estimation may be thrown off if we hand it nonstationary data, such as may exist prior to steady state magnetization
if (drop_vols > 0) {
  message("Dropping ", drop_vols, " volumes from ROI time series prior to correlation.")
  roiavgmat <- roiavgmat[-1*(1:drop_vols),]
  if (!is.null(nuisance_regressors)) { nuisance_df <- nuisance_df[-1*(1:drop_vols),] }
  if (length(censor_vols) > 0L) {
    censor_vols <- censor_vols - drop_vols #shift censor vector based on the number dropped
    censor_vols <- censor_vols[censor_vols > 0] #omit any censored volumes that may have fallen in the truncated period
  }
}

#We need to detrend before bandpass filtering, if requested
#note that we could incorporate these as additional columns of the nuisance regressors if there was no temporal filter
#then apply those regressors at the stage of ARIMA or OLS removal. But then we start
if (!is.null(detrend_order)) {
  message("Applying detrending of order: ", detrend_order, " to ROI aggregated time series.")
  roiavgmat <- apply(roiavgmat, 2, function(ts) { detrendts(ts, order=detrend_order) })
  
  if (!is.null(nuisance_regressors)) {
    message("Applying detrending of order: ", detrend_order, " to nuisance regressors.")
    nuisance_df <- apply(nuisance_df, 2, function(ts) { detrendts(ts, order=detrend_order) })
  }
}

##NB. I have developed significant trepidation about temporal filtering when computing correlations among regions.
##For reasoning, see Bright & Murphy 2017 NeuroImage; Arbabsharani et al., 2014 NeuroImage; Davey et al., 2013, NeuroImage
##Minimally, temporal filtering reduces degrees of freedom and induces autocorrelation (see Fig2a,b in Bright)
##But this becomes even more difficult if we attempt to model the autoregressive structure of the data through prewhitening prior to correlation
##In this case, the ARIMA model itself acts as a temporal filter on the data.
##Thus, if requested, we filter first, then apply ARIMA, but with a scary warning!
if (!is.null(delta_t)) {
  message("Bandpass filtering ROI time series prior to computing correlations using FIR1 filter.")
  message("Sampling frequency: dt = ", delta_t, "s, freq_low = ", freq_low, "Hz, freq_high = ", freq_high, "Hz")
  message("Be cautious about correlations among time series after temporal filtering. See Davey et al., 2013 and Bright and Murphy 2017")
  
  if (!is.null(fit_p)) {
    message("Be *really* cautious about combining temporal filtering and AR(I)MA modeling!")
    message("We will filter first, then fit ARIMA models after, but AR(I)MA is itself a filter.")
    message("See Arbabshirani et al., 2014 and Bright and Murphy 2017 for warnings")
  }
  
  filter_order <- min(length(ts), 300) #don't have more than n filter frequencies (doesn't do anything bad, but it has to be overkill)
  tictoc::tic("Filtering ROI time series")
  #roiavgmat <- apply(roiavgmat, 2, function(ts) {
  roiavgmat <- foreach(ts=iter(roiavgmat, by="column"), .inorder=TRUE, .combine=cbind, .noexport="roiavgmat") %dopar% {
    fir1Bandpass(ts, TR=delta_t, low=freq_low, high=freq_high, n=filter_order, padx=100, detrend=NULL)
  }
  tictoc::toc()
  
  if (!is.null(nuisance_regressors)) {
    message("We will apply the same filter to nuisance regressors prior to their removal.")
    nuisance_df <- apply(nuisance_df, 2, function(ts) {
      fir1Bandpass(ts, TR=delta_t, low=freq_low, high=freq_high, n=filter_order, padx=100, detrend=NULL)
    })
  }
}

##fit ARIMA model to each time series to prewhiten, if requested
##z-score nuisance regressors to help with numerical optimization (we don't care about parameter estimates)
xreg <- if (is.null(nuisance_regressors)) { NULL } else { scale(nuisance_df) }

#TODO: implement regression with spike regressors instead of chopping time points below.
if (!is.null(fit_p)) {
  message("Fitting ARIMA model to each ROI time series prior to computing correlations.")
  message("Model order: p = ", fit_p, ", d = ", fit_d, ", q = ", fit_q)
  message("Note that ARIMA model is fit prior to censoring or truncation given the importance of temporal continguity.")
  message("Mean and linear drift terms are included in ARIMA model by default.")
  
  #NB. The correct approach for estimating cross-correlation on prewhitened time series is to estimate an ARIMA model on one
  #series (x), then apply that model to the other series (y), rather than fitting a new model on the second series (y).
  #http://finzi.psych.upenn.edu/R/library/TSA/html/prewhiten.html
  #https://stats.stackexchange.com/questions/191131/testing-significance-of-cross-correlated-series/
  #https://stats.stackexchange.com/questions/194130/cross-correlation-of-two-autocorrelated-signals-removing-autocorrelation-with-a
  #https://onlinecourses.science.psu.edu/stat510/node/75/
  #http://support.sas.com/documentation/cdl/en/etsug/63348/HTML/default/viewer.htm#etsug_arima_sect033.htm
  
  #but, in the fMRI literature, people often apply a unique ARIMA to each voxel/unit and keep the residuals. When examining 
  #the CCF, this means the input and response series will not have the same filter applied, which is conceptually problematic.
  #Thus, return a set of models here for each ROI, and apply them appropriately during cross-correlation (essentially cell by cell)
  
  tictoc::tic("Fitting ARIMA models to ROI time series")
  #roi_arima_fits <- apply(roiavgmat, 2, function(ts) {
  roi_arima_fits <- foreach(ts=iter(roiavgmat, by="column"), .inorder=TRUE, .noexport="roiavgmat") %dopar% {
    tryCatch(forecast::Arima(as.vector(ts), order=c(fit_p, fit_d, fit_q), include.drift=TRUE, include.mean=TRUE, xreg=xreg, method="CSS-ML"), #fit model of specified order
      error=function(e) { print(e); forecast::Arima(as.vector(ts), order=c(fit_p, fit_d, fit_q), include.drift=TRUE, include.mean=TRUE, xreg=xreg, method="ML") }) #fall back to full ML
  }
  tictoc::toc()
  
  #There is a challenge of how to handle prewhitening in the presence of exogenous covariates (xreg) since simply passing
  #the model argument to Arima will fit the other time series with the exact coefficients for both the ARMA part and the
  #covariates. This seems problematic if the other time series has a substantially different trend or effect of the nuisance signals.
  #In prewhitening the premise is that the ARMA coefficients act as a filter on both Y and X (see Shumway & Stoffer 2017, p. 267),
  #not that the covariates are part of the filter.
  
  #The most principled thing I can think of is that we wish to estimate the structure of the *noise* in one series, which means that the ARMA
  #coefficients must be derived conditioned on the external covariates. This is a regression model with ARMA errors.
  #When we wish to examine the cross-correlation, however, the ARMA coefficients are the *filter* that changes the temporal structure (and ensures)
  #that cross-correlation is not due to the autoregressive properties of the input. BUT, we also wish to derive unique estimates of the
  #response series for intercept, drift, and exogenous covariates. Thus, we should fix the ARMA coefficients to equality between models, but
  #uniquely estimate the external regressor effects in each.
  #See bottom of this page:
  #http://support.sas.com/documentation/cdl/en/etsug/63939/HTML/default/viewer.htm#etsug_arima_sect012.htm
  
  #bg tests the null of no serial correlation *up to* the max order. Thus, no real need to test sequential lags, just the max.
  message("Checking whiteness of ARIMA residuals using Breusch-Godfrey test.")
  message("This tests the null hypothesis of serial correlation up to lag order: ", white_lags, " [hard coded for now]")
  whitechecks <- sapply(roi_arima_fits, function(mout) {
    #deprecate Ljung-Box; use Breusch-Godfrey
    #https://stats.stackexchange.com/questions/148004/testing-for-autocorrelation-ljung-box-versus-breusch-godfrey
    #test <- Box.test(mout$residuals, lag=lag, type="Ljung-Box")
    #test$p.value
    test <- lmtest::bgtest(mout$residuals ~ 1, order=white_lags)
    test$p.value
  })
  message("The percentage of non-white residuals is: ", pct_nonwhite <- round(100*sum(whitechecks < .05)/length(whitechecks), 3), "%. This should be no more than ~5% given a nominal Type I error rate.") 
  if (pct_nonwhite > 5) { warning("Whiteness checks failed in ARIMA modeling. Consider changing (increasing) your model order!!")}
} else if (!is.null(nuisance_regressors)) {
  message("Removing nuisance regressors through OLS regression. See Bright & Murphy 2017 for why this is probably suboptimal!")
  roiavgmat <- apply(roiavgmat, 2, function(ts) {
    residuals(lm(ts ~ xreg)) #remove columns in nuisance
  })
}

#keep intelligent names on time x roi matrix
colnames(roiavgmat) <- paste0("roi", maskvals)
rownames(roiavgmat) <- paste0("vol", 1:nrow(roiavgmat))

#handle roi diagnostics output, adding whiteness p-value checks if relevant
if (!is.null(roi_diagnostics_fname)) {
  tt <- as.data.frame(table(roimask_prebrainmask))
  names(tt) <- c("maskval", "nvox_total")
  tt$maskval <- as.numeric(as.character(tt$maskval))
  tt <- subset(tt, maskval != 0)
  roi_diagnostics <- data.frame(dataset=fname_rsproc, maskval=maskvals, nvox_good=nvox_good_per_roi, nvox_observed=nvox_observed_per_roi)
  roi_diagnostics <- merge(roi_diagnostics, tt, by="maskval")
  roi_diagnostics$prop_masked <- with(roi_diagnostics, 1 - nvox_observed/nvox_total)
  roi_diagnostics$prop_missing <- with(roi_diagnostics, 1 - nvox_good/nvox_observed)
  if (!is.null(roi_arima_fits)) { roi_diagnostics$bg_arima_white_pval <- whitechecks }
  write.csv(file=roi_diagnostics_fname, roi_diagnostics, row.names=FALSE)
}

##apply censoring to resulting time series
censorvec <- rep(0, nrow(roiavgmat))
goodVols <- 1:nrow(roiavgmat)

if (length(censor_vols) > 0L) {
  message("Censoring volumes ", paste0(censor_vols, collapse=", "), " based on ", fname_censor1D)
  goodVols <- goodVols[-censor_vols]
  censorvec[censor_vols] <- 1
}

roiavgmat_censored <- roiavgmat[goodVols,]

##output ts file if requested
if (nchar(ts_out_file) > 0L) {
  message("Writing time series to: ", ts_out_file)
  if(is.null(fname_censor1D)) {
    df <- roiavgmat
  } else {
    df <- cbind(censor=censorvec, roiavgmat)
  }
  write.table(df, file=ts_out_file, col.names=write_header, row.names=FALSE)
  
  #output individually prewhitened time series from ARIMA, if available
  if (!is.null(roi_arima_fits)) {
    pw_ts_file <- paste0(tools::file_path_sans_ext(ts_out_file, compression = TRUE), "_prewhitened", file_ext(ts_out_file))
    message("Writing prewhitened time series to: ", pw_ts_file)
    roiavgmat_whitened <- sapply(roi_arima_fits, residuals)
    colnames(roiavgmat_whitened) <- paste0("roi", maskvals)
    rownames(roiavgmat_whitened) <- paste0("vol", 1:nrow(roiavgmat_whitened))
    if(is.null(fname_censor1D)) {
      df <- roiavgmat_whitened
    } else {
      df <- cbind(censor=censorvec, roiavgmat_whitened)
    }
    write.table(df, file=pw_ts_file, col.names=write_header, row.names=FALSE)
  }
}

#If we only want the timeseries, quit before computing correlations
if (ts_only) {
  if (njobs > 1) { try(stopCluster(clusterobj)) }
  quit(save="no", 0, FALSE)
}

pp <- ifelse(partial==TRUE, "_partial", "") #insert 'partial' into message and file name if relevant
pw <- ifelse(is.null(roi_arima_fits), "", "_prewhitened") #insert 'partial' into message and file name if relevant
for (m in 1:length(corr_method)) {
  message("Computing", sub("_", " ", pp[m]), " correlations among ROI times series using method: ", ifelse(corr_method[m]=="auto", "robust", corr_method[m]))
  if (partial[m] == TRUE && !is.null(roi_arima_fits)) {
    message("At present, cannot combine partial correlation and ARIMA-based prewhitening since there is no single time series matrix. Skipping this combination.")
    next
  }
  cormat <- genCorrMat(as.data.frame(roiavgmat_censored), method=corr_method[m], fisherz=fisherz, partial=partial[m], roi_arima_fits=roi_arima_fits)
  this_out_file <- paste0(tools::file_path_sans_ext(out_file, compression = TRUE), "_", corr_method[m], pp[m], pw, file_ext(out_file)) #add method-specific suffix to file
  message("Writing correlations to: ", this_out_file)
  if (grepl(".*\\.gz$", this_out_file, perl=TRUE)) {
    ##write compressed
    gzf <- gzfile(this_out_file, "w")
    write.table(cormat, file=gzf, col.names=write_header, row.names=FALSE, na=na_string)
    close(gzf)
  } else {
    write.table(cormat, file=this_out_file, col.names=write_header, row.names=FALSE, na=na_string)
  }
}

if (njobs > 1) { try(stopCluster(clusterobj)) }
quit(save="no", 0, FALSE)
