#this script is not dependent on run_model_index because it reads raw data at l1
Sys.setenv(fsl_pipeline_file="/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")
Sys.setenv(run_model_index=1)

to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))
#source(file.path(fsl_model_arguments$pipeline_home, "functions", "deconvolve_funcs.R"))
#source(file.path(fsl_model_arguments$pipeline_home, "functions", "spm_funcs.R"))

library(tidyverse)
library(abind)
library(oro.nifti)
library(reshape2)
library(dependlab)
library(oro.nifti)
library(parallel)
library(foreach)
library(doParallel)
library(readr)

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model
feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model

#used for reading l1 data
#load(file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(paste(fsl_model_arguments$analysis_name, feat_run_outdir, "lvl2_inputs", sep="_"), ".RData")))

load(file.path(fsl_model_arguments$pipeline_home, "configuration_files", "MMClock_aroma_preconvolve_fse_groupfixed_sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed_lvl2_inputs.RData"))

#registerDoSEQ()
cl <- makeCluster(20) #hard code for now
registerDoParallel(cl)
#clusterExport(cl, c("sigmoid", "spm_hrf", "generate_feature", "dsigmoid", "deconvolve_nlreg", "deconvolve_nlreg_resample")) #make sure functions are available

subinfo$dir_found <- file.exists(subinfo$mr_dir)

hrf_pad <- 32 #based on spm_hrf for 1.0 second TR. Used to chop output of deconvolvefilter

#Schaefer 400
#NB: this is not an NN-interpolated file. So, fractional values get dumped by -thr -uthr below
#master <- "/gpfs/group/mnh5174/default/lab_resources/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_7Networks_order_fonov_mni152_2.3mm_ants.nii.gz"

#run once outside of R
#orig <- "/gpfs/group/mnh5174/default/lab_resources/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_7Networks_order_fonov_mni152_1mm_ants.nii.gz"
#orig=/gpfs/group/mnh5174/default/lab_resources/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_7Networks_order_fonov_mni152_1mm_ants.nii.gz
#ResampleImageBySpacing 3 $orig Schaefer2018_400Parcels_7Networks_order_fonov_2.3mm_ants_nn.nii.gz 2.3 2.3 2.3 0 0 1

#maintains RPI coordinates to match NIfTI files (use this)
#flirt -in $orig -ref $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii -out Schaefer2018_400Parcels_7Networks_order_2.3mm_flirt_nn -applyisoxfm 2.3 -interp nearestneighbour

master <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/Schaefer2018_400Parcels_7Networks_order_2.3mm_flirt_nn.nii.gz"
##57 is L primary motor for finger/hand
##18 is L V1
##215 is R V1

#system(paste0("fslmaths ", master, " -thr 57 -uthr 57 -bin /gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/masks/l_motor_2.3mm -odt char"))
#system(paste0("fslmaths ", master, " -thr 256 -uthr 256 -bin /gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/masks/r_motor_2.3mm -odt char"))
#system(paste0("fslmaths ", master, " -thr 18 -uthr 18 -bin /gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/masks/l_v1_2.3mm -odt char"))
#system(paste0("fslmaths ", master, " -thr 215 -uthr 215 -bin /gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/masks/r_v1_2.3mm -odt char"))

hippo_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise"
mask_dir <- file.path(hippo_dir, "masks")

## atlas_files <- c(
##   file.path(mask_dir, "long_axis_l_2.3mm.nii.gz"),
##   file.path(mask_dir, "long_axis_r_2.3mm.nii.gz"),
##   file.path(mask_dir, "l_motor_2.3mm.nii.gz"),
##   file.path(mask_dir, "r_motor_2.3mm.nii.gz"),
##   file.path(mask_dir, "l_v1_2.3mm.nii.gz"),
##   file.path(mask_dir, "r_v1_2.3mm.nii.gz"),
##   file.path(mask_dir, "harvardoxford-subcortical_prob_Left_Accumbens_2009c_thr20_2.3mm.nii.gz"),
##   file.path(mask_dir, "harvardoxford-subcortical_prob_Right_Accumbens_2009c_thr20_2.3mm.nii.gz"),
##   file.path(mask_dir, "hippo_dcm/masks/vmpfc_clust1_z5.7_2009c.nii") #add vmPFC from NeuroSynth
## )

#Just the new cobra-based long axis masks
atlas_files <- c(
  file.path(mask_dir, "long_axis_l_cobra_2.3mm.nii.gz"),
  file.path(mask_dir, "long_axis_r_cobra_2.3mm.nii.gz")
)

atlas_imgs <- lapply(atlas_files, readNIfTI, reorient=FALSE)

#whether to use unsmoothed data
data_type <- "smooth_in_mask" #alternatives "unsmoothed" "smoothed"
bush2015 <- FALSE #whether to use Bush 2015 algorithm or 2011 algorithm

out_dir <- file.path(hippo_dir, "deconvolved_timeseries", data_type)
if (bush2015) { out_dir <- paste0(out_dir, "_b2015") }

#use this for 1st-level results
l1_inputs <- feat_l2_inputs_df$feat_dir

TR <- 1.0 #seconds

#determine nifti files for atlas
l1_niftis <- sapply(l1_inputs, function(x) {
  fsf <- readLines(file.path(x, "design.fsf"))
  nifti <- grep("^set feat_files\\(1\\)", fsf, perl=TRUE, value=TRUE)
  stopifnot(length(nifti)==1L)
  nifti <- paste0(sub("set feat_files\\(1\\) \"([^\"]+)\"", "\\1", nifti, perl=TRUE), ".nii.gz")
  return(nifti)
})


#handle unsmoothed or smooth in mask: create truncated files that match smoothed side
if (data_type != "smoothed") {
  l1_niftis <- sub("mni_5mm_aroma", "mni_nosmooth_aroma", l1_niftis)
  if (data_type == "unsmoothed") {
    l1_niftis <- sub("nfaswuktm", "nfawuktm", l1_niftis)
  } else if (data_type=="smooth_in_mask") {
    #l1_niftis <- sub("nfaswuktm_(.*)\\.nii\\.gz", "fawuktm_\\1_hippblur.nii.gz", l1_niftis) #this is cobra
    l1_niftis <- sub("nfaswuktm_(.*)\\.nii\\.gz", "fawuktm_\\1_hippblur_harvardoxford.nii.gz", l1_niftis) #this is harvard-oxford
  }

  #strip smoothing suffix
  l1_niftis <- sub("(_clock\\d+)_5", "\\1", l1_niftis)

  novalue <- foreach(x=iter(l1_niftis), .packages="RNifti") %dopar% {
    if (!file.exists(x)) {
      orig_file <- sub("_drop\\d+(_trunc\\d+)*", "", x, perl=TRUE)
      dropn <- as.numeric(sub(".*drop(\\d+).*", "\\1", x, perl=TRUE))
      if (grepl("_trunc\\d+", x)) {
        truncn <- as.numeric(sub(".*_trunc(\\d+).*", "\\1", x, perl=TRUE))
      } else {
        truncn <- RNifti::niftiHeader(orig_file)$dim[5] #zero based indexing in fslroi
      }
      system(paste("fslroi", orig_file, x, dropn, truncn - dropn))
    }

    return(x)
  }
}

#l1_niftis should now contain the unsmoothed nifti files of the same length as the smoothed side

if (!all(fexists <- file.exists(l1_niftis))) {
  print(l1_niftis[!fexists])
  stop("Could not locate all l1 input file")
}

#l1_subset <- grep("11302|11305|11228|11366", l1_niftis, perl=TRUE)
#l1_niftis <- l1_niftis[l1_subset]
#feat_l2_inputs_df <- feat_l2_inputs_df[l1_subset,]

#Bush deconvolution settings
nev_lr <- .01 #neural events learning rate (default in algorithm)
epsilon <- .005 #convergence criterion (default)
kernel <- spm_hrf(TR)$hrf #canonical SPM difference of gammas
Nresample <- 25

#for testing
#l1_niftis <- l1_niftis[1:5]

#loop over atlas files
for (ai in 1:length(atlas_files)) {
  cat("Working on atlas: ", atlas_files[ai], "\n")
  aimg <- atlas_imgs[[ai]]
  a_indices <- which(aimg != 0, arr.ind=TRUE)
  a_coordinates <- cbind(a_indices, t(apply(a_indices, 1, function(r) { translateCoordinate(i=r, nim=aimg, verbose=FALSE) })))
  a_coordinates <- as.data.frame(a_coordinates) %>% setNames(c("i", "j", "k", "x", "y", "z")) %>%
    mutate(vnum=1:n(), atlas_value=aimg[a_indices], atlas_name=basename(atlas_files[ai])) %>%
    mutate_at(vars(x,y,z), round, 2) %>% select(vnum, atlas_value, everything())

  atlas_img_name <- basename(sub(".nii(.gz)*", "", atlas_files[ai], perl=TRUE))
  dir.create(file.path(out_dir, atlas_img_name, "deconvolved"), showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(out_dir, atlas_img_name, "original"), showWarnings=FALSE, recursive=TRUE)
  
  #loop over niftis
  ff <- foreach(si = 1:length(l1_niftis), .packages=c("dplyr", "readr", "data.table", "reshape2", "dependlab", "foreach", "iterators"),
    .export=c("l1_niftis", "kernel", "epsilon", "nev_lr")) %dopar% {
    #for (si in 1:length(l1_niftis)) {
    
    this_subj <- feat_l2_inputs_df %>% select(subid, run_num, contingency, emotion, drop_volumes) %>% dplyr::slice(si)
    out_name <- file.path(out_dir, atlas_img_name, "deconvolved", with(this_subj[1,,drop=F], paste0(subid, "_run", run_num, "_", atlas_img_name, "_deconvolved.csv.gz")))

    if (file.exists(out_name)) {
      message("Deconvolved file already exists: ", out_name)
      return(NULL)
    }

    cat("  Deconvolving subject: ", l1_niftis[si], "\n")
    dump_out <- tempfile()
    afnistat <- runAFNICommand(paste0("3dmaskdump -mask ", atlas_files[ai], " -o ", dump_out, " ", l1_niftis[si]))
    ts_out <- data.table::fread(dump_out) #read time series
    #ijk <- ts_out[, 1:3] + 1 #keep in 1-based indexing in outputs (AFNI uses zero-based)
    #names(ijk) <- c("i", "j", "k")
    #ijk$vnum <- 1:nrow(ijk)

    #to_deconvolve is a voxels x time matrix
    to_deconvolve <- as.matrix(ts_out[, -1:-3]) #remove ijk

    to_deconvolve <- t(apply(to_deconvolve, 1, scale)) #need to unit normalize for algorithm not to choke on near-constant 100-normed data
    #to_deconvolve <- t(apply(to_deconvolve, 1, function(x) { scale(x, scale=FALSE) })) #just demean, which will rescale to percent signal change around 0 (this matches Bush 2015)
    #to_deconvolve <- t(apply(to_deconvolve, 1, function(x) { x/mean(x)*100 - 100 })) #pct signal change around 0

    temp_i <- tempfile()
    temp_o <- tempfile()

    #to_deconvolve %>%  as_tibble() %>% write_delim(path=temp_i, col_names=FALSE)

    #zero pad tail end (based on various readings, but not original paper)
    #this was decided on because we see the deconvolved signal dropping to 0.5 for all voxels

    to_deconvolve %>% cbind(matrix(0, nrow=nrow(to_deconvolve), ncol=hrf_pad)) %>% as_tibble() %>% write_delim(path=temp_i, col_names=FALSE)

    #to_deconvolve %>% cbind(matrix(0, nrow=nrow(to_deconvolve), ncol=hrf_pad)) %>% as_tibble() %>%
    #  write_delim(path="/gpfs/group/mnh5174/default/lab_resources/deconvolution-filtering/cpp/bin/break.txt", col_names=FALSE)

    #vv <- to_deconvolve[117,]
    #vv <- vv - min(vv)
    #vv <- vv/max(vv)
    #K <- hrf_pad
    #dummy <- runif(ncol(to_deconvolve)+K)*2e-9 - 1e-9
    #activation_indices <- seq(K+1, K+ncol(to_deconvolve))
    #dummy[activation_indices] <- log(vv/(1-vv))

    #fixed <- as.matrix(read.table("/gpfs/group/mnh5174/default/lab_resources/deconvolution-filtering/cpp/bin/break_decon.txt"))
    #fixed <- fixed[,c(-1*1:hrf_pad, -1*365:396)]

    #test1 <- deconvolve_nlreg(to_deconvolve[117,], kernel=kernel, nev_lr=nev_lr, epsilon=epsilon)
    #test2 <- deconvolve_nlreg(to_deconvolve[118,], kernel=kernel, nev_lr=nev_lr, epsilon=epsilon)

    #cor(test1, fixed[1,])
    #all.equal(test1, fixed[1,])
    #cor(test2, fixed[2,])
    #summary(test1)
    #summary(fixed[1,])

    if (bush2015) {
      #use R implementation of Bush 2015 algorithm
      alg_input <- as.matrix(read.table(temp_i))
      deconv_mat <- foreach(vox_ts=iter(alg_input, by="row"), .combine="rbind", .packages=c("dependlab")) %do% {
        reg <- tryCatch(deconvolve_nlreg_resample(as.vector(vox_ts), kernel=kernel, nev_lr=nev_lr, epsilon=epsilon, Nresample=Nresample),
          error=function(e) { cat("Problem deconvolving: ", l1_niftis[si], as.character(e), "\n", file="deconvolve_errors", append=TRUE); return(rep(NA, length(vox_ts))) })

        if (is.list(reg)) { reg <- reg$NEVmean } #just keep the mean resampled events vector
        return(reg)
      }
    } else {
      #use C++ implementation of Bush 2011 algorithm

      #fo is 1/TR, I think...
      #Looking at spm_hrf, it generates a vector of 33 values for the HRF. deconvolvefilter pads the time series at the beginning by this length
      #if you don't return a convolved result, it doesn't do the trimming for you...
      res <- system(paste0("/gpfs/group/mnh5174/default/lab_resources/bin/deconvolvefilter -i=", temp_i, " -o=", temp_o, " -convolved=0 -fo=1 -thread=2"), intern=FALSE)
      if (res != 0) {
        cat("Problem deconvolving: ", l1_niftis[si], "\n", file="deconvolve_errors_compiled", append=TRUE)
        deconv_mat <- matrix(NA, nrow=nrow(to_deconvolve), ncol=ncol(to_deconvolve))
      } else {
        #readr may be faster, but it fails on parsing some outputs...
        #deconv_mat <- as.matrix(read_table(temp_o, col_names=FALSE)) %>% unname() #remove names to avoid confusion in melt
        #deconv_mat <- scan(temp_o, what="character", sep="\n") %>% strsplit("\\s+")
        #if (length(unique(sapply(deconv_mat, length))) != 1L) { browser() } 

        deconv_mat <- as.matrix(read.table(temp_o, header=FALSE)) %>% unname()  #remove names to avoid confusion in melt

        #NB. 17Apr2019. I modified the compiled C++ program to chop the leading zeros itself for all outputs (rather than leaving the leading hrf_pad)
        #deconv_mat <- deconv_mat[,c(-1*1:hrf_pad, seq(-ncol(deconv_mat), -ncol(deconv_mat)+hrf_pad-1))] #trim leading and trailing padding

      }

    }

    #trip hrf end-padding for both C++ and R variants
    deconv_mat <- deconv_mat[,c(seq(-ncol(deconv_mat), -ncol(deconv_mat)+hrf_pad-1))] #trim trailing padding added above

    #to_deconvolve_augment <- cbind(matrix(0, nrow=nrow(to_deconvolve), ncol=hrf_pad), to_deconvolve)
    #deconv_mat <- foreach(vox_ts=iter(to_deconvolve_augment[1:10,], by="row"), .combine="rbind") %dopar% {
    #  reg <- tryCatch(deconvolve_nlreg(as.vector(vox_ts), kernel=kernel, nev_lr=nev_lr, epsilon=epsilon),
    #    error=function(e) { cat("Problem deconvolving: ", l1_niftis[si], "\n", file="deconvolve_errors", append=TRUE); return(rep(NA, length(vox_ts))) })
    #  return(reg)
    #}

    #14Apr2019: Yes, looking at the results from the local function versus the result of the compiled code,
    # there is a 32-column padding on the front of the results from the deconvolvefilter output.
    # Thus, we need to chop this off posthoc to get the alignment to be proper
    #all_ccf <- lapply(1:10, function(ii) {
    #  vv <- ccf(deconv_mat_loc[ii,], deconv_mat[ii,], plot=FALSE)
    #  vv$acf[vv$lag==0]
    #})
    
    #deconv_mat <- t(parApply(cl=cl, X=to_deconvolve, MARGIN=1, FUN=function(vox_ts) { #transpose to maintain voxels x time
    #  deconvolve_nlreg(vox_ts, kernel=kernel, nev_lr=nev_lr, epsilon=epsilon)
    #  #res <- tryCatch(deconvolve_nlreg(vox_ts, kernel=kernel, nev_lr=nev_lr, epsilon=epsilon),
    #  #  error=function(e) { browser(); return(rep(NA, length(ts))) }) #return NA in event of failure
    #  #return(res)
    #}))

    ## deconv_mat <- t(apply(X=to_deconvolve, 1, FUN=function(ts) { #transpose to maintain voxels x time
    ##   browser()
    ##   tryCatch(deconvolve_nlreg(ts, kernel=kernel, nev_lr=nev_lr, epsilon=epsilon),
    ##     error=function(e) { print(e); return(rep(NA, length(ts))) }) #return NA in event of failure
    ## }))

    #melt this for combination
    deconv_melt <- reshape2::melt(deconv_mat, value.name="decon", varnames=c("vnum", "time"))
    to_deconvolve_melt <- reshape2::melt(to_deconvolve, value.name="BOLD_z", varnames=c("vnum", "time"))
    
    deconv_df <- deconv_melt %>% mutate(vnum=as.numeric(vnum)) %>% left_join(a_coordinates, by="vnum") %>%
      mutate(nifti=sub("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", "", l1_niftis[si], fixed=TRUE)) %>%
      cbind(feat_l2_inputs_df %>% select(subid, run_num, contingency, emotion, drop_volumes) %>% dplyr::slice(si)) %>% #add metadata
      mutate(time=time+drop_volumes-1) %>% select(-drop_volumes, -nifti) %>% #dropping nifti for now to save on file size
      select(subid, run_num, contingency, emotion, time, atlas_name, atlas_value, vnum, x, y, z, decon) #omitting i,j,k for now

    orig_df <- to_deconvolve_melt %>% mutate(vnum=as.numeric(vnum)) %>% left_join(a_coordinates, by="vnum") %>%
      mutate(nifti=sub("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", "", l1_niftis[si], fixed=TRUE)) %>%
      cbind(feat_l2_inputs_df %>% select(subid, run_num, contingency, emotion, drop_volumes) %>% dplyr::slice(si)) %>% #add metadata
      mutate(time=time+drop_volumes-1) %>% select(-drop_volumes, -nifti) %>% #dropping nifti for now to save on file size
      select(subid, run_num, contingency, emotion, time, atlas_name, atlas_value, vnum, x, y, z, BOLD_z)

    out_name <- file.path(out_dir, atlas_img_name, "deconvolved", with(this_subj[1,,drop=F], paste0(subid, "_run", run_num, "_", atlas_img_name, "_deconvolved.csv.gz")))
    write_csv(deconv_df, path=out_name)
    
    out_name <- file.path(out_dir, atlas_img_name, "original", with(this_subj[1,,drop=F], paste0(subid, "_run", run_num, "_", atlas_img_name, "_original.csv.gz")))
    write_csv(orig_df, path=out_name)
    
  }
}
