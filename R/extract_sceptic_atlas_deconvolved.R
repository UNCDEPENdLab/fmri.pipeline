#this script is not dependent on run_model_index because it reads raw data at l1
Sys.setenv(fsl_pipeline_file="/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")
Sys.setenv(run_model_index=1)

to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

source("voxelwise_deconvolution.R")
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

#####
#vestiges of things that should be passed in

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model
feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model

#used for reading l1 data
#load(file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(paste(fsl_model_arguments$analysis_name, feat_run_outdir, "lvl2_inputs", sep="_"), ".RData")))

load(file.path(fsl_model_arguments$pipeline_home, "configuration_files", "MMClock_aroma_preconvolve_fse_groupfixed_sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed_lvl2_inputs.RData"))

#registerDoSEQ()
#cl <- makeCluster(2) #hard code for now
#registerDoParallel(cl)
#clusterExport(cl, c("sigmoid", "spm_hrf", "generate_feature", "dsigmoid", "deconvolve_nlreg", "deconvolve_nlreg_resample")) #make sure functions are available

subinfo$dir_found <- file.exists(subinfo$mr_dir)

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

hippo_dir <- "/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise"
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
## atlas_files <- c(
##   file.path(mask_dir, "long_axis_l_cobra_2.3mm.nii.gz"),
##   file.path(mask_dir, "long_axis_r_cobra_2.3mm.nii.gz")
## )

atlas_files <- Sys.getenv("atlas")
if (atlas_files=="") {
  stop("Did not receive atlas as environment variable")
}

## atlas_files <- c("/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_136_2.3mm.nii.gz",
##                  "/proj/mnhallqlab/projects/clock_analysis/fmri/vmpfc_ah/vmpfc_ba_max_2.3mm_filled.nii.gz",
##                  "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/pauli_da_midbrain_2.3mm.nii.gz",
##                  "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/bilateral_striatum_tight_7Networks_2.3mm.nii.gz",
##                  "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/pauli_combined_integermask_2.3mm.nii.gz")

#whether to use unsmoothed data
data_type <- "smoothed" #alternatives "unsmoothed" "smoothed"
bush2015 <- FALSE #whether to use Bush 2015 algorithm or 2011 algorithm

out_dir <- file.path(hippo_dir, "deconvolved_timeseries", data_type)
if (bush2015) { out_dir <- paste0(out_dir, "_b2015") }

#use this for 1st-level results
l1_inputs <- feat_l2_inputs_df$feat_dir

TR <- 1.0 #seconds

l1_inputs <- sub("/gpfs/group/mnh5174/default", "/proj/mnhallqlab/studies", l1_inputs)

#determine nifti files for atlas
l1_niftis <- sapply(l1_inputs, function(x) {
  fsf <- readLines(file.path(x, "design.fsf"))
  nifti <- grep("^set feat_files\\(1\\)", fsf, perl=TRUE, value=TRUE)
  stopifnot(length(nifti)==1L)
  nifti <- paste0(sub("set feat_files\\(1\\) \"([^\"]+)\"", "\\1", nifti, perl=TRUE), ".nii.gz")
  return(nifti)
})

l1_niftis <- sub("/gpfs/group/mnh5174/default", "/proj/mnhallqlab/studies", l1_niftis)


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

  novalue <- foreach(x=iter(l1_niftis), .packages="RNifti") %do% {
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

#for testing
#l1_niftis <- l1_niftis[1:5]

#call voxelwise_deconvolution here
metadata <- feat_l2_inputs_df %>% dplyr::select(subid, run_num, contingency, emotion)

voxelwise_deconvolution(l1_niftis, metadata, out_dir="/proj/mnhallqlab/users/michael/sceptic_decon", TR=1, time_offset=2.0, atlas_files=atlas_files, mask=NULL, nprocs=18, save_original_ts=FALSE,
  out_file_expression=expression(paste0("sub", this_subj$subid, "_run", this_subj$run_num, "_", atlas_img_name)))
