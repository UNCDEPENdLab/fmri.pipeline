# This script saves the cluster beta estimates for SCEPTIC fMRI based on an a priori ROI or atlas

#used for testing: group fixed entropy
#Sys.setenv(fsl_pipeline_file="/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")
#Sys.setenv(run_model_index=1)


#load the master configuration file

to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

extract_z <- Sys.getenv("extract_z") #whether to grab copes or zstats
if (nchar(extract_z) == 0L) {
  extract_z <- FALSE #default to betas
} else {
  if (extract_z %in% c("1", "TRUE", "true")) {
    extract_z <- TRUE
  } else {
    extract_z <- FALSE
  }
} 

load(to_run)

#load("/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData") #current arguments
#run_model_index <- 1

source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))

#1) load spatial maps RData object
#2) rebuild into 4d cube (where fourth dimension is run/subject)
#3) clusterize each effect of interest in stats outputs using 3dclust and generating mask
#library(ggplot2)
library(tidyverse)
library(abind)
library(oro.nifti)
library(reshape2)
library(robust)
library(car)
library(dependlab)
library(oro.nifti)
library(parallel)
library(foreach)
library(doParallel)

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model

feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model
n_l1_copes <- fsl_model_arguments$n_l1_copes[run_model_index] #number of l1 copes determines number of FEAT LVL3 analyses to run (1 per LVL1 cope)
l1_cope_names <- fsl_model_arguments$l1_cope_names[[run_model_index]] #names of l1 copes (used for folder naming)

load(file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(paste(fsl_model_arguments$analysis_name, feat_run_outdir, "lvl2_inputs", sep="_"), ".RData")))

#registerDoSEQ()
cl <- makeCluster(fsl_model_arguments$n_cluster_beta_cpus)
registerDoParallel(cl)

subinfo$dir_found <- file.exists(subinfo$mr_dir)

feat_lvl2_dirname <- "FEAT_LVL2_runtrend.gfeat" #should populate this to the structure at some point

#whether to extract from beta series (very slow)
calculate_beta_series <- FALSE

atlas_files <- c("/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/masks/long_axis_l_2.3mm.nii.gz",
  "/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/masks/long_axis_r_2.3mm.nii.gz")

atlas_imgs <- lapply(atlas_files, readNIfTI, reorient=FALSE)

# For atlas extraction, we have 1) lvl1 inputs, and 2) lvl2 inputs
#   LVL1 will give us run-level decision signals irrespective of condition (those these could be merged in later)
#   LVL2 will give us subject-level contrasts such as scram > fear
#
# We need, however, to identify the inputs to the analyses in both cases.
# Previously, we used the l3 analysis to identify subject inputs.
# Here, we need to go back to original inputs prior to even lvl2 analysis

#use this for 1st-level results
l1_inputs <- feat_l2_inputs_df$feat_dir

l1_niftis <- sapply(l1_inputs, function(x) {
  fsf <- readLines(file.path(x, "design.fsf"))
  nifti <- grep("^set feat_files\\(1\\)", fsf, perl=TRUE, value=TRUE)
  stopifnot(length(nifti)==1L)
  nifti <- paste0(sub("set feat_files\\(1\\) \"([^\"]+)\"", "\\1", nifti, perl=TRUE), ".nii.gz")
  return(nifti)
})

##rework using subinfo structure as the authoritative guide (rather than repeated searches)
copedf <- c()
for (s in 1:nrow(subinfo)) {
  for (cope in 1:n_l1_copes) {
    expectdir <- file.path(subinfo[s,"mr_dir"], fsl_model_arguments$expectdir, feat_run_outdir, feat_lvl2_dirname, paste0("cope", cope, ".feat"))
    if (dir.exists(expectdir)) {
      copedf <- rbind(copedf, data.frame(id=subinfo[s,"id"], model=feat_run_outdir, cope=cope, fsldir=expectdir, stringsAsFactors=FALSE))
    } else {
      message("could not find expected directory: ", expectdir)
    }
  }
}

mdf <- merge(subinfo, copedf, by="id", all.y=TRUE)

#remove bad ids
mdf <- mdf %>% filter(!id %in% fsl_model_arguments$badids)
mdf <- arrange(mdf, id, model, cope)

#this differs from l3 approach where we use the inputs to FLAME

for (l1 in 1:n_l1_copes) {
  l1_contrast_name <- l1_cope_names[l1]
  model_output_dir <- file.path(feat_lvl3_outdir, l1_contrast_name)

  l1_subinfo <- mdf %>% filter(cope==l1) %>% mutate(numid=1:n())
  subject_inputs <- l1_subinfo$fsldir
  
  #generate separate files for each l1 contrast (reset here)
  all_rois <- list()
  all_beta_series <- list()
  
  #figure out what the l2 contrasts are and read relevant statistics for each
  l2fsf <- file.path(dirname(subject_inputs[1]), "design.fsf")
  stopifnot(file.exists(l2fsf))
  l2_syntax <- readLines(l2fsf)

  l2_contrast_info <- grep("\\s*set fmri\\(conname_real\\.\\d+\\).*", l2_syntax, value=TRUE, perl=TRUE)
  l2_contrast_nums <- as.numeric(sub("\\s*set fmri\\(conname_real\\.(\\d+)\\).*", "\\1", l2_contrast_info, perl=TRUE))
  l2_contrast_names <- sub("\\s*set fmri\\(conname_real\\.\\d+\\)\\s*\"?([^\"]+)\"?.*", "\\1", l2_contrast_info, perl=TRUE)

  n_l2_contrasts <- max(l2_contrast_nums)
  l2_contrast_names <- l2_contrast_names[order(l2_contrast_nums)] #order l2 contrast names in ascending order to match l2 loop below

  #hard coding location of beta series analysis for now
  beta_series_inputs <- sub(paste0("^(", fsl_model_arguments$fmri_dir, "/", fsl_model_arguments$idregex, "/", fsl_model_arguments$expectdir,
    ")/.*"), "\\1/sceptic-clock_bs-feedback-preconvolve_fse_groupfixed", subject_inputs)
  
  #loop over l2 contrasts
  #l2_loop_outputs <- list()
  #for (l2 in 1:n_l2_contrasts) {

  #parallel version
  l2_loop_outputs <- foreach(l2=iter(1:n_l2_contrasts), .packages=c("oro.nifti", "dplyr")) %dopar% {
    l2_loop_rois <- list()
    l2_loop_bs <- list()
    
    l2_contrast_name <- l2_contrast_names[l2] #current l2 contrast
    if (extract_z) {
      copefiles <- file.path(subject_inputs, "stats", paste0("zstat", l2, ".nii.gz"))
    } else {
      copefiles <- file.path(subject_inputs, "stats", paste0("cope", l2, ".nii.gz"))
    }
    
    imgdims <- dim(oro.nifti::readNIfTI(copefiles[1], read_data=FALSE))
    
    #generate concatenated cope file image of l2 images (one per subject)
    copeconcat <- array(0, dim=c(imgdims, length(copefiles)))
    #copefiles <- copefiles[file.exists(copefiles)] #don't suggest enabling this in general because it can lead to a mismatch with id
    for (i in 1:length(copefiles)) { copeconcat[,,,i] <- readNIfTI(copefiles[i], reorient=FALSE)@.Data }
    
    atlas_df <- list()

    #should not need l3 in case of atlas
    for (ai in 1:length(atlas_imgs)) {
      aimg <- atlas_imgs[[ai]]
      a_indices <- which(aimg != 0, arr.ind=TRUE)
      a_coordinates <- cbind(a_indices, t(apply(a_indices, 1, function(r) { translateCoordinate(i=r, nim=aimg, verbose=FALSE) })))
      a_coordinates <- as.data.frame(a_coordinates) %>% setNames(c("i", "j", "k", "x", "y", "z")) %>%
        mutate(vnum=1:n(), atlas_value=aimg[a_indices], atlas_name=basename(atlas_files[ai])) %>% select(vnum, atlas_value, everything())
      
      #get betas for each voxel
      beta_mat <- apply(a_coordinates[,c("i", "j", "k")], 1, function(r) { copeconcat[r["i"], r["j"], r["k"], ] })
      
      #reshape into data.frame with beta, numeric numid, and vnum
      if (extract_z) {
        beta_df <- reshape2::melt(beta_mat, varnames=c("numid", "vnum"), value.name="zstat")
      } else {
        beta_df <- reshape2::melt(beta_mat, varnames=c("numid", "vnum"), value.name="beta")
      }
      beta_df <- beta_df %>% inner_join(a_coordinates, by="vnum") %>% select(-i, -j, -k) %>% inner_join(l1_subinfo %>% select(numid, id, fsldir), by="numid")
      
      atlas_df[[ai]] <- beta_df
    }

    atlas_df <- do.call(rbind, atlas_df) %>% mutate(l1_contrast=l1_contrast_name, l2_contrast=l2_contrast_name)       

    #NOT IMPLEMENTED YET
    #handle beta series extraction (NB. beta_series_inputs should be in same order as subject_inputs based on use of sub above)
    if (calculate_beta_series) {          
      beta_series_df <- get_beta_series(beta_series_inputs, roimask, n_bs=50)

      #for identification, add cluster information to beta series from ROI data.frame
      beta_series_df <- beta_series_df %>% left_join(roi_df, by=c("feat_input_id")) %>%
        select(feat_input_id, run, trial, everything())

      #beta_series_df %>% group_by(feat_input_id) %>% summarize(mean(bs_value), mean(cope_value))
    } else {
      beta_series_df <- data.frame()
    }
    
    l2_loop_rois[[paste(l1, l2, sep=".")]] <- atlas_df
    l2_loop_bs[[paste(l1, l2, sep=".")]] <- beta_series_df

    #l2_loop_outputs[[l2]] <- list(rois=l2_loop_rois, beta_series=l2_loop_bs) #serial version
    return(list(rois=l2_loop_rois, beta_series=l2_loop_bs))
  }

  #tack on roi betas from this l2 contrast to the broader set
  all_rois <- bind_rows(rlang::flatten(lapply(l2_loop_outputs, "[[", "rois")))
  all_beta_series <- bind_rows(rlang::flatten(lapply(l2_loop_outputs, "[[", "beta_series")))
  
  #organize models intelligently
  all_rois <- all_rois %>% arrange(l1_contrast, l2_contrast)

  #readr::write_csv(x=all_rois, file.path(model_output_dir, paste0(l1_contrast_name, "_atlas_betas.csv.gz")))
  if (extract_z) {
    readr::write_csv(x=all_rois, file.path(model_output_dir, paste0(l1_contrast_name, "_atlas_zstats.csv.gz")))
  } else {
    readr::write_csv(x=all_rois, file.path(model_output_dir, paste0(l1_contrast_name, "_atlas_betas.csv.gz")))
  }

  if (calculate_beta_series) {
    all_beta_series <- all_beta_series %>% arrange(l1_contrast, l2_contrast, run, trial)
    readr::write_csv(x=all_beta_series, file.path(model_output_dir, paste0(l1_contrast_name, "_atlas_beta_series.csv.gz")))
  }

}


#not uniquely useful at present (CSVs have it all)
#save(all_rois, dmat, file=file.path(model_output_dir, "sceptic_clusters.RData"))

try(stopCluster(cl)) #cleanup pool upon exit of this script
