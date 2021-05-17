# This script saves the cluster beta estimates for SCEPTIC fMRI based on an a priori ROI or atlas

#used for testing: group fixed entropy
#Sys.setenv(fsl_pipeline_file="/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")
#Sys.setenv(run_model_index=1)

#' Function to extract subject-level coefficients from an FSL group analysis
#' 
#' @param atlas_files vector of filenames for NIfTI-based atlases. The function will loop over each and extract coefficients
#' @param extract_z extract lower-level z-statistics if TRUE, or copes (betas) if FALSE
#' @param extract_beta_series whether to extract coefficients from parallel beta series analyses that have already been run (very slow)
#' @param aggregate whether to take the average (or other central tendency measure) of voxels within a given mask value
#' @importFrom checkmate test_integerish
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @export
extract_glm_betas_atlas <- function(
  config_file, atlas_files, run_model_index=1, extract_z=FALSE,
  extract_beta_series=FALSE, ncpus=1, aggregate=TRUE, aggFUN=mean) {

  #TODO: compare NIfTI header of atlas against maps from which you're extracting

  library(tidyverse)
  library(oro.nifti)
  library(reshape2)
  #library(robust)
  #library(car)
  #library(dependlab)
  library(oro.nifti)
  library(parallel)
  library(foreach)
  library(doParallel)

  if (is.character(config_file)) {
    if (!file.exists(config_file)) { stop("Cannot locate configuration file", config_file) }
    load(config_file)
  } else {
    stopifnot(inherits(config_file, "list")) #should eventually be more specific
  }

  checkmate::assert_list(config_file)
  checkmate::assert_data_frame(config_file$subject_data)
  checkmate::assert_numeric(config_file$n_l1_copes)
  checkmate::assert_list(config_file$l1_cope_names)
  stopifnot(length(config_file$l1_cope_names) == length(config_file$n_l1_copes))
  sapply(atlas_files, checkmate::check_file_exists)

  subinfo <- fsl_model_arguments$subject_data
  feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model

  feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model
  n_l1_copes <- fsl_model_arguments$n_l1_copes[run_model_index] #number of l1 copes determines number of FEAT LVL3 analyses to run (1 per LVL1 cope)
  l1_cope_names <- fsl_model_arguments$l1_cope_names[[run_model_index]] #names of l1 copes (used for folder naming)

  load(file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(paste(fsl_model_arguments$analysis_name, feat_run_outdir, "lvl2_inputs", sep="_"), ".RData")))

  ## tmp fix
  new_root <- "/proj/mnhallqlab/studies"
  old_root <- "/gpfs/group/mnh5174/default"

  feat_l2_inputs_df$feat_dir <- sub(old_root, new_root, feat_l2_inputs_df$feat_dir)

  #overriding argument for now
  #ncpus <- fsl_model_arguments$n_cluster_beta_cpus

  if (ncpus > 1L) {
    cl <- makeCluster(ncpus, type="FORK")
    registerDoParallel(cl)
    on.exit(try(stopCluster(cl))) #cleanup pool upon exit of this function
  } else {
    registerDoSEQ()
  }

  #verify that subject-level analysis directories exist
  subinfo$dir_found <- file.exists(subinfo$mr_dir)

  feat_lvl2_dirname <- "FEAT_LVL2_runtrend.gfeat" #should populate this to the structure at some point

  #use this for 1st-level results
  l1_inputs <- feat_l2_inputs_df$feat_dir

  #read NIfTI inputs from .feat directories
  #not needed for basic extraction
  #l1_niftis <- read_feat_inputs(l1_inputs)

  ##rework using subinfo structure as the authoritative guide (rather than repeated searches) for crawling over expected cope directories
  copedf <- c()
  for (s in seq_len(nrow(subinfo))) {
    for (cope in 1:n_l1_copes) {
      expectdir <- file.path(subinfo[s, "mr_dir"], fsl_model_arguments$expectdir, feat_run_outdir, feat_lvl2_dirname, paste0("cope", cope, ".feat"))
      if (dir.exists(expectdir)) {
        copedf <- rbind(copedf, data.frame(id=subinfo[s, "id"], model=feat_run_outdir, cope=cope, fsldir=expectdir, stringsAsFactors=FALSE))
      } else {
        message("could not find expected directory: ", expectdir)
      }
    }
  }

  mdf <- merge(subinfo, copedf, by="id", all.y=TRUE)

  #remove bad ids
  mdf <- mdf %>% filter(!id %in% fsl_model_arguments$bad_ids)
  mdf <- arrange(mdf, id, model, cope)

  atlas_imgs <- lapply(atlas_files, readNIfTI, reorient=FALSE)

  ### Core loop over levels
  for (l1 in 1:n_l1_copes) {
    l1_contrast_name <- l1_cope_names[l1]
    model_output_dir <- file.path(feat_lvl3_outdir, l1_contrast_name)

    l1_subinfo <- mdf %>% dplyr::filter(cope==l1) %>% dplyr::mutate(numid=1:n())
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
    beta_series_inputs <- sub(
      paste0(
        "^(", fsl_model_arguments$fmri_dir, "/", fsl_model_arguments$idregex, "/", fsl_model_arguments$expectdir,
        ")/.*"
      ), "\\1/sceptic-clock_bs-feedback-preconvolve_fse_groupfixed", subject_inputs
    )

    #loop over l2 contrasts
    l2_loop_outputs <- list()
    for (l2 in 1:n_l2_contrasts) {

    ##parallel version
    #l2_loop_outputs <- foreach(l2=iter(1:n_l2_contrasts), .packages=c("oro.nifti", "dplyr")) %do% {
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
      for (i in seq_len(copefiles)) {
        copeconcat[, , , i] <- readNIfTI(copefiles[i], reorient = FALSE)@.Data
      }

      atlas_df <- list()

      #should not need l3 in case of atlas
      for (ai in seq_len(atlas_imgs)) {
        aimg <- atlas_imgs[[ai]]
        atlas_name <- sub("(\\.nii)*(\\.gz)*$", "", basename(atlas_files[ai]))
        a_indices <- which(aimg != 0, arr.ind=TRUE)
        m_vals <- unique(aimg[a_indices])
        is_int_mask <- all(sapply(m_vals, checkmate::test_integerish))
        a_coordinates <- cbind(a_indices, t(apply(a_indices, 1, function(r) { translateCoordinate(i=r, nim=aimg, verbose=FALSE) })))
        a_coordinates <- as.data.frame(a_coordinates) %>%
          setNames(c("i", "j", "k", "x", "y", "z")) %>%
          mutate(vnum = 1:n(), atlas_value = aimg[a_indices], atlas_name = basename(atlas_files[ai])) %>%
          select(vnum, atlas_value, everything())

        #get coefficients for each voxel in the mask for each subject
        coef_mat <- apply(a_coordinates[, c("i", "j", "k")], 1, function(r) { copeconcat[r["i"], r["j"], r["k"], ] })

        #reshape into data.frame with beta, numeric numid, and vnum
        dv_name <- ifelse(isTRUE(extract_z), "zstat", "beta")
        coef_df <- reshape2::melt(coef_mat, varnames = c("numid", "vnum"), value.name = "value") %>%
          inner_join(a_coordinates, by = "vnum") %>%
          select(-i, -j, -k)


        if (is_int_mask && aggregate) {
          coef_df <- coef_df %>%
            group_by(numid, atlas_value) %>%
            summarize(atlas_name = atlas_name[1], value = aggFUN(value)) %>%
            ungroup()
        }

        coef_df <- coef_df %>%
          inner_join(l1_subinfo %>%
          dplyr::select(numid, id, fsldir), by="numid") %>%
          dplyr::rename(!!dv_name := value)

        atlas_df[[ai]] <- coef_df
      }

      atlas_df <- do.call(rbind, atlas_df) %>%
        mutate(l1_contrast = l1_contrast_name, l2_contrast = l2_contrast_name)

      #NOT IMPLEMENTED YET
      #handle beta series extraction (NB. beta_series_inputs should be in same order as subject_inputs based on use of sub above)
      if (extract_beta_series) {
        beta_series_df <- get_beta_series(beta_series_inputs, roimask, n_bs=50)

        #for identification, add cluster information to beta series from ROI data.frame
        beta_series_df <- beta_series_df %>%
          left_join(roi_df, by = c("feat_input_id")) %>%
          select(feat_input_id, run, trial, everything())

        #beta_series_df %>% group_by(feat_input_id) %>% summarize(mean(bs_value), mean(cope_value))
      } else {
        beta_series_df <- data.frame()
      }

      l2_loop_rois[[paste(l1, l2, sep=".")]] <- atlas_df
      l2_loop_bs[[paste(l1, l2, sep=".")]] <- beta_series_df

      l2_loop_outputs[[l2]] <- list(rois=l2_loop_rois, beta_series=l2_loop_bs) #serial version
      #return(list(rois=l2_loop_rois, beta_series=l2_loop_bs))
    }

    #tack on roi betas from this l2 contrast to the broader set
    all_rois <- bind_rows(rlang::flatten(lapply(l2_loop_outputs, "[[", "rois")))
    all_beta_series <- bind_rows(rlang::flatten(lapply(l2_loop_outputs, "[[", "beta_series")))

    #organize models intelligently
    all_rois <- all_rois %>% arrange(l1_contrast, l2_contrast)

    #fsuffix <- ifelse(isTRUE(extract_z), paste0("_", atlas_name, "_zstats.csv.gz"), paste0("_", atlas_name, "_betas.csv.gz"))
    fsuffix <- ifelse(isTRUE(extract_z), "_atlas_zstats.csv.gz", "_betas.csv.gz")
    message("writing: ", file.path(model_output_dir, paste0(l1_contrast_name, fsuffix)))
    readr::write_csv(x=all_rois, file.path(model_output_dir, paste0(l1_contrast_name, fsuffix)))

    if (extract_beta_series) {
      all_beta_series <- all_beta_series %>% arrange(l1_contrast, l2_contrast, run, trial)
      readr::write_csv(x=all_beta_series, file.path(model_output_dir, paste0(l1_contrast_name, "_", atlas_name, "_beta_series.csv.gz")))
    }
  }

}
