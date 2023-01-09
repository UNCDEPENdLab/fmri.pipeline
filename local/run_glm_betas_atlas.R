
#setwd("/proj/mnhallqlab/users/michael/fmri.pipeline/example_files")
#load("MMClock_aroma_preconvolve_fse_groupfixed_sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed_lvl2_inputs.RData")
#load("MMClock_aroma_preconvolve_fse_groupfixed.RData")
setwd("/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/configuration_files")
load("MMClock_aroma_preconvolve_fse_groupfixed.RData")

## new_root <- "/proj/mnhallqlab/studies"
## old_root <- "/gpfs/group/mnh5174/default"

## feat_l2_inputs_df$feat_dir <- sub(old_root, new_root, feat_l2_inputs_df$feat_dir)

## fsl_model_arguments$pipeline_home <- sub(old_root, new_root, fsl_model_arguments$pipeline_home)
## fsl_model_arguments$group_output_dir <- sub(old_root, new_root, fsl_model_arguments$group_output_dir)
## fsl_model_arguments$subject_data$mr_dir <- sub(old_root, new_root, fsl_model_arguments$subject_data$mr_dir)
## fsl_model_arguments$fmri_dir <- sub(old_root, new_root, fsl_model_arguments$fmri_dir, fixed=TRUE)

mean_nz <- function(x) {
    nzvals <- abs(x - 0) > .Machine$double.eps*2
    x <- mean(x[nzvals], na.rm=TRUE)
    return(x)
}

atlas_files <- c("/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_136_2.3mm.nii.gz",
                 "/proj/mnhallqlab/projects/clock_analysis/fmri/vmpfc_ah/vmpfc_ba_max_2.3mm_filled.nii.gz",
                 "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/pauli_da_midbrain_2.3mm.nii.gz",
                 "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/bilateral_striatum_tight_7Networks_2.3mm.nii.gz",
                 "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/pauli_combined_integermask_2.3mm.nii.gz")

#pe max z > 4.41 clusters (break up posterior mega-cluster)
atlas_files <- "/proj/mnhallqlab/studies/MMClock/group_analyses/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_z4p417_k50_clusters.nii.gz"


#u_chosen_quantile
#extract_glm_betas_atlas(fsl_model_arguments, atlas_files, run_model_index=1, extract_z=FALSE,
#  extract_beta_series=FALSE, ncpus=4, aggregate=TRUE, aggFUN=mean_nz)

#v_chosen
#extract_glm_betas_atlas(fsl_model_arguments, atlas_files, run_model_index=2, extract_z=FALSE,
#  extract_beta_series=FALSE, ncpus=4, aggregate=TRUE, aggFUN=mean_nz)

#v_entropy
#extract_glm_betas_atlas(fsl_model_arguments, atlas_files, run_model_index=3, extract_z=FALSE,
#  extract_beta_series=FALSE, ncpus=4, aggregate=TRUE, aggFUN=mean_nz)

#rt_vmax_change
#extract_glm_betas_atlas(fsl_model_arguments, atlas_files, run_model_index=4, extract_z=FALSE,
#  extract_beta_series=FALSE, ncpus=4, aggregate=TRUE, aggFUN=mean_nz)

#pe_max at present
extract_glm_betas_atlas(fsl_model_arguments, atlas_files, run_model_index=5, extract_z=FALSE,
  extract_beta_series=FALSE, ncpus=4, aggregate=TRUE, aggFUN=mean_nz)



# #not uniquely useful at present (CSVs have it all)
# #save(all_rois, dmat, file=file.path(model_output_dir, "sceptic_clusters.RData"))
# 

