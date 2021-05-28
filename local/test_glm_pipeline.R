library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(lgr)
library(foreach)
library(doParallel)
library(dependlab)
library(emmeans)

setwd("/proj/mnhallqlab/users/michael/fmri.pipeline/R")
source("setup_glm_pipeline.R")
source("finalize_pipeline_configuration.R")
source("glm_helper_functions.R")
source("fsl_l1_model.R")
source("setup_l1_models.R")
source("specify_contrasts.R")
source("build_l2_models.R")
source("glm_helper_functions.R")
source("lookup_nifti_inputs.R")
#source("build_design_matrix.R")

trial_df <- readRDS("/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/mmy3_trial_df_selective_groupfixed.rds") %>%
  mutate(rt_sec=rt_csv/1000) %>%
  select(-isi_onset, -iti_onset)

#l1_models <- build_l1_model(trial_df, value_cols=c("pe_max", "v_chosen", "v_entropy"))
#saveRDS(l1_models, "l1_model_cache.rds")
l1_models <- readRDS("../local/l1_model_cache.rds")

subject_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/example_files/mmclock_subject_data.rds") %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder

#test naming collisions
trial_df <- trial_df %>% rename(subid = id) %>% mutate(id=subid)
subject_df <- subject_df %>% rename(subid = id) %>% mutate(id=subid)

gpa <- setup_glm_pipeline(analysis_name="testing", scheduler="slurm", working_directory = tempdir(),
  subject_data=subject_df, run_data=NULL, trial_data=trial_df,
  tr=1.0,
  vm=c(id="subid", run_number="run"),
  fmri_file_regex="nfaswuktm_clock[1-8]_5\\.nii\\.gz",
  fmri_path_regex="clock[0-9]",
  run_number_regex=".*clock([0-9]+)_5.*",
  n_expected_runs=8,
  l1_models=l1_models, bad_ids=c(10637), #test exclusion
  confound_settings=list(
    confound_file="nuisance_regressors.txt",
    confound_columns = c("csf", "dcsf", "wm", "dwm"),
    l1_confound_regressors = c("csf", "dcsf", "wm", "dwm")
  )
)



gpa <- setup_l1_models(gpa)

build_l2_models(gpa$run_data)
