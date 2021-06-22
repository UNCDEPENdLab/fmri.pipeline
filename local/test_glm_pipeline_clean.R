library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(lgr)
library(foreach)
library(doParallel)
library(dependlab)
library(emmeans)

#library(devtools)
#install_github("UNCDEPENdLab/fmri.pipeline")
library(fmri.pipeline)

trial_df <- readRDS(system.file("example_files/mmy3_trial_df_selective_groupfixed.rds", package="fmri.pipeline")) %>%
  mutate(rt_sec = rt_csv / 1000) %>%
  select(-isi_onset, -iti_onset) %>%
  mutate(session=1) %>%
  dplyr::rename(run_number="run")

# if you want to specify before the gpa object is setup
#l1_models <- build_l1_models(trial_data=trial_df, value_cols=c("pe_max", "v_chosen", "v_entropy"))

subject_df <- readRDS(system.file("example_files/mmclock_subject_data.rds", package="fmri.pipeline")) %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder

run_df <- readRDS(system.file("example_files/mmclock_run_data.rds", package = "fmri.pipeline"))

#test naming collisions
trial_df <- trial_df %>% rename(subid = id) %>% mutate(id=subid)
run_df <- run_df %>% rename(subid = id, run=run_number) %>% mutate(id=subid)
subject_df <- subject_df %>% rename(subid = id) %>% mutate(id=subid)

gpa <- setup_glm_pipeline(analysis_name="testing", scheduler="slurm",
  working_directory = "/proj/mnhallqlab/users/michael/fmri_test",
  subject_data=subject_df, run_data=run_df, trial_data=trial_df,
  tr=1.0, l1_models="prompt", l2_models=NULL, l3_models=NULL,
  vm=c(id="subid", run_number="run"),
  fmri_file_regex="nfaswuktm_clock[1-8]_5\\.nii\\.gz",
  fmri_path_regex="clock[0-9]",
  run_number_regex=".*clock([0-9]+)_5.*",
  n_expected_runs=8,
  bad_ids=c(10637), #test exclusion
  confound_settings=list(
    motion_params_file = "motion.par",
    confound_input_file="nuisance_regressors.txt",
    confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
    l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
    exclude_run = "max(FD) > 4 | sum(FD > .5)/length(FD) > .15", #this must evaluate to a scalar per run
    exclude_subject = "n_good_runs < 4"
  )
)

rm(trial_df)

### model definition

#if you want to build l1 models after model is setup
#gpa <- build_l1_models(gpa)

# interactive model builder for level 2
gpa <- build_l2_models(gpa)

# interactive model builder for level 3
gpa <- build_l3_models(gpa)


### setup of model ingredients: fsf files and timing files

# create all FSF files for level one runs
gpa <- setup_l1_models(gpa)

#### execute all L1 models in parallel by splitting across jobs

jobs <- run_feat_sepjobs(gpa, level=1L)

### setup L2 model ingredients

gpa <- setup_l2_models(gpa)

### execute all L2 models (must follow l1 model execution)

jobs <- run_feat_sepjobs(gpa, level=2)

# L3 is not yet complete

# eventually where this is going
#push_pipeline(gpa, l1_model_set = c("pe_only"), l2_model_set = "with_run_number")