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
source("setup_l2_models.R")
source("fsl_l2_model.R")
source("build_l1_models.R")
source("glm_helper_functions.R")
source("lookup_nifti_inputs.R")
source("get_l1_confounds.R")
source("run_feat_sepjobs.R")
source("cluster_job_submit.R")
# source("build_design_matrix.R")

trial_df <- readRDS("/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/mmy3_trial_df_selective_groupfixed.rds") %>%
  mutate(rt_sec = rt_csv / 1000) %>%
  select(-isi_onset, -iti_onset) %>%
  mutate(session=1) %>%
  dplyr::rename(run_number="run")

# l1_models <- build_l1_models(trial_data=trial_df, value_cols=c("pe_max", "v_chosen", "v_entropy"))
# saveRDS(l1_models, "l1_model_cache.rds")
l1_models <- readRDS("../local/l1_model_cache.rds")
l1_models$models[[1]]$signals <- l1_models$models[[1]]$model_signals
l1_models$models[[1]]$model_signals <- NULL
l1_models$models[[2]]$signals <- l1_models$models[[2]]$model_signals
l1_models$models[[2]]$model_signals <- NULL
l1_models$models[[1]]$regressors <- l1_models$models[[1]]$model_regressors
l1_models$models[[1]]$model_regressors <- NULL
l1_models$models[[2]]$regressors <- l1_models$models[[2]]$model_regressors
l1_models$models[[2]]$model_regressors <- NULL
l1_models$events <- lapply(l1_models$events, function(ee) {
  ee %>% mutate(session=1, run_number=run)
})

l1_models$signals <- lapply(l1_models$signals, function(ss) {
  if (is.data.frame(ss$value)) {
    ss$value <- ss$value %>%
      mutate(session = 1, run_number = run)
  }
  return(ss)
})

#rename from _dt to d_ for derivatives to match bdm
l1_models$models$pe_only$regressors <- c("clock", "feedback", "pe", "d_pe" )
rownames(l1_models$models$pe_only$contrasts) <- colnames(l1_models$models$pe_only$contrasts) <- c("clock", "feedback", "pe", "d_pe")

subject_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/example_files/mmclock_subject_data.rds") %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder

run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/example_files/mmclock_run_data.rds")
#gpa$run_data$..id.. <- NULL
#saveRDS(gpa$run_data, file = "/proj/mnhallqlab/users/michael/fmri.pipeline/example_files/mmclock_run_data.rds")

#test naming collisions
trial_df <- trial_df %>% rename(subid = id) %>% mutate(id=subid)
run_df <- run_df %>% rename(subid = id, run=run_number) %>% mutate(id=subid)
subject_df <- subject_df %>% rename(subid = id) %>% mutate(id=subid)

gpa <- setup_glm_pipeline(analysis_name="testing", scheduler="slurm",
  working_directory = "/proj/mnhallqlab/users/michael/fmri_test",
  subject_data=subject_df, run_data=run_df, trial_data=trial_df,
  tr=1.0,
  vm=c(id="subid", run_number="run"),
  fmri_file_regex="nfaswuktm_clock[1-8]_5\\.nii\\.gz",
  fmri_path_regex="clock[0-9]",
  run_number_regex=".*clock([0-9]+)_5.*",
  n_expected_runs=8,
  l1_models=l1_models, bad_ids=c(10637), #test exclusion
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

#gpa <- build_l1_models(gpa)

# interactive model builder for level 2
gpa <- build_l2_models(gpa)

# interactive model builder for level 3
gpa <- build_l3_models(gpa)


#### setup


# create all FSF files for level one runs
gpa <- setup_l1_models(gpa)

#### execution

jobs <- run_feat_sepjobs(gpa, level=1L)

# todo
# gpa <- verify_lv1_runs(gpa)

# load("test_gpa.RData")

#new nomenclature
# gpa$l1_model_setup$fsl <- gpa$l1_model_setup$fsl %>%
#   dplyr::rename(l1_feat_fsf=feat_file) %>%
#   dplyr::mutate(l1_feat_dir=sub(".fsf", ".feat", l1_feat_fsf, fixed=TRUE))

#save(gpa, file="test_gpa.RData")

# gpa$parallel$fsl$l2_feat_time <- "1:00:00"
# gpa$parallel$fsl$l2_feat_memgb <- "20"
# gpa$parallel$fsl$l3_feat_time <- "1:00:00"
# gpa$parallel$fsl$l3_feat_memgb <- "20"

# setup of l2 models (should follow l1)

gpa <- setup_l2_models(gpa)

jobs <- run_feat_sepjobs(gpa, level=2)

push_pipeline(gpa, l1_model_set = c("pe_only"), l2_model_set = "with_run_number")