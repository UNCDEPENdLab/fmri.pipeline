library(dplyr)
library(emmeans)
library(fmri.pipeline)
library(readr)

source("/proj/mnhallqlab/projects/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
trial_df <- get_trial_data(repo_directory = "/proj/mnhallqlab/projects/clock_analysis", dataset = "mmclock_fmri", groupfixed = TRUE) %>%
  dplyr::rename(run_number = "run")

subject_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_subject_data.rds") %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder

# run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_run_data.rds")

gpa <- setup_glm_pipeline(analysis_name="mmclock_nov2021", scheduler="slurm",
  output_directory = "/proj/mnhallqlab/users/michael/mmclock_pe",
  subject_data=subject_df, trial_data=trial_df, # run_data=run_df,
  tr=1.0, drop_volumes = 2,
  n_expected_runs=8,
  l1_models=NULL, l2_models=NULL, l3_models=NULL,
  fmri_file_regex="nfaswuktm_clock[1-8]_5\\.nii\\.gz",
  fmri_path_regex="clock[0-9]",
  run_number_regex=".*clock([0-9]+)_5.*",
  confound_settings=list(
    motion_params_file = "motion.par",
    confound_input_file="nuisance_regressors.txt",
    confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
    l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
    exclude_run = "max(FD) > 5 | sum(FD > .9)/length(FD) > .10", #this must evaluate to a scalar per run
    exclude_subject = "n_good_runs < 4",
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL
  ),
  parallel=list(
    fsl=list(
      l1_feat_alljobs_time="144:00:00"

    ),
    

  )
)

gpa <- build_l1_models(gpa, from_spec_file = "pe_all.yaml")

gpa$parallel$fsl$l1_feat_alljobs_time <- "144:00:00"




gpa <- build_l2_models(gpa)
gpa <- build_l3_models(gpa)

run_glm_pipeline(gpa)




# trial_df$dummy <- rnorm(nrow(trial_df))
# stats_df <- trial_df %>%
#   group_by(id, run_number) %>%
#   do({
#     df <- .
#     res <- lm(dummy ~ rew_om_c + abspe, df)
#     vv <- tryCatch(car::vif(res), error = function(e) {
#       return(c(rew_om_c=NA, abspe_c=NA, `rew_om_c:abspec_c`=NA))
#     })
#     data.frame(run_number=df$run_number[1], rewFunc=df$rewFunc[1], vif=vv, eff=names(vv))
#   }) %>%
#   ungroup()

# library(afex)
# aov_ez(id = "id", data = stats_df, dv = "vif", within = "rewFunc")

# mm <- lmer(vif ~ rewFunc*eff + (1|id), stats_df)
# emmeans(mm, ~eff | rewFunc)

# vif check
# trial_df %>% dplyr::filter(id == 10638) -> test
# test$dummy <- rnorm(400)
# car::vif(lm(dummy ~ abspe * rew_om, test)) #BAD
# car::vif(lm(dummy ~ abspe_c * rew_om_c, test)) #GOOD