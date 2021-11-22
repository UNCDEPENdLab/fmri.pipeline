library(dplyr)
library(emmeans)
library(fmri.pipeline)
library(readr)

source("/proj/mnhallqlab/projects/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
trial_df <- get_trial_data(repo_directory = "/proj/mnhallqlab/projects/clock_analysis", dataset = "mmclock_fmri", groupfixed = TRUE) %>%
  arrange(id, run, trial) %>%
  select(-isi_onset, -iti_onset) %>%
  mutate(
    session = 1, rew_om = score_csv > 0, rew_om_c = rew_om - 0.5,
  log_kld3_cum2 = log(kld3_cum2 + .00001), log_kld3 = log(kld3 + .00001)) %>%
  dplyr::rename(run_number = "run") %>%
  group_by(id, run_number) %>%
  mutate(abs_pe_c = abs_pe - mean(abs_pe, na.rm=TRUE), abspexrew = abs_pe_c*rew_om_c) %>%
  ungroup()


# trial_df_new %>%
#   select(trial, v_entropy, v_entropy_wi, v_entropy_wi_change) %>%
#   head()

# trial_df %>%
#   select(trial, v_entropy, v_entropy_wi, v_entropy_change) %>%
#   head()


# N.B. Looks like there is a discrepancy between old and new trial_df

# old have entropy_change as entropy(t) - entropy(t-1) -- so, change on this trial compared to prior
# trial_df <- readRDS("/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/mmy3_trial_df_selective_groupfixed.rds") %>%
#   group_by(id, run) %>%
#   mutate(v_entropy_wi = as.vector(scale(v_entropy))) %>%
#   ungroup() %>%
#   arrange(id, run, trial) %>%
#   mutate(rt_sec = rt_csv / 1000) %>%
#   select(-isi_onset, -iti_onset) %>%
#   mutate(session = 1, abspe = abs(pe_max), rew_om_c = rew_om - 0.5) %>%
#   dplyr::rename(run_number = "run") %>%
#   group_by(id, run_number) %>%
#   mutate(abspe_c = abspe - mean(abspe, na.rm=TRUE), abspexrew = abspe_c*rew_om_c) %>%
#   ungroup() %>%
#   dplyr::select(
#     -starts_with("pe_trial_fixed"), -starts_with("v_trial_fixed"),
#     -pe_1h, -pe_2h, -v_entropy_1h, -v_entropy_2h,
#     -kld_newlearn, -kld_forget, -intrinsic_discrepancy
#   )

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

subject_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_subject_data.rds") %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder

# run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_run_data.rds")

# kld_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/local/kld_df_aug2021.rds") %>%
#   dplyr::select(id, run, trial, starts_with("kld")) %>%
#   dplyr::rename(run_number = run) %>%
#   mutate(
#     log_kld3_cum2 = log(kld3_cum2 + .00001),
#     log_kld3 = log(kld3 + .00001)
#   )
#
#trial_df <- trial_df %>% left_join(kld_df, by = c("id", "run_number", "trial"))

gpa <- setup_glm_pipeline(analysis_name="mmclock_nov2021_fixed", scheduler="slurm",
  output_directory = "/proj/mnhallqlab/users/michael/abspe",
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
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)", # large FD after last offset, or time after the task ended
    spike_volumes = NULL
  )
)


gpa <- build_l1_models(gpa)

#gpa$parallel

#save(gpa, file="gpa_alll1_10Nov2021.RData")
#load(file="abspe_gpa_13Oct2021.RData")
gpa <- build_l2_models(gpa)
gpa <- build_l3_models(gpa)
save(gpa, file="all_sceptic_17Nov2021.RData")
#gpa <- finalize_pipeline_configuration(gpa)
#gpa$drop_volumes <- 2
#gpa$parallel$finalize_time <- "5:00:00"
run_glm_pipeline(gpa)
#gpa <- setup_l1_models(gpa)
