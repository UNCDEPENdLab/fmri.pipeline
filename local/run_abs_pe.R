library(dplyr)
library(fmri.pipeline)

trial_df <- readRDS("/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/mmy3_trial_df_selective_groupfixed.rds") %>%
  mutate(rt_sec = rt_csv / 1000) %>%
  select(-isi_onset, -iti_onset) %>%
  mutate(session = 1, abspe=abs(pe_max), absxrew=abspe*rew_om) %>%
  dplyr::rename(run_number = "run") %>%
  dplyr::select(
    -starts_with("pe_trial_fixed"), -starts_with("v_trial_fixed"),
    -pe_1h, -pe_2h, -v_entropy_1h, -v_entropy_2h,
    -kld_newlearn, -kld_forget, -intrinsic_discrepancy
  )

subject_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_subject_data.rds") %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder

run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_run_data.rds")

# kld_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/local/kld_df_aug2021.rds") %>%
#   dplyr::select(id, run, trial, starts_with("kld")) %>%
#   dplyr::rename(run_number = run) %>%
#   mutate(
#     log_kld3_cum2 = log(kld3_cum2 + .00001),
#     log_kld3 = log(kld3 + .00001)
#   )

# trial_df <- trial_df %>% left_join(kld_df, by = c("id", "run_number", "trial"))

gpa <- setup_glm_pipeline(analysis_name="mmclock_aug2021", scheduler="slurm",
  output_directory = "/proj/mnhallqlab/users/michael/abspe",
  subject_data=subject_df, run_data=run_df, trial_data=trial_df,
  tr=1.0,
  n_expected_runs=8,
  l1_models=NULL, l2_models=NULL, l3_models=NULL,
  confound_settings=list(
    motion_params_file = "motion.par",
    confound_input_file="nuisance_regressors.txt",
    confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
    l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
    exclude_run = "max(FD) > 5 | sum(FD > .5)/length(FD) > .15", #this must evaluate to a scalar per run
    exclude_subject = "n_good_runs < 4"
  )
)


gpa <- build_l1_models(gpa)
gpa <- build_l2_models(gpa)
gpa <- build_l3_models(gpa)

run_glm_pipeline(gpa)
