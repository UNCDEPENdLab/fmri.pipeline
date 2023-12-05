library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(lgr)
library(foreach)
library(doParallel)
library(emmeans)
library(glue)
library(oro.nifti)
library(RNifti)
# library(fmri.pipeline)

# setwd("/proj/mnhallqlab/users/michael/fmri.pipeline/R")
# setwd("~/Data_Analysis/r_packages/fmri.pipeline/R")
file.sources <- list.files("~/Data_Analysis/r_packages/fmri.pipeline/R",
  pattern = "*.R$", full.names = TRUE,
  ignore.case = TRUE
)
sapply(file.sources, source, .GlobalEnv)

# source("setup_glm_pipeline.R")
# source("finalize_pipeline_configuration.R")
# source("glm_helper_functions.R")
# source("fsl_helper_functions.R")
# source("fsl_l1_model.R")
# source("fsl_l2_model.R")
# source("fsl_l3_model.R")
# source("setup_l1_models.R")
# source("setup_l2_models.R")
# source("setup_l3_models.R")
# source("specify_contrasts.R")
# source("build_l2_models.R")
# source("build_l1_models.R")
# source("glm_helper_functions.R")
# source("lookup_nifti_inputs.R")
# source("get_l1_confounds.R")
# source("truncate_runs.R")
source("run_feat_sepjobs.R")
# source("cluster_job_submit.R")
# source("insert_df_sqlite.R")
# source("read_df_sqlite.R")
source("get_output_directory.R")
source("populate_l1_from_yaml.R")
source("l1_helper_functions.R")
 source("build_design_matrix.R")
source("fmri_utility_fx.R")
source("build_fwe_correction.R")
source("ptfce_spec.R")
source("R_batch_job.R")
source("voxelwise_deconvolution.R")
source("event_lock_timeseries.R")
source("fmri_ts.R")

d <- data.table::fread("/proj/mnhallqlab/projects/medusa_simulation/simulated_data/batch026/rw_ntrial30_samp5_smid5_tr1_snr1_nC5_nT765_cnv0/rep01/medusa_preanalyzed/outputFile_basicTrial_samp100_tr0dot5_snr1_nC5_nT150_rep1_subject_001_COMP_MASK/deconvolved/sub002_deconvolved.csv.gz")

decon_dat <- d %>%
  select(vnum, time, decon, atlas_value) %>%
  mutate(atlas_value = round(atlas_value) - 1) %>%
  data.table() # remove weird rounding error issues from matlab

event_data_tot <- get(load("/proj/mnhallqlab/projects/medusa_simulation/simulated_data/batch026/rw_ntrial30_samp5_smid5_tr1_snr1_nC5_nT765_cnv0/task_designs/rep01/td_rw_ntrial30_samp5_smid5_tr1_snr1_nC5_nT765_cnv0_subject_002.RData"))
event_data <- event_data_tot$ev_times %>%
  mutate(id = "002") %>%
  select(id, event, run, trial, start_vol) %>%
  spread(event, start_vol) %>%
  tibble()
ev_vol <- event_data_tot$events$feedback$act_func$time
time_be <- ev_vol[1]
time_af <- ev_vol[length(ev_vol)]

fmri_event_data <- fmri_ts$new(
  ts_data = decon_dat, event_data = event_data, tr = 1.0,
  vm = list(value = c("decon"), key = c("vnum", "atlas_value"))
)

interp_dt <- get_medusa_interpolated_ts(fmri_event_data,
  event = "feedback", time_before = time_be, time_after = time_af,
  collide_before = "feedback", collide_after = NULL,
  pad_before = -1.5, pad_after = 1.5, output_resolution = metadata$TR,
  group_by = c("atlas_value", "trial")
)


ff <- system.file("example_files/mmy3_trial_df_selective_groupfixed.rds", package = "fmri.pipeline")
#trial_df <- readRDS("/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/mmy3_trial_df_selective_groupfixed.rds") %>%
trial_df <- readRDS(ff) %>%
  mutate(rt_sec = rt_csv / 1000) %>%
  select(-isi_onset, -iti_onset)
  #mutate(session=1) %>%
  #dplyr::rename(run_number="run")

# l1_models <- build_l1_models(trial_data=trial_df, value_cols=c("pe_max", "v_chosen", "v_entropy"))
# saveRDS(l1_models, "l1_model_cache.rds")
# l1_models <- readRDS("../local/l1_model_cache.rds")
# l1_models$models[[1]]$signals <- l1_models$models[[1]]$model_signals
# l1_models$models[[1]]$model_signals <- NULL
# l1_models$models[[2]]$signals <- l1_models$models[[2]]$model_signals
# l1_models$models[[2]]$model_signals <- NULL
# l1_models$models[[1]]$regressors <- l1_models$models[[1]]$model_regressors
# l1_models$models[[1]]$model_regressors <- NULL
# l1_models$models[[2]]$regressors <- l1_models$models[[2]]$model_regressors
# l1_models$models[[2]]$model_regressors <- NULL
# l1_models$events <- lapply(l1_models$events, function(ee) {
#   ee %>% mutate(session=1, run_number=run)
# })

# l1_models$signals <- lapply(l1_models$signals, function(ss) {
#   if (is.data.frame(ss$value)) {
#     ss$value <- ss$value %>%
#       mutate(session = 1, run_number = run)
#   }
#   return(ss)
# })

#rename from _dt to d_ for derivatives to match bdm
# l1_models$models$pe_only$regressors <- c("clock", "feedback", "pe", "d_pe" )
# rownames(l1_models$models$pe_only$contrasts) <- colnames(l1_models$models$pe_only$contrasts) <- c("clock", "feedback", "pe", "d_pe")

ff <- system.file("example_files/mmclock_subject_data.rds", package = "fmri.pipeline")
subject_df <- readRDS(ff) %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder

ff <- system.file("example_files/mmclock_run_data.rds", package = "fmri.pipeline")
run_df <- readRDS(ff)
#gpa$run_data$..id.. <- NULL
#saveRDS(gpa$run_data, file = "/proj/mnhallqlab/users/michael/fmri.pipeline/example_files/mmclock_run_data.rds")

#test naming collisions
# trial_df <- trial_df %>% rename(subid = id) %>% mutate(id=subid)
# run_df <- run_df %>% rename(subid = id) %>% mutate(id=subid)
# subject_df <- subject_df %>% rename(subid = id) %>% mutate(id=subid)

run_df$tr <- sample(c(1, 2), nrow(run_df), replace=TRUE)

subject_df <- subject_df %>% dplyr::slice(1:8)
trial_df <- trial_df %>% filter(id %in% subject_df$id)

gpa <- setup_glm_pipeline(analysis_name="testing", scheduler="slurm",
  output_directory = "/proj/mnhallqlab/users/michael/fmri_test",
  subject_data=subject_df, run_data=run_df, trial_data=trial_df,
  tr=NULL,
  vm=c(id="subid"),
  fmri_file_regex="nfaswuktm_clock[1-8]_5\\.nii\\.gz",
  fmri_path_regex="clock[0-9]",
  run_number_regex=".*clock([0-9]+)_5.*",
  n_expected_runs=8,
  l1_models="prompt", #bad_ids=c(10637), #test exclusion
  drop_volumes=2,
  confound_settings=list(
    motion_params_file = "motion.par",
    confound_input_file="nuisance_regressors.txt",
    confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
    l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
    exclude_run = "max(FD) > 5 | sum(FD > .9)/length(FD) > .10", #this must evaluate to a scalar per run
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)", # time, volume, last_onset, last_offset
    exclude_subject = "n_good_runs < 4",
    spike_volumes = "FD > 0.9"
  ),
  # output_locations=list(
  #   feat_l1_directory="{gpa$output_directory}/id-{id}"
  # )
)

#gpa$confound_settings$truncate_run <- "FD > 0.9 & time > last_offset"
#gpa$confound_settings$run_truncation_columns <- c("FD", "time", "last_offset")
#rm(trial_df)

# cached above
gpa <- build_l1_models(gpa)


gpa2 <- finalize_pipeline_configuration(gpa)


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

# gpa$parallel$fsl$l2_feat_time <- "1:00:00"
# gpa$parallel$fsl$l2_feat_memgb <- "20"
# gpa$parallel$fsl$l3_feat_time <- "1:00:00"
# gpa$parallel$fsl$l3_feat_memgb <- "20"

# setup of l2 models (should follow l1)

gpa <- setup_l2_models(gpa)

#save(gpa, file="gpa_tmp_9Jul2021.RData")
load(file="gpa_tmp_9Jul2021.RData")
gpa$l2_models$models$l2$by_subject$cope_list <- lapply(
  gpa$l2_models$models$l2$by_subject$cope_list, function(xx) {
    xx %>% dplyr::rename(l2_cope_number = l2_cope, l2_cope_name = l2_cope_names)
  }
)


jobs <- run_feat_sepjobs(gpa, level=2)

gpa$parallel$fsl$l3_feat_cpusperjob <- 16

gpa <- setup_l3_models(gpa)

jobs <- run_feat_sepjobs(gpa, level = 3)

push_pipeline(gpa, l1_model_set = c("pe_only"), l2_model_set = "with_run_number")


####

#afni clustsim tests

setwd("/proj/mnhallqlab/studies/MMClock/MR_Proc/10637_20140304/mni_5mm_aroma/sceptic_vchosen_ventropy_dauc_pemax_vtime_preconvolve")
res4d_files <- list.files(pattern = "res4d.nii.gz", getwd(), full.names = T, recursive = T)
fwhmx_mask_files <- list.files(pattern = "mask.nii.gz", getwd(), full.names = T, recursive = T)

mytest <- clustsim_spec$new(
  fwhmx_input_files = res4d_files, fwhmx_mask_files = fwhmx_mask_files, scheduler = "slurm", prefix = "test", out_dir = "/proj/mnhallqlab/users/michael/fmri.pipeline/local",
  clustsim_mask = "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii", ncpus=8
)
mytest$run()




# look at old/new problems with regressor scaling, outlier detection in FLAME
old <- "/proj/mnhallqlab/studies/MMClock/MR_Proc/11338_20141213/mni_5mm_aroma/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/designmatrix.RData"
new <- "/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/feat_l1/sub-11338/entropy/entropy_bdm_setup.RData"
load(old)
load(new)

# clock regressor
old_clock <- do.call(rbind, d$design[, "clock"])
new_clock <- do.call(rbind, d_obj$design[, "clock"])
summary(old_clock[, "onset"] - new_clock[, "onset"])
summary(old_clock[, "duration"] - new_clock[, "duration"])

# feedback regressor
old_feedback <- do.call(rbind, d$design[, "feedback"])
new_feedback <- do.call(rbind, d_obj$design[, "feedback"])
summary(old_feedback[, "onset"] - new_feedback[, "onset"])
summary(old_feedback[, "duration"] - new_feedback[, "duration"])

# entropy regressor
old_v_entropy <- data.frame(do.call(rbind, d$design[, "v_entropy"]))
new_v_entropy <- data.frame(do.call(rbind, d_obj$design[, "entropy_clock"]))

# small trial gap at trial 123? Oh, this is a missed trial (rt > 4!) -- so, new is good on this
which(diff(new_v_entropy[, "trial"]) == 2)
old_v_entropy[120:125, ]
new_v_entropy[120:125, ]

old_v_entropy <- old_v_entropy[-123, ]

summary(old_v_entropy[, "onset"] - new_v_entropy[, "onset"])
summary(old_v_entropy[, "duration"] - new_v_entropy[, "duration"])
summary(old_v_entropy[, "value"] - new_v_entropy[, "value"]) # good

# convolved regressors
str(d_obj$design_convolved)
str(d$design_convolved)

library(dplyr)
# old had the excess shortening by one TR (bad truncation calculation) -- need to drop 1 more volume in new to match
new_convolved <- do.call(bind_rows, lapply(d_obj$design_convolved, function(x) x[-nrow(x), ])) # drop last row
old_convolved <- do.call(bind_rows, d$design_convolved)

sapply(1:ncol(new_convolved), function(x) { cor(new_convolved[,x], old_convolved[,x])})

# looks 100% fine
sapply(1:ncol(new_convolved), function(x) {
  summary(new_convolved[, x] - old_convolved[, x])
})

# basically no change.
summary(new_convolved[, "entropy_clock"] - old_convolved[, "v_entropy"])

# convolved inputs to FEAT
new_evs <- lapply(
  list.files(path = "/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/feat_l1/sub-11338/entropy/FEAT_LVL1_run1.feat/custom_timing_files", full.names = T), 
  function(x) {
    as.numeric(readLines(x))
  }
)

new_evs <- lapply(new_evs, function(x) {
  x[-length(x)]
}) # drop last element for alignment

old_evs <- lapply(
  list.files(path = "/proj/mnhallqlab/studies/MMClock/MR_Proc/11338_20141213/mni_5mm_aroma/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/FEAT_LVL1_run1.feat/custom_timing_files", full.names = T),
  function(x) {
    as.numeric(readLines(x))
  }
)

# no meaningful differences in the inputs to FEAT LVL1 run1
cor(old_evs[[1]], new_evs[[1]])
summary(old_evs[[1]] - new_evs[[1]])
summary(old_evs[[1]])

cor(old_evs[[2]], new_evs[[2]])
summary(old_evs[[2]] - new_evs[[2]])
summary(old_evs[[2]])

cor(old_evs[[3]], new_evs[[3]])
summary(old_evs[[3]] - new_evs[[3]])
summary(old_evs[[3]])

# check confounds -- 100% match

old_confounds <- read.table("/proj/mnhallqlab/studies/MMClock/MR_Proc/11338_20141213/mni_5mm_aroma/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/FEAT_LVL1_run1.feat/confoundevs.txt")
new_confounds <- read.table("/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/feat_l1/sub-11338/entropy/FEAT_LVL1_run1.feat/confoundevs.txt")[-268,]

lapply(names(old_confounds), function(x) cor(old_confounds[[x]], new_confounds[[x]]))
