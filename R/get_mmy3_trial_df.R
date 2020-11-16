get_mmy3_trial_df <- function(cache_dir="/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline", model="selective", groupfixed=TRUE, allow_cache=TRUE) {
  require(tidyverse)

  expect_file <- file.path(cache_dir, paste0("mmy3_trial_df_", model, "_", ifelse(groupfixed, "groupfixed", "subjspecific"), ".rds"))
  if (allow_cache) {
    if (file.exists(expect_file)) {
      message("Reading cached trial_df from: ", expect_file)
      trial_df <- readRDS(expect_file)
      return(trial_df)
    } else {
      message("Generating trial_df from scratch")
    }
  }
  
  ###
  # SCEPTIC MMClock Y3

  if (model=="fixed") {
    #fixed LR V analysis (not selective maintenance)
    if (groupfixed) {
      #trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_fixed_fixedparams_fmri_ffx_trial_statistics.csv.gz")
      trial_df <- read.csv("/proj/mnhallqlab/users/michael/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_fixed_fixedparams_fmri_ffx_trial_statistics.csv.gz")
    } else {
      trial_df <- read.csv("/proj/mnhallqlab/users/michael/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_fixed_mfx_trial_statistics.csv.gz")
    }
  } else if (model=="fixed_uv") {
    #this is the fixed_uv analysis
    if (groupfixed) {
      trial_df <- read.csv("/proj/mnhallqlab/users/michael/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_fixed_uv_ureset_fixedparams_fmri_ffx_trial_statistics.csv.gz")
    } else {
      trial_df <- read.csv("/proj/mnhallqlab/users/michael/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_fixed_uv_ureset_mfx_trial_statistics.csv.gz")
    }
    
    trial_df <- trial_df %>% 
      mutate( #just for u model
        d_auc=0, #uv reset has no d statistics
        u_chosen_sqrt=sqrt(u_chosen)
      ) %>% group_by(id, run) %>%
      dplyr::mutate(
        u_chosen_quantile = if_else(run_trial==1, NA_real_, u_chosen_quantile),
        u_chosen_z = as.vector(scale(u_chosen)), #z-scored trial
      )
  } else if (model=="selective") {
    #factorized, selective maintenance, equal basis-generalization width
    if (groupfixed) {
      trial_df <- read.csv("/proj/mnhallqlab/users/michael/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz")
    } else {
      trial_df <- read.csv("/proj/mnhallqlab/users/michael/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz") #may have gotten moved

      #old naming scheme
      #trial_df <- read.csv("/proj/mnhallqlab/users/michael/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_mfx_trial_statistics.csv.gz")
      #trial_df <- read.csv("/proj/mnhallqlab/users/michael/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_mfx_trial_statistics.csv.gz")
    }
  }

  trial_df <- trial_df %>%
  mutate(
    run_trial=case_when(
      trial >= 1 & trial <= 50 ~ trial,
      trial >= 51 & trial <= 100 ~ trial - 50L, #dplyr/rlang has gotten awfully picky about data types!!
      trial >= 101 & trial <= 150 ~ trial - 100L,
      trial >= 151 & trial <= 200 ~ trial - 150L,
      trial >= 201 & trial <= 250 ~ trial - 200L,
      trial >= 251 & trial <= 300 ~ trial - 250L,
      trial >= 301 & trial <= 350 ~ trial - 300L,
      trial >= 351 & trial <= 400 ~ trial - 350L,
      TRUE ~ NA_integer_),
    v_entropy_no5=if_else(run_trial <= 5, NA_real_, v_entropy),
    #d_auc_sqrt=if_else(d_auc > 0, NA_real_, sqrt(-1*d_auc)), #only compute the sqrt of d_auc for negative (i.e., reasonable) observations
    v_entropy_sqrt=sqrt(v_entropy),
    rew_om=if_else(score_vba > 0, 1, 0)
  ) %>% #for win/loss maps
  group_by(id, run) %>%
  dplyr::mutate(   #compute rt_swing within run and subject
    rt_vmax_lag = dplyr::lag(rt_vmax, 1, order_by=trial),
    rt_vmax_change = abs(rt_vmax - rt_vmax_lag),
    rt_vmax_change_dir = rt_vmax - rt_vmax_lag,
    v_entropy_lag = dplyr::lag(v_entropy, 1, order_by=trial),
    v_entropy_change = v_entropy - v_entropy_lag, #change in entropy
    v_entropy_change_pos = v_entropy_change*(v_entropy_change > 0),
    v_entropy_change_neg = abs(v_entropy_change*(v_entropy_change < 0)),
    rt_swing = abs( c(NA, diff(rt_csv)))/1000,
    rt_swing_sqrt=sqrt(rt_swing),
    pe_1h=if_else(run_trial <= 25, pe_max, 0), #first-half, second-half split for pe and entropy
    pe_2h=if_else(run_trial > 25, pe_max, 0),
    v_entropy_1h=if_else(run_trial <= 25, v_entropy, 0),
    v_entropy_2h=if_else(run_trial > 25, v_entropy, 0)
  ) %>% ungroup()
  
  #verify the within-run z-scoring
  #library(skimr)
  #trial_df %>% group_by(id, run) %>% select(u_chosen_z) %>% skim()

  lrates <- c(0.05, 0.1, 0.15, 0.2)
  lrates_char <- sub("0.", "p", sprintf("%.2f", lrates), fixed=TRUE)
  
  #add generic action-value under a fixed learning rate (following Badre/Frank)
  trial_df <- trial_df %>% arrange(id, trial) %>% group_by(id) %>% do({
    this_subj <- .
    for (lc in lrates_char) { #initialize zero values for PE and V at each learning rate
      this_subj[[ paste0("v_trial_fixed_", lc) ]] <- 0
      this_subj[[ paste0("pe_trial_fixed_", lc) ]] <- 0
    }

    for (ll in 1:length(lrates)) {
      vcol <- paste0("v_trial_fixed_", lrates_char[ll])
      pecol <- paste0("pe_trial_fixed_", lrates_char[ll])
      
      for (i in 1:nrow(this_subj)) {
        if (i > 1)  { this_subj[ i, vcol ] <- this_subj[i-1, vcol] + lrates[ll]*(this_subj[i-1, "score_csv"] - this_subj[i-1, vcol]) }
        this_subj[i, pecol] <- this_subj[i, "score_csv"] - this_subj[i,vcol]
      }
    }

    this_subj
  }) %>% ungroup()

  saveRDS(trial_df, file=expect_file)
  return(trial_df)
}
