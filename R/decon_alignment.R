library(tidyverse)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)

getMainDir <- function(institution="UNC") {
  if (institution == "UNC") {
    return("/proj/mnhallqlab")
  } else if (institution == "PSU") {
    return("/gpfs/group/mnh5174/default")
  }
}

source("get_mmy3_trial_df.R")
source("event_lock_timeseries.R")

#setwd(file.path(getMainDir(), "projects", "clock_analysis", "fmri", "hippo_voxelwise"))
setwd(file.path(getMainDir(), "users", "michael", "sceptic_decon"))

trial_df <- get_mmy3_trial_df(model="selective", groupfixed=TRUE) %>%
  mutate(rt_time=clock_onset + rt_csv/1000, #should be pretty redundant with isi_onset, but still
    rt_vmax=rt_vmax/10, #to put into seconds
    rt_vmax_cum=clock_onset + rt_vmax)

## trial_df <- read_csv(file.path(getMainDir(), "clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz"))
## trial_df <- trial_df %>%
##   group_by(id, run) %>%  
##   dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), #compute rt_swing within run and subject
##                 rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
##                 rt_lag = lag(rt_csv) ,
##                 omission_lag = lag(score_csv==0),
##                 rt_vmax_lag = lag(rt_vmax),
##                 v_entropy_wi = scale(v_entropy),
##                 run_trial=case_when(
##                   trial >= 1 & trial <= 50 ~ trial,
##                   trial >= 51 & trial <= 100 ~ trial - 50, #dplyr/rlang has gotten awfully picky about data types!!
##                   trial >= 101 & trial <= 150 ~ trial - 100,
##                   trial >= 151 & trial <= 200 ~ trial - 150,
##                   trial >= 201 & trial <= 250 ~ trial - 200,
##                   trial >= 251 & trial <= 300 ~ trial - 250,
##                   trial >= 301 & trial <= 350 ~ trial - 300,
##                   trial >= 351 & trial <= 400 ~ trial - 350,
##                   TRUE ~ NA_real_)) %>% ungroup() %>%
##   dplyr::mutate(rt_csv=rt_csv/1000, rt_vmax=rt_vmax/10) %>% 
##       mutate(rt_vmax_cum=clock_onset + rt_vmax)
    
## base_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise"

## atlas_dirs <- list.dirs(file.path(base_dir, "deconvolved_timeseries","smooth_in_mask_b2015"), recursive=FALSE)
## #atlas_dirs <- c(file.path(base_dir, "deconvolved_timeseries", "long_axis_l_2.3mm"),
## #                file.path(base_dir, "deconvolved_timeseries", "long_axis_r_2.3mm"))

base_dir <- file.path(getMainDir(), "users/michael/sceptic_decon")

#all atlases in use
atlas_dirs <- list.dirs(base_dir, recursive=FALSE)
atlases <- file.path("/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks", paste0(basename(atlas_dirs), ".nii.gz"))

atlas_dirs <- grep("Schaefer", atlas_dirs, value=TRUE)
atlases <- grep("Schaefer", atlases, value=TRUE)

#hard code DAN for a bit
#atlas_dirs <- atlas_dirs[7]
#atlases <- atlases[7]

ncores <- Sys.getenv("ncores")
ncores <- ifelse(ncores == "", 20, as.numeric(ncores))

cl <- makeCluster(ncores)
registerDoParallel(cl)
on.exit(stopCluster(cl))

events <- c("clock_onset", "feedback_onset", "clock_long", "feedback_long", "rt_long", "rt_vmax_cum")
nbins <- 12 #splits along axis

for (a in 1:length(atlas_dirs)) {
  aname <- basename(atlas_dirs[a])

  mask <- oro.nifti::readNIfTI(atlases[a], reorient=FALSE)
  mi <- which(mask > 0, arr.ind = TRUE)
  maskvals <- mask[mi]

  if (isTRUE(checkmate::test_integerish(maskvals))) {
    continuous <- FALSE
    message("Atlas: ", aname, " has integer values. Using unique mask values for aggregation in event-locking.")
    bin_cuts <- NULL
  } else {
    continuous <- TRUE
    message("Atlas: ", aname, "has continuous values. Dividing into ", nbins, "bins")
    bin_cuts <- seq(min(mask[mi]) - 1e-5, max(mask[mi]) + 1e-5, length.out=nbins+1)
  }
  
  afiles <- list.files(file.path(atlas_dirs[a], "deconvolved"), full.names = TRUE)

  for (e in events) {
    if (e == "clock_long") {
      evt_col <- "clock_onset"
      time_before=-5
      time_after=10
    } else if (e == "feedback_long") {
      evt_col <- "feedback_onset"
      time_before=-1
      time_after=10
    } else if (e == "rt_long") {
      evt_col <- "rt_time"
      time_before=-4
      time_after=7
    } else {
      evt_col <- e
      time_before=-3
      time_after=3
    }

    out_name <- paste0(aname, "_", e, "_decon_locked.csv.gz")
    if (file.exists(out_name)) {
      message("Output file already exists: ", out_name)
      next
    }

    elist <- foreach(fname=iter(afiles), .packages = c("dplyr", "readr")) %dopar% {

      #add sub and run for now since I screwed this up in the outputs...
      id <- as.numeric(sub("^.*/sub(\\d+)_.*", "\\1", fname))
      run <- as.numeric(sub("^.*/sub\\d+_run(\\d+).*", "\\1", fname))
      d <- read_csv(fname) %>% select(-atlas_name, -x, -y, -z)
      if (all(is.na(d$decon))) { browser() } #weird problem

      #discretize atlas value into bins (continuous) or unique values (integer mask)
      if (isTRUE(continuous)) {
        d <- d %>% mutate(atlas_value=cut(atlas_value, bin_cuts))
      } else {
        d <- d %>% mutate(atlas_value=factor(atlas_value))
      }
      
      dsplit <- split(d, d$atlas_value)
      rm(d) #garbage cleanup

      subj_df <- trial_df %>% filter(id==!!id & run==!!run)
      
      subj_lock <- tryCatch(event_lock_decon(dsplit, subj_df, event = evt_col, time_before=time_before, time_after=time_after),
        error=function(err) {
          cat("Problems with event locking ", fname, " for event: ", e, "\n  ",
            as.character(err), "\n\n", file="evtlockerrors.txt", append=TRUE); return(NULL)
        }
      )
      return(subj_lock)
    }
    
    all_e <- bind_rows(elist)
    all_e$atlas <- aname
    message("Writing output: ", out_name)
    write_csv(all_e, path=out_name)
  }
}

#Times in the deconvolved files should reflect the +2 seconds for the dropped volumes
#So, this should align with the rt_csv columns appropriately.


