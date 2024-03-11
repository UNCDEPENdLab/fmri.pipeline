library(dplyr)
library(fmri.pipeline)
# setwd(file.path(getMainDir(), "projects", "clock_analysis", "fmri", "hippo_voxelwise"))
setwd(file.path(getMainDir(), "users", "michael", "sceptic_decon"))

trial_df <- get_mmy3_trial_df(model="selective", groupfixed=TRUE) %>%
  mutate(rt_time=clock_onset + rt_csv/1000, #should be pretty redundant with isi_onset, but still
    rt_vmax=rt_vmax/10, #to put into seconds
    rt_vmax_cum=clock_onset + rt_vmax)

## trial_df <- read_csv(file.path(getMainDir(), "clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz"))
    
## base_dir <- "/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise"

## atlas_dirs <- list.dirs(file.path(base_dir, "deconvolved_timeseries","smooth_in_mask_b2015"), recursive=FALSE)
## #atlas_dirs <- c(file.path(base_dir, "deconvolved_timeseries", "long_axis_l_2.3mm"),
## #                file.path(base_dir, "deconvolved_timeseries", "long_axis_r_2.3mm"))

base_dir <- file.path(getMainDir(), "users/michael/sceptic_decon")

#all atlases in use
#atlas_dirs <- list.dirs(base_dir, recursive=FALSE)

# schaefer 200 -> 400 MEDuSA remapping process
atlas_dirs <- "/proj/mnhallqlab/users/michael/sceptic_decon/Schaefer_400_remap"
atlases <- file.path("/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks", paste0(basename(atlas_dirs), ".nii.gz"))

atlas_dirs <- grep("Schaefer", atlas_dirs, value=TRUE)
atlases <- grep("Schaefer", atlases, value=TRUE)

#nbins <- 12 # splits along axis

# atlas_dirs <- list.dirs(base_dir, recursive = FALSE)[c(1,3,4)]
# atlases <- file.path("/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks", paste0(basename(atlas_dirs), ".nii.gz"))

# these are all ROIs in the Schaefer 400 that are in the DorsAttn (7) or DorsAttnA + DorsAttnB (17)
# node numbers are based on the 7-network order, consistent with the 200-400 remapping in concatenate_decon_remap.R
# NB. Jul 2022: This is missing 3 nodes that were retained in the N=47
# dan_400_7net_rois <- c(
#   1L, 8L, 28L, 31L, 46L, 58L, 69L, 70L, 71L, 72L, 73L, 74L, 75L,
#   76L, 77L, 78L, 79L, 80L, 81L, 82L, 83L, 84L, 85L, 86L, 87L, 88L,
#   89L, 90L, 91L, 94L, 201L, 209L, 216L, 227L, 230L, 247L, 267L,
#   271L, 272L, 273L, 274L, 275L, 276L, 277L, 278L, 279L, 280L, 281L,
#   282L, 283L, 284L, 285L, 286L, 287L, 288L, 289L, 290L, 291L, 292L,
#   293L
# )

dan_labels <- readxl::read_excel("/proj/mnhallqlab/projects/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH DAN Labels 400 Good Only 47 parcels.xlsx")

dan_400_7net_rois <- dan_labels$roi7_400

alignments <- list(
  clock_long = list(
    evt_col = "clock_onset",
    time_before = -5,
    time_after = 10,
    collide_before = "iti_onset", # censor data if we bump into the end of the prior trial
    collide_after = "clock_onset" # censor data it we hit the next trial
  ),
  rt_long = list(
    evt_col = "rt_time",
    time_before = -4,
    time_after=7,
    collide_before = "iti_onset",
    collide_after = "clock_onset"
  )
)

run_decon_alignment(atlases, base_dir, trial_df,
  alignments = alignments, aggregate_by = "roi_400", atlas_subset = dan_400_7net_rois,
  overwrite = TRUE, tr = 1.0, ncpus = 16, mem_per_cpu="12g", walltime = "20:00:00", scheduler = "slurm"
)

# long axis hippo

atlas_dirs <- c(
  "/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/deconvolved_timeseries/smooth_in_mask/long_axis_l_2.3mm",
  "/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/deconvolved_timeseries/smooth_in_mask/long_axis_r_2.3mm"
)

atlases <- c(
  "/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/masks/long_axis_l_2.3mm.nii.gz",
  "/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/masks/long_axis_r_2.3mm.nii.gz"
)

base_dir <- "/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/deconvolved_timeseries/smooth_in_mask"

run_decon_alignment(atlases, base_dir, trial_df,
  alignments = alignments, nbins=12,
  overwrite = TRUE, tr = 1.0, ncpus = 16, mem_per_cpu = "12g", walltime = "20:00:00", scheduler = "slurm"
)


##### OLD

library(tidyverse)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(data.table)
library(checkmate)

getMainDir <- function(institution="UNC") {
  if (institution == "UNC") {
    return("/proj/mnhallqlab")
  } else if (institution == "PSU") {
    return("/gpfs/group/mnh5174/default")
  }
}

setwd("/proj/mnhallqlab/users/michael/fmri.pipeline/R")
source("get_mmy3_trial_df.R")
source("event_lock_timeseries.R")
source("fmri_ts.R")


# hard code DAN for a bit
atlas_dirs <- atlas_dirs[3]
atlases <- atlases[3]

#events <- c("clock_onset", "feedback_onset", "clock_long", "feedback_long", "rt_long", "rt_vmax_cum")
events <- c("whole_trial", "rt_to_rt", "rt8")


ncores <- Sys.getenv("ncores")
ncores <- ifelse(ncores == "", 20, as.numeric(ncores))

cl <- makeCluster(ncores)
registerDoParallel(cl)
on.exit(stopCluster(cl))

# registerDoSEQ()


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
  
  d_files <- list.files(file.path(atlas_dirs[a], "deconvolved"), full.names = TRUE)

  pad_before <- -1.5 #seconds before earliest event, to aid in interpolation
  pad_after <- 1.5
  collide_before <- collide_after <- NULL #default to no censoring
  
  for (e in events) {
    if (e == "clock_long") {
      evt_col <- "clock_onset"
      time_before <- -5
      time_after <- 10
      collide_before <- "clock_onset" #censor data at trial boundaries (prior clock onset, next clock onset)
      collide_after <- "clock_onset"
    } else if (e == "feedback_long") {
      evt_col <- "feedback_onset"
      time_before=-1
      time_after=10
    } else if (e == "whole_trial") { #from clock onset through ITI
      evt_col <- "clock_onset"
      time_before=0
      time_after=12
      collide_before <- NULL #irrelevant
      collide_after <- "clock_onset"
    } else if (e == "rt_to_rt") { #from current response to next one
      evt_col <- "rt_time"
      time_before=0
      time_after=12
      collide_before <- NULL #irrelevant
      collide_after <- "rt_time"
    } else if (e == "rt_long") {
      evt_col <- "rt_time"
      time_before=-4
      time_after=7
    } else if (e == "rt8") { #symmetric 16-second windows around rt, no collision
      evt_col <- "rt_time"
      time_before=-8
      time_after=8
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
    
    elist <- foreach(fname=iter(d_files), .packages = c("dplyr", "readr", "data.table")) %dopar% {

      #add sub and run for now since I screwed this up in the outputs...
      id <- as.numeric(sub("^.*/sub(\\d+)_.*", "\\1", fname))
      run <- as.numeric(sub("^.*/sub\\d+_run(\\d+).*", "\\1", fname))
      d <- read_csv(fname) %>% select(-atlas_name, -x, -y, -z)
      if (all(is.na(d$decon))) { browser() } #weird problem

      #discretize atlas value into bins (continuous) or unique values (integer mask)
      if (isTRUE(continuous)) {
        d <- d %>% mutate(atlas_value=cut(atlas_value, bin_cuts))
      } else {
        #d <- d %>% mutate(atlas_value=factor(atlas_value))
        d <- d %>% mutate(atlas_value=as.numeric(atlas_value))
      }
      
      #run data
      subj_df <- trial_df %>% filter(id==!!id & run==!!run)

      tsobj <- fmri_ts$new(ts_data=d, event_data=subj_df, tr=1.0,
        vm=list(value=c("decon"), key=c("vnum", "atlas_value")))

      ## dsplit <- split(d, d$atlas_value)
      ## rm(d) #garbage cleanup

      ## subj_lock <- tryCatch(event_lock_decon(dsplit, subj_df, event = evt_col, time_before=time_before, time_after=time_after),
      ##   error=function(err) {
      ##     cat("Problems with event locking ", fname, " for event: ", e, "\n  ",
      ##       as.character(err), "\n\n", file="evtlockerrors.txt", append=TRUE); return(NULL)
      ##   }
      ## )

      subj_lock <- tryCatch(get_medusa_interpolated_ts(tsobj, event=evt_col,
        time_before=time_before, time_after=time_after,
        collide_before=collide_before, collide_after=collide_after,
        pad_before=pad_before, pad_after=pad_after, output_resolution=1.0,
        group_by=c("atlas_value", "trial")), #one time series per region and trial
        
        error=function(err) {
          cat("Problems with event locking ", fname, " for event: ", e, "\n  ",
            as.character(err), "\n\n", file="evtlockerrors.txt", append=TRUE); return(NULL)
        }
      )
      
      if (!is.null(subj_lock)) { 
        #tack on run and id
        subj_lock[,id := id]
        subj_lock[,run := run]
      }

      compress <- tryCatch(get_medusa_compression_score(tsobj, event=evt_col,
        time_before=time_before, time_after=time_after,
        collide_before=collide_before, collide_after=collide_after,
        group_by=c("atlas_value", "trial")),
        
        error=function(err) {
          cat("Problems with compression of ", fname, " for event: ", e, "\n  ",
            as.character(err), "\n\n", file="evtlockerrors.txt", append=TRUE); return(NULL)
        }
      )

      if (!is.null(compress)) { 
        #tack on run and id
        compress[, id := id]
        compress[, run := run]
      }
      
      return(list(ts=subj_lock, compress=compress))
    }
    
    all_e <- bind_rows(lapply(elist, "[[", "ts"))
    all_e$atlas <- aname
    message("Writing output: ", out_name)
    write_csv(all_e, path=out_name)

    all_c <- bind_rows(lapply(elist, "[[", "compress"))
    all_c$atlas <- aname
    com_out <- sub("_decon_locked.csv.gz", "_compression.csv.gz", out_name, fixed=TRUE)
    message("Writing output: ", com_out)
    write_csv(all_c, path=com_out)
  }
}

#Times in the deconvolved files should reflect the +2 seconds for the dropped volumes
#So, this should align with the rt_csv columns appropriately.
