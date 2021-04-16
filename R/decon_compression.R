library(readr)
library(dplyr)
library(tidyr)
library(data.table)

trial_df <- readRDS("/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/mmy3_trial_df_selective_groupfixed.rds") %>%
  mutate(rt_time=clock_onset + rt_csv/1000, #should be pretty redundant with isi_onset, but still
    rt_vmax=rt_vmax/10, #to put into seconds
    rt_vmax_cum=clock_onset + rt_vmax)

#example <- read_csv("deconvolved/sub11350_run1_Schaefer_DorsAttn_2.3mm_deconvolved.csv.gz")
example <- read_csv("/proj/mnhallqlab/users/michael/sceptic_decon/Schaefer_DorsAttn_2.3mm/deconvolved/sub11325_run7_Schaefer_DorsAttn_2.3mm_deconvolved.csv.gz")
example <- example %>% select(vnum, time, decon, atlas_value)
example$decon2 <- example$decon

run_df <- trial_df %>% filter(id==11325 & run==7)

x <- fmri_ts$new(ts_data=example, event_data=run_df, tr=1.0,
  vm=list(value=c("decon", "decon2"), key=c("vnum", "atlas_value")))

rm(example)

#vv <- x$get_ts(orig_names=TRUE)

test <- event_lock_ts(x, event="rt_time", collide_after="rt_time", time_before=0, time_after=12)

int_dt2 <- get_medusa_interpolated_ts(x, event="clock_onset", time_before=-3.0, time_after=3.0,
  collide_before="clock_onset", collide_after=NULL,
  pad_before=-1.5, pad_after=1.5, output_resolution=1.0,
  group_by=c("atlas_value", "trial"))

compress <- get_medusa_compression_score(x, event="clock_onset", time_before=-3.0, time_after=3.0,
  collide_before="clock_onset", collide_after=NULL, group_by=c("atlas_value", "trial"))

compress[, .(tcor = cor(decon_compress_0.9, trial)), by=atlas_value]



test <- lapply(1:length(test), function(vv) { test[[vv]]$trial <- names(test)[vv]; return(test[[vv]]) })
testdf <- bind_rows(test)

rle_test <- rle_dt$new(testdf, keys=c("key1", "key2", "trial", "time"), optimize=FALSE)
rle_df <- rle_test$get()
#  mutate(id=11350, run=1)

ts_data <- list(metadata=list(id=11350, run=1), timedata=example)
run_df <- trial_df %>% filter(id==11350 & run==1)

test <- event_lock_decon(ts_data=example, run_df=run_df, event="clock_onset", collide_before="clock_onset")

test <- example %>% filter(vnum==1)

time_df <- example %>% filter(atlas_value==31) %>% select(-atlas_value)
xtabs(~vnum + time, time_df)

time_df_wide <- time_df %>% pivot_wider(id_cols=c("time"), names_from="vnum", values_from="decon", names_prefix="v") %>%
  select(-time) %>% as.matrix() #don't really need it for now

decomp <- svd(scale(time_df_wide))
cc1 <- cumsum(decomp$d^2/sum(decomp$d^2))
decomp2 <- prcomp(time_df_wide, scale=TRUE)
cc2 <- cumsum(decomp2$sdev^2/sum(decomp2$sdev^2))

#voxels x time vs. time x voxels
decomp <- svd(scale(t(time_df_wide)))

cc1 <- cumsum(decomp$d^2/sum(decomp$d^2))

min(which(cc1 > 0.95))
min(which(cc2 > 0.95))

#for the compressibility -- we're talking about small stretches of decon, like 10 seconds
#so, we want a time x measures/voxels matrix to match general applications
#but, we will end up with many more variables than timepoints
example <- scale(time_df_wide[10:20,])
xx <- eigen(cor(example))



decomp <- svd(example)
cumsum(decomp$d^2/sum(decomp$d^2))


decomp <- svd(t(example))
cumsum(decomp$d^2/sum(decomp$d^2))

cc1 <- cumsum(decomp$d^2/sum(decomp$d^2))

min(which(cc1 > 0.95))

#number of components needed to exceed 0.9
np9 <- approx(y=1:length(decomp$d), x=cc1, xout=0.9)$y
1-(3/11)




#external compression
write.csv(round(example, 3), file="test_decon.csv", row.names=F, col.names=F)


### build glm
trial_df <- readRDS("/proj/mnhallqlab/projects/clock_analysis/fmri/fsl_pipeline/mmy3_trial_df_selective_groupfixed.rds") %>% mutate(rt_sec=rt_csv/1000) %>%
  select(-isi_onset, -iti_onset)
result <- build_l1_model(trial_df, value_cols=c("pe_max", "v_chosen", "v_entropy"))

#test modification
result2 <- build_l1_model(trial_df, l1_model_set=result, value_cols=c("pe_max", "v_chosen", "v_entropy"))
