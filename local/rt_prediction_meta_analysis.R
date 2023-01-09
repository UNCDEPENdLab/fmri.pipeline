# get_all_coef_df <- function(dir, pattern=".*rt_rs.*RT", emt_number=1, trend_match="trend") {
#   coef_files <- list.files(dir, pattern=pattern, full.names = TRUE)
  
#   all_coef <- lapply(coef_files, function(x) {
#     df <- readRDS(x)
#     data <- df$emtrends_list[[paste0("emt_", emt_number)]]
#     dupe_names <- duplicated(names(data)) #getting duplicate model_name element
#     data <- data[, !dupe_names, with=FALSE] #data.table subset
#     data %>% 
#       #dplyr::select(-outcome...8, -df, -emt_number, -emt_label, -z.ratio) %>%
#       dplyr::select(-df, -emt_number, -emt_label, -z.ratio)
#       #dplyr::rename(rew_om=outcome...2)
#   })

#   all_coef <- all_coef %>%
#     rbindlist() %>%
#     droplevels()

#   all_coef$Pow_fac <- factor(all_coef$Pow)
#   all_coef <- all_coef %>% dplyr::rename(rew_om = outcome...2)
#   split_coef <- split(all_coef, by = c("Time", "Freq"))
#   lens <- sapply(split_coef, length) # there are empty splits due to something about how the original files were written (per AYD)
#   split_coef <- split_coef[lens > 0]

#   return(split_coef)
# }

# split_df_rs <- get_all_coef_df("/Users/hallquist/OneDrive/collected_letters/papers/meg/plots/wholebrain/rt/output")
# save(split_df_rs, file="split_df_rs.RData")

# library(doFuture)

# sub_df <- split_df_rs[[1]]

#   sub_df <- split_coef[["0.714.f_11.892"]]
#   sub_df <- split_coef[["0.714.11.892"]]
#   browser()
#   m_results <- foreach(sub_df = iter(split_coef), .packages = "brms") %dopar% {
#     br_model <- brm(
#       rt_csv_sc.trend | se(std.error, sigma = FALSE) ~ 1 + Pow_fac * rew_om + (1 + Pow_fac + rew_om + Pow_fac:rew_om | Sensor),
#       sub_df,
#       chains = 1, cores = 1, iter = 15000,
#       control = list(
#         adapt_delta = 0.999,
#         max_treedepth = 13
#       )
#     )

#     # does not converge easily (probably because of all of the ~0 ranefs)
#     br_model_ffx <- brm(
#       rt_csv_sc.trend ~ 1 + Pow_fac * rew_om + (1 + Pow_fac + rew_om + Pow_fac:rew_om | Sensor),
#       sub_df,
#       chains = 1, cores = 1, iter = 15000,
#       control = list(
#         adapt_delta = 0.999,
#         max_treedepth = 13
#       )
#     )

#     posterior <- posterior_samples(br_model)

#     posterior_means <- emmeans(br_model, ~ Pow_fac | rew_om)

#     vv <- gather_emmeans_draws(posterior_means)
#     posterior_diff <- pairs(posterior_means) # hi - lo power for reward and omission
#     # also need the difference of differences
#     return(list(model = br_model, posterior = posterior, emms = posterior_means, emdiffs = posterior_diff))
#   }


# test <- get_all_coef_df("/Users/hallquist/OneDrive/collected_letters/papers/meg/plots/wholebrain/rt/output")

###

# future.apply::future_lapply(some.values, function(value) {
#   compute_something(value)
# })

# makeClusterFunctionsSlurm(
#   template = "slurm-simple",
#   array.jobs = TRUE,
#   nodename = "localhost",
#   scheduler.latency = 1,
#   fs.latency = 65
# )


load("/proj/mnhallqlab/users/michael/split_df_rs.RData")

# due to weird file naming issues in generation (per Alex), some splits are empty. Drop these up front.
has_zero <- sapply(split_df_rs, nrow) == 0
split_df_rs <- split_df_rs[!has_zero]

library(doFuture)
library(doRNG)
library(foreach)
library(future.batchtools)
library(iterators)
registerDoFuture()
# plan(future.batchtools::batchtools_slurm)

split_df_rs <- lapply(split_df_rs, function(xx) {
  xx$Pow_fac <- factor(xx$Pow_fac, levels = c("-2", "2"), labels = c("Lo", "Hi"))
  xx
})

chunk_size <- 10

#https://tdhock.github.io/blog/2019/future-batchtools/
future::plan(
  future.batchtools::batchtools_slurm,
  template = "slurm-simple",
  resources = list(
    walltime = 40*60*chunk_size, # 40 minutes (specified in seconds), 10 chunks
    memory = 8000, # 8 GB
    ncpus = 4,
    chunks.as.arrayjobs = FALSE
  )
)

#start with simpler problem
#split_df_rs <- split_df_rs[1:12]

res <- foreach(
  df = iter(split_df_rs), .packages=c("brms", "tidybayes", "emmeans", "glue", "dplyr", "ggdist", "data.table"),
  .options.future = list(chunk.size = chunk_size)) %dorng% {

  metadata <- df %>%
    slice(1) %>%
    dplyr::select(.filename, Time, Freq, rhs)
  
  br_model <- brm(
    rt_csv_sc.trend | se(std.error, sigma = FALSE) ~ 1 + Pow_fac * rew_om + (1 + Pow_fac + rew_om + Pow_fac:rew_om | Sensor),
    df,
    chains = 4, cores = 4, iter = 15000,
    control = list(
      adapt_delta = 0.999,
      max_treedepth = 13
    )
  ) # sample_prior=TRUE

  #base_dir <- "/proj/mnhallqlab/projects/Clock_MEG/rt_prediction_brms/model_objs"
  #fname <- file.path(base_dir, make.names(glue("brms_{as.character(df$Freq[1])}_{as.character(df$Time[1])}.RData")))
  #posterior <- posterior_samples(br_model)

  posterior_df <- br_model %>%
    tidy_draws() %>%
    dplyr::select(-ends_with("__"))

  posterior_summaries <- rbindlist(lapply(posterior_df, median_hdci), idcol=TRUE)

  # each object is 1.1 GB, and we have 2851 of them... not worth the disk space right now!
  #save(br_model, posterior_df, file=fname) # detailed model and statistics for later

  # cell means
  condition_means <- emmeans(br_model, ~ Pow_fac * rew_om)

  # don't need these details
  #vv <- gather_emmeans_draws(posterior_means)

  # hi - lo power for reward and omission
  by_reward_diffs <- pairs(condition_means, by="rew_om")

  # difference in rew - omission for high versus low power (2nd order)
  int_contrast <- contrast(condition_means, interaction = "consec")

  # to get posterior predicted values
  # xx <- predict(br_model)

  # to get sensor-specific predictions for all four levels
  # pred_grid <- crossing(
  #   rew_om = unique(df$rew_om),
  #   Pow_fac = unique(df$Pow_fac),
  #   Sensor = unique(df$Sensor),
  #   std.error=0
  # )

  # yy <- posterior_epred(br_model, newdata=pred_grid)

  # yy <- pred_grid %>%
  #   add_epred_draws(br_model, ndraws=1e3)

  # useful for getting linfct
  # dumb <- lm(rt_csv_sc.trend ~ 1 + Pow_fac * rew_om, df)

  # for now, I can't figure out how to combine the group=X approach of brms hypothesis
  # with the ease and convenience of emmeans. The hypothesis function will compute cluster-
  # level estimates of each contrast. For now, let's do that.

  # high power - low power at omission
  hilo_omission_by_sensor <- hypothesis(br_model, "Pow_facHi = 0", group="Sensor", scope="coef")$hypothesis

  # high power - low power at omission
  hilo_reward_by_sensor <- hypothesis(br_model, "Pow_facHi + Pow_facHi:rew_omReward = 0", group="Sensor", scope="coef")$hypothesis

  # high - low, reward - omission
  int_contrast_by_sensor <- hypothesis(br_model, "Pow_facHi:rew_omReward = 0", group = "Sensor", scope = "coef")$hypothesis

  return(list(
    metadata = metadata,
    posterior_summaries = posterior_summaries,
    condition_means = as.data.frame(condition_means),
    by_reward_diffs = as.data.frame(by_reward_diffs),
    int_contrast = as.data.frame(int_contrast),
    hilo_omission_by_sensor = hilo_omission_by_sensor,
    hilo_reward_by_sensor = hilo_reward_by_sensor,
    int_contrast_by_sensor = int_contrast_by_sensor
  ))
}

saveRDS(res, file="mega_brms_results.rds")

# brm_loop <- function(coef_list, rhs=NULL, outcome=NULL {
#   if (is.null(rhs) { 
#     rt_csv_sc.trend ~ 1 + Pow_fac * rew_om + (1 + Pow_fac + rew_om + Pow_fac:rew_om | Sensor),
#   }
# }

