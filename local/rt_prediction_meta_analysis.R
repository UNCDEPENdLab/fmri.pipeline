get_all_coef_df <- function(dir, pattern=".*rt_rs.*RT", emt_number=1, trend_match="trend") {
  coef_files <- list.files(dir, pattern=pattern, full.names = TRUE)
  
  all_coef <- lapply(coef_files, function(x) {
    df <- readRDS(x)
    data <- df$emtrends_list[[paste0("emt_", emt_number)]]
    dupe_names <- duplicated(names(data)) #getting duplicate model_name element
    data <- data[, !dupe_names, with=FALSE] #data.table subset
    data %>% 
      #dplyr::select(-outcome...8, -df, -emt_number, -emt_label, -z.ratio) %>%
      dplyr::select(-df, -emt_number, -emt_label, -z.ratio)
      #dplyr::rename(rew_om=outcome...2)
  })

  all_coef <- all_coef %>%
    rbindlist() %>%
    droplevels()

  all_coef$Pow_fac <- factor(all_coef$Pow)
  all_coef <- all_coef %>% dplyr::rename(rew_om = outcome...2)
  split_coef <- split(all_coef, by = c("Time", "Freq"))
  lens <- sapply(split_coef, length) # there are empty splits due to something about how the original files were written (per AYD)
  split_coef <- split_coef[lens > 0]

  return(split_coef)
}

split_df_rs <- get_all_coef_df("/Users/hallquist/OneDrive/collected_letters/papers/meg/plots/wholebrain/rt/output")
save(split_df_rs, file="split_df_rs.RData")

library(doFuture)

sub_df <- split_df_rs[[1]]

  sub_df <- split_coef[["0.714.f_11.892"]]
  sub_df <- split_coef[["0.714.11.892"]]
  browser()
  m_results <- foreach(sub_df = iter(split_coef), .packages = "brms") %dopar% {
    br_model <- brm(
      rt_csv_sc.trend | se(std.error, sigma = FALSE) ~ 1 + Pow_fac * rew_om + (1 + Pow_fac + rew_om + Pow_fac:rew_om | Sensor),
      sub_df,
      chains = 1, cores = 1, iter = 15000,
      control = list(
        adapt_delta = 0.999,
        max_treedepth = 13
      )
    )

    # does not converge easily (probably because of all of the ~0 ranefs)
    br_model_ffx <- brm(
      rt_csv_sc.trend ~ 1 + Pow_fac * rew_om + (1 + Pow_fac + rew_om + Pow_fac:rew_om | Sensor),
      sub_df,
      chains = 1, cores = 1, iter = 15000,
      control = list(
        adapt_delta = 0.999,
        max_treedepth = 13
      )
    )

    post <- posterior_samples(br_model)

    posterior_means <- emmeans(br_model, ~ Pow_fac | rew_om)

    vv <- gather_emmeans_draws(posterior_means)
    posterior_diff <- pairs(posterior_means) # hi - lo power for reward and omission
    # also need the difference of differences
    return(list(model = br_model, posterior = post, emms = posterior_means, emdiffs = posterior_diff))
  }


test <- get_all_coef_df("/Users/hallquist/OneDrive/collected_letters/papers/meg/plots/wholebrain/rt/output")

brm_loop <- function(coef_list, rhs=NULL, outcome=NULL {
  if (is.null(rhs) { 
    rt_csv_sc.trend ~ 1 + Pow_fac * rew_om + (1 + Pow_fac + rew_om + Pow_fac:rew_om | Sensor),
  }
}
