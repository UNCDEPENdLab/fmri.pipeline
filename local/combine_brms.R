mega <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/local/mega_brms_results.rds")

library(data.table)

get_df_from_mega <- function(mega, nm) {
  res <- rbindlist(lapply(mega, function(el) {
    df <- el[[nm]]
    cbind(el$metadata, df)
  }))

  res[, Evid.Ratio := NULL]
  res[, Post.Prob := NULL]
  res[, Hypothesis := NULL]
  res[, rhs := NULL]
  res[, .filename := NULL]
  res[, Sensor := NULL] # this was in metadata and only reflect first row (wrong)
  setnames(res, "Group", "Sensor") # 'Group' is the emmeans word for the ranef blocking variable (sensor)

  return(res)
}

hilo_omission_by_sensor <- get_df_from_mega(mega, "hilo_omission_by_sensor")
hilo_reward_by_sensor <- get_df_from_mega(mega, "hilo_reward_by_sensor")
int_contrast_by_sensor <- get_df_from_mega(mega, "int_contrast_by_sensor")


fwrite(hilo_omission_by_sensor, file="hilo_omission_by_sensor.csv.gz")
fwrite(hilo_reward_by_sensor, file = "hilo_reward_by_sensor.csv.gz")
fwrite(int_contrast_by_sensor, file = "int_contrast_by_sensor.csv.gz")
