library(lme4)
library(dplyr)
library(brms)
library(foreach)
library(data.table)

#rds <- readRDS("/Users/michael/Downloads/meg_tf_rdf_wholebrain_rt_ri_single_sensor_RT1165")
rds <- readRDS("/Users/hallquist/OneDrive/collected_letters/papers/meg/plots/wholebrain/rt/output/meg_tf_rdf_wholebrain_rt_rs_single_sensor_RT1165")
names(rds$emtrends_list$emt_1)[16] <- "model_name_dummy"
library(ggplot2)

rds$emtrends_list$emt_1$Pow_fac <- as.factor(rds$emtrends_list$emt_1$Pow)
pdf("crazy_sensor.pdf", width=50, height=50)

ggplot(rds$emtrends_list$emt_1, aes(x=outcome...2, y=rt_csv_sc.trend, linetype=Pow_fac, group=Pow_fac)) +
  geom_errorbar(aes(ymin=rt_csv_sc.trend - std.error, ymax=rt_csv_sc.trend + std.error)) +
  geom_point() + geom_line(size=2) + facet_wrap(~Sensor)
dev.off()

pdf("crazy_sensor_smoothagg.pdf", width=5, height=5)

ggplot(rds$emtrends_list$emt_1, aes(x=outcome...2, y=rt_csv_sc.trend, linetype=Pow_fac, group=Pow_fac)) +
  stat_smooth(method="lm")
dev.off()


#lmer(rt_csv_sc.trend ~ Pow_fac) #outcome...2


library(dplyr)
library(brms)
library(emmeans)
ff2 <- brm(
  rt_csv_sc.trend | se(std.error, sigma=TRUE) ~ 1 + Pow_fac + (1 + Pow_fac | Sensor),
  rds$emtrends_list$emt_1 %>% dplyr::filter(outcome...2 == "Reward"), chains= 4, cores=4, iter = 15000,
  control = list(
    adapt_delta = 0.999,
    max_treedepth = 13
  )
)

post <- posterior_samples(ff2)

test <- emmeans(ff2, ~ Pow_fac)
pairs(test)


library(tidybayes)
xx <- contrast(test, "tukey")

vv <- gather_emmeans_draws(xx)

str(rds$emtrends_list)



ff3 <- brm(
  rt_csv_sc.trend | se(std.error, sigma=TRUE) ~ 1 + Pow_fac*outcome...2 + (1 + Pow_fac + outcome...2 | Sensor),
  rds$emtrends_list$emt_1, chains= 4, cores=4, iter = 15000,
  control = list(
    adapt_delta = 0.999,
    max_treedepth = 13
  )
)

emm <- emmeans(ff3, ~ Pow_fac | outcome...2)
p1 <- pairs(emm)
p2 <- pairs(p1)

x <- rds$emtrends_list$emt_1 %>%
  dplyr::select(Sensor, Pow_fac, rt_csv_sc.trend, std.error, outcome...2) %>%
  pivot_wider(values_from=c(rt_csv_sc.trend, std.error), names_from=c(Pow_fac), id_cols=c(Sensor, outcome...2)) %>%
  mutate(zhigh = rt_csv_sc.trend_2/std.error_2, zlow=`rt_csv_sc.trend_-2`/`std.error_-2`, zdiff=zhigh - zlow)


# 
# x <- rds$emtrends_list$emt_1 %>%
#   group_by(Pow_fac, outcome...2) %>%
#   mutate(gnum=cur_group_id()) %>%
#   ungroup() %>%
#   dplyr::select(gnum, Pow_fac, rt_csv_sc.trend, std.error, outcome...2) %>%
#   pivot_wider(values_from=c(rt_csv_sc.trend, std.error), names_from=c(Pow_fac), id_cols=c(rnum, outcome...2)) %>%
#   mutate(zpow = (rt_csv_sc.trend_2/std.error_2) - (`rt_csv_sc.trend_-2`/`std.error_-2`))
# 



lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)

library(rFSA)
totest <- ddf %>%
    filter(term == "v_max_wi") %>%
    select(estimate, side, evt_time, stream) %>%
    mutate(across(c(side, evt_time, stream), factor))

xx <- FSA(estimate ~ 1, data = totest, m = 2, interactions = TRUE)
abc <- lm(estimate ~ stream*side*evt_time, totest)


library(brms)
library(emmeans)
library(dplyr)
library(ggplot2)
library(data.table)
ddf <- ddf %>% mutate(
    time_fac = factor(evt_time)
)

ggplot(ddf %>% filter(term == "v_max_wi"), aes(x = time_fac, y = stream, fill = estimate)) +
    geom_tile()

ddf$splitid <- with(ddf, paste0(stream, time_fac, side))

fnull <- brm(
    estimate | se(std.error, sigma=FALSE) ~ 1 + (1 | splitid),
    ddf %>% filter(term == "v_max_wi"), chains= 4, cores=4, iter = 12000
)


ff <- brm(
    estimate | se(std.error, sigma=TRUE) ~ (stream + time_fac + side)^2 + (1 | splitid),
    ddf %>% filter(term == "v_max_wi"), chains= 4, cores=4, iter = 12000
)

ff <- add_criterion(ff, criterion = "loo")
fnull <- add_criterion(fnull, criterion = "loo")

loo_compare(fnull, ff) %>%
    print(simplify = F)


(mw <- model_weights(fnull, ff))

summary(ff)
test <- emmeans(ff, ~ stream | time_fac)

ff2 <- brm(
    estimate | se(std.error, sigma=TRUE) ~ 1 + (1 | stream) + 
    (1 | time_fac) + (1 | side) + (1 | stream:time_fac),
    ddf %>% filter(term == "v_max_wi"), chains= 4, cores=4, iter = 15000,
    control = list(
        adapt_delta = 0.999,
        max_treedepth = 13
    )
)


pairs(test)


library(tidybayes)
xx <- contrast(test, "tukey")

vv <- gather_emmeans_draws(xx)

pdf("brms_meta_demo_vmaxwi.pdf", width=12, height=8)
ggplot(
    vv,
    aes(y = time_fac, x = .value)
) +
    stat_halfeyeh() +
        facet_wrap(~contrast) +
        geom_vline(xintercept = 0, color = "red", lty = 2) +
        coord_flip()
dev.off()
