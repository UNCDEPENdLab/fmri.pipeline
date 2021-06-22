
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