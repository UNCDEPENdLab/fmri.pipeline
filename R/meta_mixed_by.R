# meta_mixed_by <- function(coef_df, )

setwd("/Users/michael/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode/")

#for example
load("/Users/michael/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode/rt_decode_output_streams_mixed_by.Rdata")
setattr(ddf, "split_variables", c("stream", "evt_time", "side"))

#' Function to run Bayesian random-effects meta-regression on coefficients from mixed_by
#' 
#' @param coef_df A tidy data.frame of coefficients produced by mixed_by.
#' @param terms Which model terms to fit in meta-regression. If 'all', then all
#'   unique values of the \code{term} column of \code{coef_df} will be fit in separate
#'   meta-regression models. Alternatively, a character vector of terms can be specified
#'   to run meta-regressions on a smaller subset
#' @param fit_subsets If \code{TRUE} or \code{"all"}, then all subsets of the split
#'   variables will be tested in separate meta-regression models. This can be useful for
#'   model comparison tests to examine whether a given factor has an overall effect. If
#'   \code{FALSE} or \code{"none"}, no split variable subsets will be fit. If \code{"individual"},
#'   then the split variables are added in subsets individually, but their interactions are not.
#'   
#' @importFrom brms brm
meta_mixed_by <- function(coef_df, terms = "all", fit_subsets = "all", max_order=3,
                          brms_args = list(chains = 4, cores = 4, iter = 12000)) {
    checkmate::assert_data_frame(coef_df)
    checkmate::assert_character(terms)

    checkmate::assert_character(attr(coef_df, "split_on"))
    split_on <- attr(coef_df, "split_on")
    checkmate::assert_subset(split_on, names(coef_df)) #verify that all split variables are in coefs
    if (length(terms) == 1L && terms[1L] == "all") terms <- split_on

    if (is.logical(fit_subsets))

    coef_subset <- coef_df[term %in% terms, ] # only fit models to requested terms
    split_coef <- split(coef_subset, by = terms)

    # build formula list
    if (isTRUE(fit_subsets)) {
        for (tt in seq_len(terms)) {
            combs <- gtools::combinations(n = length(terms), r = tt, v = terms)
            for (cc in combs) {
                #and loop over model order
            }
        }
    }


    for (ff in split_coef) {

    }
    

fnull <- brm(
    estimate | se(std.error, sigma = FALSE) ~ 1 + (1 | splitid),
    ddf %>% filter(term == "v_max_wi"),
    chains = 4, cores = 4, iter = 12000
)


}

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