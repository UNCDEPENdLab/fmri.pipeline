# for example
# setwd("/Users/michael/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode/")
# load("/Users/michael/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode/rt_decode_output_streams_mixed_by.Rdata")
# setattr(ddf, "split_variables", c("stream", "evt_time", "side"))

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
#' @export
meta_mixed_by <- function(coef_df, terms = "all", fit_subsets = "all", max_order = 3, 
                          outcome=NULL, fixef= rhs=NULL,
                          brms_args = list(chains = 4, cores = 4, iter = 12000)) {
  
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package \"brms\" must be installed to use meta_mixed_by. Use install.packages('brms') first.",
      call. = FALSE
    )
  }

  checkmate::assert_data_frame(coef_df)
  checkmate::assert_character(terms)
  checkmate::assert_subset(outcome, names(coef_df))
  

  checkmate::assert_character(attr(coef_df, "split_on"))
  split_on <- attr(coef_df, "split_on")
  checkmate::assert_subset(split_on, names(coef_df)) # verify that all split variables are in coefs
  if (length(terms) == 1L && terms[1L] == "all") terms <- split_on

  if (is.logical(fit_subsets)) {
    coef_subset <- coef_df[term %in% terms, ]
  } # only fit models to requested terms
  split_coef <- split(coef_subset, by = terms)
  
  for (sub_df in split_coef) {
    
  }

  # build formula list
  if (isTRUE(fit_subsets)) {
    for (tt in seq_len(terms)) {

      combs <- combn(terms, m = tt) # elements x combinations matrix
      for (cc in combs) {
        # and loop over model order
      }
    }
  }

# 
#   for (ff in split_coef) {
# 
#   }
# 
# 
#   fnull <- brm(
#     estimate | se(std.error, sigma = FALSE) ~ 1 + (1 | splitid),
#     ddf %>% filter(term == "v_max_wi"),
#     chains = 4, cores = 4, iter = 12000
#   )
  
  
}

#generally support:
# draws on emmeans of brm objects
# stacking split coef_dfs, allowing for new things like estimate ~ type where type could be selective versus full etc.