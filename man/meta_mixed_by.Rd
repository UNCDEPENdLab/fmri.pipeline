% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_mixed_by.R
\name{meta_mixed_by}
\alias{meta_mixed_by}
\title{Function to run Bayesian random-effects meta-regression on coefficients from mixed_by}
\usage{
meta_mixed_by(
  coef_df,
  terms = "all",
  fit_subsets = "all",
  max_order = 3,
  outcome = NULL,
  fixef = NULL,
  rhs = NULL,
  brms_args = list(chains = 4, cores = 4, iter = 12000)
)
}
\arguments{
\item{coef_df}{A tidy data.frame of coefficients produced by mixed_by.}

\item{terms}{Which model terms to fit in meta-regression. If 'all', then all
unique values of the \code{term} column of \code{coef_df} will be fit in separate
meta-regression models. Alternatively, a character vector of terms can be specified
to run meta-regressions on a smaller subset}

\item{fit_subsets}{If \code{TRUE} or \code{"all"}, then all subsets of the split
variables will be tested in separate meta-regression models. This can be useful for
model comparison tests to examine whether a given factor has an overall effect. If
\code{FALSE} or \code{"none"}, no split variable subsets will be fit. If \code{"individual"},
then the split variables are added in subsets individually, but their interactions are not.}
}
\description{
Function to run Bayesian random-effects meta-regression on coefficients from mixed_by
}
