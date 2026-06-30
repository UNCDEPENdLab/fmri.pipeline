# Function to run Bayesian random-effects meta-regression on coefficients from mixed_by

Function to run Bayesian random-effects meta-regression on coefficients
from mixed_by

## Usage

``` r
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
```

## Arguments

- coef_df:

  A tidy data.frame of coefficients produced by mixed_by.

- terms:

  Which model terms to fit in meta-regression. If 'all', then all unique
  values of the `term` column of `coef_df` will be fit in separate
  meta-regression models. Alternatively, a character vector of terms can
  be specified to run meta-regressions on a smaller subset

- fit_subsets:

  If `TRUE` or `"all"`, then all subsets of the split variables will be
  tested in separate meta-regression models. This can be useful for
  model comparison tests to examine whether a given factor has an
  overall effect. If `FALSE` or `"none"`, no split variable subsets will
  be fit. If `"individual"`, then the split variables are added in
  subsets individually, but their interactions are not.

- max_order:

  Maximum subset order for meta-regression combinations.

- outcome:

  Optional outcome column override.

- fixef:

  Optional fixed-effect specification.

- rhs:

  Optional right-hand-side formula specification.

- brms_args:

  Arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
