# Mixed by runs a set of mixed-effects models for each combination of a set of factors. Its primary use is to run the same model on different splits of the data.

Mixed by runs a set of mixed-effects models for each combination of a
set of factors. Its primary use is to run the same model on different
splits of the data.

## Usage

``` r
mixed_by(
  data,
  outcomes = NULL,
  rhs_model_formulae = NULL,
  model_formulae = NULL,
  split_on = NULL,
  external_df = NULL,
  external_merge_by = NULL,
  padjust_by = "term",
  padjust_method = "BY",
  outcome_transform = NULL,
  scale_predictors = NULL,
  ncores = 1L,
  cl = NULL,
  refit_on_nonconvergence = 3,
  tidy_args = list(effects = "fixed", conf.int = TRUE),
  lmer_control = lmerControl(optimizer = "nloptwrap"),
  calculate = c("parameter_estimates_reml", "parameter_estimates_ml", "fit_statistics"),
  return_models = FALSE,
  emmeans_spec = NULL,
  emtrends_spec = NULL
)
```

## Arguments

- data:

  A data.frame or data.table object containing stacked data for each
  combination of the `split_on` variables. The function will run
  separate mixed-effect models for each combination. Alternatively, a
  vector of filenames can be passed, which will be read in sequentially
  and fit (.rds, .csv, .dat, and .txt supported at present).

- outcomes:

  A character vector of outcome variables to be analyzed

- rhs_model_formulae:

  A named list of lme4-format formula specifying the exact model to be
  run for each data split.

- model_formulae:

  Alternative to the outcome + rhs_model_formulae approach. This is a
  list of lme4-format formulae that includes the outcome on the
  left-hand side. This is useful if the outcomes change from one model
  to the next, but you don't want the Cartesian product (combinations)
  of outcomes and rhs_model_formulae

- split_on:

  A character vector of columns in `data` used to split the analyses
  into separate models.

- external_df:

  An optional data.frame/data.table containing external data that should
  be joined with `data` prior to model fitting. Useful if `data`
  contains external time series data and `external_df` is a dataset of
  behavioral variables that do not vary by neural sensor/region.
  Optionally, this can be a single filename to a .rds file if you want
  to pass in the filename, not the data itself.

- external_merge_by:

  A character vector specifying which columns of `external_df` and
  `data` should be used for creating a combined dataset.

- padjust_by:

  A character vector or list consisting of one or more variables over
  which an adjusted p-value should be calculated. This defaults to
  adjusting by each term (fixed effect) in the model, which will adjust
  for all tests for a given term across variables in `split_on`.

- padjust_method:

  The adjustment method (see
  [`?p.adjust`](https://rdrr.io/r/stats/p.adjust.html)) for adjusting
  p-values. Multiple values can be passed as a character vector, in
  which case multiple corrections will be added as distinct columns.

- outcome_transform:

  A vectorized function that will be applied to the outcome variable
  prior to running the model. For example,
  `outcome_transform=function(x) { as.vector(scale(x)) } `.

- scale_predictors:

  An optional vector of predictor names that should be z-scored (unit
  normalized) prior to entering into the lmer model. Note that this
  scaling will be applied to the predictor regardless of which rhs model
  is being run as long as that predictor is in the model.

- ncores:

  The number of compute cores to be used in the computation. Defaults to
  1.

- cl:

  An optional external cl object (created by a variant of makeCluster)
  used for computation. Can save the overhead of starting and stopping
  many workers in a loop context.

- refit_on_nonconvergence:

  The number of times a model should be refit if it does not converge.
  Final estimates from one iteration are used as starting values for the
  next.

- tidy_args:

  A list of arguments passed to tidy.merMod for creating the coefficient
  data.frame. By default, the function only returns the fixed effects
  and computes confidence intervals using the Wald method.

- lmer_control:

  An lmerControl object specifying any optimization settings to be
  passed to lmer()

- calculate:

  A character vector specifying what calculations should be returned by
  the function. The options are: "parameter_estimates_reml",
  "parameter_estimates_ml", "fit_statistics", "fitted", and "residuals".

- return_models:

  A boolean indicating whether to return fitted model objects, which can
  be used for post hoc contrasts, visualization and statistics. Note
  that model objects can get very large, so be careful with this option
  since it could generate a massive data object.

- emmeans_spec:

  A named list of emmeans calls to be run for each model to obtain model
  predictions. Any arguments that are valid for emmeans can be passed
  through the list structure.

- emtrends_spec:

  A named list of emtrends calls to be run for each model to obtain
  model-predicted slopes. Any arguments that are valid for emtrends can
  be passed through the list structure.

## Value

A data.table object containing all coefficients for each model estimated
separately by `split_on`

## Details

In general, restricted maximum likelihood (REML) should be used for
making inferences about parameter estimates, whereas ML should be used
for model comparisons based on log-likelihood (e.g., AIC). If
"parameter_estimates_reml" are requested in `calculate`, then models
will be fitted with REML. If "parameter_estimates_ml" and/or
"fit_statistics" are requested, models will be fitted with ML. Note that
if both REML and ML are requested, each model is fit twice since the
estimators each have advantages and disadvantages noted above.

Example of emmeans_spec usage: mixed_by(data, emmeans_spec=list(
em1=list(specs = ~ memory \| noise_level, adjust = "sidak", weights =
"cells"), em2=list(specs = ~ memory \* noise_level, weights = "equal"),
em3=list(specs = ~ memory) ))
