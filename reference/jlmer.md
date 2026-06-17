# Fit a linear mixed model with Julia MixedModels

`jlmer()` is an
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)-style wrapper
around Julia's `MixedModels.fit()`. It fits the model in Julia, then
converts the result back to an R `lmerMod` object through JellyMe4 so
downstream R tooling such as `broom.mixed`, `emmeans`, and `lme4`
methods can be used at the call site.

## Usage

``` r
jlmer(
  formula,
  data,
  REML = TRUE,
  JULIA_HOME = NULL,
  na.action = stats::na.omit,
  lmer_test = TRUE,
  lmer_test_tol = 1e-08,
  julia_setup_args = list(),
  julia_packages = c("MixedModels", "RCall", "JellyMe4"),
  fit_args = list(),
  keep_julia_model = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- formula:

  A two-sided mixed model formula.

- data:

  A data frame containing variables used in `formula`.

- REML:

  Logical scalar. If `TRUE`, estimate with restricted maximum
  likelihood; otherwise use maximum likelihood.

- JULIA_HOME:

  Optional path to the Julia installation. Passed to
  [`JuliaCall::julia_setup()`](https://rdrr.io/pkg/JuliaCall/man/julia_setup.html).

- na.action:

  Optional missing-data action for formula variables. Supports
  [`stats::na.omit`](https://rdrr.io/r/stats/na.fail.html),
  [`stats::na.exclude`](https://rdrr.io/r/stats/na.fail.html),
  [`stats::na.fail`](https://rdrr.io/r/stats/na.fail.html), or `NULL`.
  Omitted or excluded rows are removed before sending data to Julia.

- lmer_test:

  Logical scalar. If `TRUE`, convert the returned `lmerMod` object to
  `lmerTest`'s `lmerModLmerTest` class. This matches the default
  `mixed_by(engine = "lme4")` call site, which uses
  [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html).

- lmer_test_tol:

  Numeric tolerance passed to
  [`lmerTest::as_lmerModLmerTest()`](https://rdrr.io/pkg/lmerTest/man/as_lmerModLmerTest.html)
  when `lmer_test = TRUE`.

- julia_setup_args:

  Named list of additional arguments passed to
  [`JuliaCall::julia_setup()`](https://rdrr.io/pkg/JuliaCall/man/julia_setup.html).

- julia_packages:

  Character vector of Julia packages to load. Defaults to
  `c("MixedModels", "RCall", "JellyMe4")`.

- fit_args:

  Named list of additional keyword arguments passed to Julia's
  `fit(MixedModel, ...)` call.

- keep_julia_model:

  Logical scalar. If `TRUE`, the temporary Julia model variable is left
  in Julia's `Main` module and its name is stored in
  `attr(result, "jlmer")`. By default all temporary Julia variables are
  cleared after conversion to R.

- verbose:

  Logical scalar. If `TRUE`, print the Julia fit call.

- ...:

  Additional named keyword arguments passed to Julia's
  `fit(MixedModel, ...)` call. These are combined with `fit_args`.

## Value

An `lmerMod` object converted by JellyMe4, with a `"jlmer"` attribute
containing backend metadata. By default, the object is further converted
to `lmerTest`'s `lmerModLmerTest` class.

## Details

This function requires Julia packages `MixedModels`, `RCall`, and
`JellyMe4` in the active Julia environment. It does not install Julia
packages automatically.
