# Compile an FWE specification into an executable plan

Resolve an FWE specification against completed GLM outputs and compile
it through a method adapter. Planning is side-effect free: output
directories and result files are described but not created. The pTFCE
adapter creates one task per selected z-statistic while retaining shared
GFEAT-level inputs in an analysis table.

## Usage

``` r
plan_fwe_correction(
  gpa,
  spec,
  output_directory = NULL,
  require_inputs = TRUE,
  source = c("auto", "setup", "cache", "filesystem"),
  cache_dir = NULL,
  refresh_status = FALSE,
  lg = NULL
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` containing group-analysis output metadata.

- spec:

  an `fwe_spec` or equivalent list.

- output_directory:

  root directory for this correction. The default is
  `<gpa$output_directory>/fwe/<spec$name>`.

- require_inputs:

  if `TRUE`, error when a required statistic, mask, or FSL smoothness
  file is missing. If `FALSE`, return a plan with blocked tasks for
  inspection.

- source, cache_dir, refresh_status:

  passed to
  [`resolve_fwe_targets()`](https://hallquistlab.github.io/fmri.pipeline/reference/resolve_fwe_targets.md).

- lg:

  optional
  [`lgr::Logger`](https://s-fleck.github.io/lgr/reference/Logger.html).

## Value

an object of class `fwe_plan` containing `analyses`, `tasks`, and
expected `artifacts` tables.
