# Resolve an FWE specification to concrete group-analysis targets

Resolve semantic model and contrast selectors against the canonical FSL
output inventory returned by
[`lookup_feat_outputs()`](https://uncdependlab.github.io/fmri.pipeline/reference/lookup_feat_outputs.md).
Resolution never constructs FEAT paths independently.

## Usage

``` r
resolve_fwe_targets(
  gpa,
  spec,
  require_existing = TRUE,
  require_complete = TRUE,
  source = c("auto", "setup", "cache", "filesystem"),
  cache_dir = NULL,
  refresh_status = FALSE,
  lg = NULL
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` with level-3 FSL setup metadata.

- spec:

  an `fwe_spec` or equivalent list.

- require_existing:

  if `TRUE`, fail when a selected statistic image is absent. If `FALSE`,
  unresolved output paths are retained for planning and preview.

- require_complete:

  if `TRUE`, fail when a selected group analysis is not complete.

- source, cache_dir, refresh_status:

  passed to
  [`lookup_feat_outputs()`](https://uncdependlab.github.io/fmri.pipeline/reference/lookup_feat_outputs.md).

- lg:

  optional
  [`lgr::Logger`](https://s-fleck.github.io/lgr/reference/Logger.html).

## Value

a data.frame of class `fwe_target_set`, with one row per selected group
contrast and stable `target_key` and `target_id` columns.
