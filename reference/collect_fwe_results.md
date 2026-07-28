# Collect FWE results and persist the artifact manifest

Refresh expected artifacts, join them to task and semantic target
metadata, and optionally write an atomic CSV manifest plus an RDS result
snapshot. Repeated collection updates the same canonical files, making
this suitable for collecting scheduler results after jobs finish.

## Usage

``` r
collect_fwe_results(
  x,
  output_directory = NULL,
  require_complete = FALSE,
  persist = TRUE
)
```

## Arguments

- x:

  an `fwe_plan`, `fwe_run`, `fwe_result`, saved RDS result, or result
  directory.

- output_directory:

  correction output directory. Defaults to the plan's output directory.

- require_complete:

  if `TRUE`, fail unless every required artifact exists.

- persist:

  if `TRUE`, atomically write `manifest/artifacts.csv` and
  `manifest/result.rds`.

## Value

an `fwe_result` containing the refreshed plan and normalized artifact
manifest.
