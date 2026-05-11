# Visualize design matrix, including event onset times and run boundaries

Visualize design matrix, including event onset times and run boundaries

## Usage

``` r
visualize_design_matrix(
  d,
  outfile = NULL,
  run_boundaries = NULL,
  events = NULL,
  include_baseline = FALSE
)
```

## Arguments

- d:

  a concatenated design matrix created by `build_design_matrix` and
  passed to `concat_design_runs`

- outfile:

  a filename used to export the design matrix using `ggsave`

- run_boundaries:

  a named vector of positions in the time series where a run boundary
  occurs (used for plotting)

- events:

  a named list of vectors containing the times of each event

- include_baseline:

  whether to display the baseline regressors in the design.
