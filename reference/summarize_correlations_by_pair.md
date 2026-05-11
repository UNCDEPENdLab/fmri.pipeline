# Summarize correlations by regressor pair across all subjects

Creates a summary table showing correlation statistics for each
regressor pair across all subjects and runs.

## Usage

``` r
summarize_correlations_by_pair(collin_summary, by_model = FALSE)
```

## Arguments

- collin_summary:

  An l1_collinearity_summary object from extract_l1_collinearity

- by_model:

  Logical. If TRUE, summarize separately by L1 model. Default FALSE.

## Value

A data.frame with columns: regressor1, regressor2, mean_cor, sd_cor,
min_cor, max_cor, mean_abs_cor, n_runs, n_flagged, pct_flagged
