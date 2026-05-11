# Summarize VIF by regressor across all subjects

Creates a summary table showing VIF statistics for each regressor across
all subjects and runs.

## Usage

``` r
summarize_vif_by_regressor(collin_summary, by_model = FALSE)
```

## Arguments

- collin_summary:

  An l1_collinearity_summary object from extract_l1_collinearity

- by_model:

  Logical. If TRUE, summarize separately by L1 model. Default FALSE.

## Value

A data.frame with columns: regressor, mean_vif, sd_vif, min_vif,
max_vif, n_runs, n_flagged, pct_flagged
