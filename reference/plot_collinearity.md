# Create a visualization of collinearity diagnostics

Generates plots summarizing VIF values and correlations across subjects.

## Usage

``` r
plot_collinearity(collin_summary, plot_type = "both")
```

## Arguments

- collin_summary:

  An l1_collinearity_summary object from extract_l1_collinearity

- plot_type:

  Character string specifying the type of plot: "vif" for VIF boxplots,
  "correlation" for correlation heatmap, "both" for both plots. Default
  is "both".

## Value

A ggplot object or list of ggplot objects
