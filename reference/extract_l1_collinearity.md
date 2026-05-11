# Extract collinearity diagnostics from all L1 models

This function reads the saved build_design_matrix objects (d_obj) from
all subjects and models, extracting collinearity diagnostics
(correlations and VIFs) into a structured format suitable for review and
quality control.

## Usage

``` r
extract_l1_collinearity(
  gpa,
  l1_model_names = NULL,
  vif_threshold = 5,
  cor_threshold = 0.8
)
```

## Arguments

- gpa:

  A glm_pipeline_arguments object that has been processed through
  setup_l1_models

- l1_model_names:

  Character vector of L1 model names to extract. If NULL (default),
  extracts from all available models.

- vif_threshold:

  Numeric threshold for flagging high VIF values. Default is 5.

- cor_threshold:

  Numeric threshold for flagging high correlations. Default is 0.8.

## Value

A list with class "l1_collinearity_summary" containing:

- `vif_summary`: A data.frame with VIF values for each regressor,
  subject, session, run, and model. Includes a flag for VIFs exceeding
  the threshold.

- `correlation_summary`: A data.frame with pairwise correlations between
  regressors for each subject, session, run, and model. Includes a flag
  for correlations exceeding the threshold.

- `high_vif`: A data.frame subset of vif_summary where VIF exceeds the
  threshold.

- `high_correlation`: A data.frame subset of correlation_summary where
  \|r\| exceeds the threshold.

- `summary_stats`: Aggregate statistics across all subjects/runs.

- `thresholds`: The VIF and correlation thresholds used.

## Details

This function iterates through all subject/model combinations in the gpa
object, loading the cached design matrix objects and extracting
collinearity information. The resulting summaries can be used to
identify problematic regressors or runs that may have estimation issues
due to multicollinearity.

Common causes of high collinearity include:

- Parametric regressors that are highly correlated with task indicator
  regressors

- Events that occur too close together in time

- Insufficient variation in parametric modulators within a run

## Examples

``` r
if (FALSE) { # \dontrun{
  # After running setup_l1_models
  collin_summary <- extract_l1_collinearity(gpa)
  
  # View high VIF regressors
  print(collin_summary$high_vif)
  
  # View highly correlated regressor pairs
  print(collin_summary$high_correlation)
  
  # Get overall summary
  print(collin_summary)
} # }
```
