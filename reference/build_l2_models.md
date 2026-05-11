# Interactive function to build an l2 model specification for setup_glm_pipeline

Interactive function to build an l2 model specification for
setup_glm_pipeline

## Usage

``` r
build_l2_models(gpa, from_spec_file = NULL, regressor_cols = NULL)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object containing an analysis pipeline to
  which \$l2_models should be added. If \$l2_models is already present,
  these will be amended.

- from_spec_file:

  optional YAML or JSON file containing settings to populated into l2
  models

- regressor_cols:

  an optional character vector of columns in `data` that should be
  considered as possible regressors

## Value

a `hi_model_set` object containing a list of higher-level regression
models

## Author

Michael Hallquist
