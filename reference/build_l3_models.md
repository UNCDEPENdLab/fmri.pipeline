# Interactive function to build an l3 model specification for setup_glm_pipeline

Interactive function to build an l3 model specification for
setup_glm_pipeline

## Usage

``` r
build_l3_models(gpa, from_spec_file = NULL, regressor_cols = NULL)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object.

- from_spec_file:

  optional YAML or JSON file containing L3 model settings.

- regressor_cols:

  optional character vector of columns to consider as regressors.
