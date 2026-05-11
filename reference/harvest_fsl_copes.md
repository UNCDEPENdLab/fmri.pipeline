# Compatibility wrapper for legacy internal callers

Compatibility wrapper for legacy internal callers

## Usage

``` r
harvest_fsl_copes(
  gpa,
  l1_model_names = NULL,
  l2_model_names = NULL,
  l3_model_names = NULL
)
```

## Arguments

- gpa:

  a glm_pipeline_arguments object

- l1_model_names:

  optional character vector of L1 models to include

- l2_model_names:

  optional character vector of L2 models to include

- l3_model_names:

  optional character vector of L3 models to process

## Value

a list of harvested FSL inputs for downstream L3 consumers
