# Harvest standardized L3 input tables from upstream producer backends

Harvest standardized L3 input tables from upstream producer backends

## Usage

``` r
harvest_l3_inputs(
  gpa,
  l3_backend,
  l1_model_names = NULL,
  l2_model_names = NULL,
  l3_model_names = NULL,
  lg = NULL
)
```

## Arguments

- gpa:

  a glm_pipeline_arguments object

- l3_backend:

  target L3 backend specification or backend name

- l1_model_names:

  optional character vector of L1 models to include

- l2_model_names:

  optional character vector of L2 models to include

- l3_model_names:

  optional character vector of L3 models to process

- lg:

  optional logger

## Value

a list of data.frames, one per L1/L2 contrast combination, for each L3
model
