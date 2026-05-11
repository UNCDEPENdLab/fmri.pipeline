# helper function to guide user through process of choosing which models to run in GLM pipeline

helper function to guide user through process of choosing which models
to run in GLM pipeline

## Usage

``` r
choose_glm_set(
  gpa,
  l1_model_names = NULL,
  l2_model_names = NULL,
  l3_model_names = NULL,
  lg = NULL
)
```

## Arguments

- l1_model_names:

  a character vector of level 1 model names

- l2_model_names:

  a character vector of level 2 model names

- l3_model_names:

  a character vector of level 3 model names

## Value

a named list containing all models that were selected along with
additional information about whether to rerun existing models
