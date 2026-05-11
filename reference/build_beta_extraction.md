# helper function to prompt user to choose combinations of l1, l2, and l3 models for beta extraction.

helper function to prompt user to choose combinations of l1, l2, and l3
models for beta extraction.

## Usage

``` r
build_beta_extraction(
  gpa,
  extract_l1 = "prompt",
  extract_l2 = "prompt",
  extract_l3 = "prompt",
  lg = NULL
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object

- lg:

  a logger object for messages

## Value

a list of data frames (l1, l2, l3) consisting of the chosen combinations
of l1, l2, and l3 models
