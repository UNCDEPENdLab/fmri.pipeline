# Estimate a level 2 (subject) model using FSL FEAT with fixed effects integration of runs

Estimate a level 2 (subject) model using FSL FEAT with fixed effects
integration of runs

## Usage

``` r
is_l2_passthrough_estimable(dmat, cmat, tol = sqrt(.Machine$double.eps))
```

## Arguments

- l1_df:

  a data.frame containing all runs for a single subject and a single l1
  model. This data.frame defines the inputs for the L2 analysis (i.e.,
  which runs to combine).

- l2_model:

  a model string in gpa\$l2_models containing the L2 model to setup

- gpa:

  a `glm_pipeline_arguments` object containing model specification

## Author

Michael Hallquist
