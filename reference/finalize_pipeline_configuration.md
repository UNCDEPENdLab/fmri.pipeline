# This is a small helper function to validate the glm_model_arguments list structure. It adds a few details such as the output directory to make it less burdensome for to setup a pipeline N.B. gpa is a shorthand abbreviation for glm_model_arguments, to save typing

This is a small helper function to validate the glm_model_arguments list
structure. It adds a few details such as the output directory to make it
less burdensome for to setup a pipeline N.B. gpa is a shorthand
abbreviation for glm_model_arguments, to save typing

## Usage

``` r
finalize_pipeline_configuration(gpa, refinalize = FALSE)
```

## Arguments

- gpa:

  A `glm_pipeline_arguments` object setup by `setup_glm_pipeline`

- refinalize:

  A logical indicating whether to force checks and finalize steps on an
  object that was previously finalized.
