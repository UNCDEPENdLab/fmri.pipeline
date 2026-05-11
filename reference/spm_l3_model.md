# Setup a third-level SPM model using lower-level con images

Setup a third-level SPM model using lower-level con images

## Usage

``` r
spm_l3_model(l3_df = NULL, gpa, execute_spm = FALSE, model_type = NULL)
```

## Arguments

- l3_df:

  a data.frame containing one row per subject with fields including id,
  session, l1_model, l3_model, l1_cope_name, and con_file

- gpa:

  a gpa (glm_pipeline_arguments) object containing model specification

- execute_spm:

  whether to execute SPM setup/estimation immediately. Default: FALSE

- model_type:

  optional override for SPM L3 model type. Options:
  "flexible_factorial", "mreg", "one_sample"

## Author

Michael Hallquist
