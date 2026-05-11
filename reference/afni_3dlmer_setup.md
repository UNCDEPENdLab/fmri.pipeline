# Internal function to setup AFNI 3dLMEr longitudinal models

Internal function to setup AFNI 3dLMEr longitudinal models

## Usage

``` r
afni_3dlmer_setup(
  gpa,
  backend,
  lg,
  l1_model_names,
  l2_model_names,
  l3_model_names,
  l2_l3_pairs,
  subj_df,
  requires_l2
)
```

## Arguments

- gpa:

  a glm_pipeline_arguments object

- backend:

  the afni backend specification

- lg:

  the current logger

- l1_model_names:

  L1 models to process

- l2_model_names:

  L2 models to process

- l3_model_names:

  L3 models to process

- l2_l3_pairs:

  data.frame of compatible L2/L3 model pairs

- subj_df:

  filtered subject data.frame

- requires_l2:

  whether the resolved upstream producer emits required inputs at level
  2

## Value

a setup list with metadata and status data.frame
