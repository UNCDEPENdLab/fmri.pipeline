# function to submit a set of jobs on a cluster to estimate many Feat level 1 models

function to submit a set of jobs on a cluster to estimate many Feat
level 1 models

## Usage

``` r
run_feat_sepjobs(
  gpa,
  level = 1L,
  model_names = NULL,
  rerun = FALSE,
  wait_for = ""
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object containing model specification

- level:

  GLM level to run.

- model_names:

  an optional level 1 model name used to subset the runs to submit to
  feat. If not provided, all level 1 models in `gpa` will be submitted
  for feat estimation.

- rerun:

  a logical indicating whether to re-run an existing directory. Default:
  FALSE

- wait_for:

  a parent job that should complete before these jobs commence

## Value

a vector of job ids for all files that were submitted to the cluster
