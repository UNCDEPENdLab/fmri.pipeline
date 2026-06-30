# function to submit 3dLMEr jobs on a cluster

function to submit 3dLMEr jobs on a cluster

## Usage

``` r
run_3dlmer_sepjobs(
  gpa,
  level = 3L,
  model_names = NULL,
  rerun = FALSE,
  wait_for = ""
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object

- level:

  GLM level to run; currently only level 3 is supported.

- model_names:

  optional subset of L3 model names to run

- rerun:

  logical, whether to rerun existing models

- wait_for:

  job ID(s) to wait for

## Value

vector of submitted job IDs
