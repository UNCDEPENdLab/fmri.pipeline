# Submit SPM jobs to the cluster for unattended MATLAB/Octave execution

Submit SPM jobs to the cluster for unattended MATLAB/Octave execution

## Usage

``` r
run_spm_sepjobs(
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

  level of analysis (1 or 3)

- model_names:

  optional model names used to subset SPM runs

- rerun:

  logical indicating whether to rerun existing directories

- wait_for:

  optional parent job ids that should complete before these jobs
  commence

## Value

a vector of job ids for all scripts that were submitted
