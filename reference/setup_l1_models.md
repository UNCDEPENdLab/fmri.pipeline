# Function for setting up level 1 model specifications and corresponding files

Function for setting up level 1 model specifications and corresponding
files

## Usage

``` r
setup_l1_models(gpa, l1_model_names = NULL)
```

## Arguments

- gpa:

  A glm_model_arguments function setup by the `setup_glm_pipeline`
  function

- l1_model_names:

  A character vector of model names within `gpa$l1_models` whose l1
  inputs should be generated. If omitted, all models within
  `gpa$l1_models` will be generated.

## Details

The `l1_model_names` argument allows the creation of multiple l1 models
to be parallelized at a superordinate level or, even better, to be
spread across independent jobs on a cluster. This function already
provides the option to parallelize over subjects for a single model if
`gpa$l1_setup_cpus` is greater than 1.

## Author

Michael Hallquist
