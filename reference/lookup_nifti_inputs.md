# helper function to locate nifti inputs for each run of data

helper function to locate nifti inputs for each run of data

## Usage

``` r
lookup_nifti_inputs(gpa, add_run_volumes = TRUE, add_dim = TRUE)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object created by setup_glm_pipeline

- add_run_volumes:

  whether to lookup and add the number of volumes in each run_nifti in
  the run_volumes column

- add_dim:

  whether to lookup and add the dimensions of each run_nifti in the
  dim_x, dim_y, and dim_z columns. This will also add an `nvoxels`
  column used in FEAT fsfs, which is x \* y \* z \* t.

## Value

a modified version of the `gpa` object that contains a \$run_nifti field
and relevant attribute columns (e.g., dim_x) the \$run_data data.frame.

## Author

Michael Hallquist
