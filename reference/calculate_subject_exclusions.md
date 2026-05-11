# small helper function to populate subject exclusions in \$run_data and \$subject_data

small helper function to populate subject exclusions in \$run_data and
\$subject_data

## Usage

``` r
calculate_subject_exclusions(gpa)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object that already has the `$exclude_run`
  column populated in `$l1_model_setup`.

## Value

a modified copy of `gpa` where \$exclude_subject has been added to
\$run_data and \$subject_data. This function also adds \$n_good_runs to
\$subject_data which is useful if we want to enforce a lower bound on
the number of runs used to define subject exclusion.

## Author

Michael Hallquist
