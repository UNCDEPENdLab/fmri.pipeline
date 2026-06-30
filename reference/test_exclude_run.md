# calculate whether to retain or exclude this run

calculate whether to retain or exclude this run

## Usage

``` r
test_exclude_run(
  gpa,
  id,
  session,
  run_number,
  generate_run_exclusion,
  all_confounds,
  lg,
  confound_path = NA_character_,
  motion_path = NA_character_
)
```

## Arguments

- gpa:

  A `glm_pipeline_arguments` object.

- id:

  subject id.

- session:

  session number.

- run_number:

  run number.

- generate_run_exclusion:

  whether exclusion should be evaluated.

- all_confounds:

  Combined confound data.

- lg:

  Logger object.

- confound_path:

  Path to the raw confound file.

- motion_path:

  Path to the raw motion file.
