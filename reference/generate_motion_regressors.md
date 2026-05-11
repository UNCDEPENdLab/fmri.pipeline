# helper function for generating motion regressors from raw 6-parameter motion coregistration

helper function for generating motion regressors from raw 6-parameter
motion coregistration

## Usage

``` r
generate_motion_regressors(
  motion_params_file = "motion.par",
  col.names = c("rx", "ry", "rz", "tx", "ty", "tz"),
  demean = FALSE,
  drop_volumes = 0L,
  last_volume = NULL,
  rot_units = "rad",
  tra_units = "mm",
  na.strings = "NA",
  lg = NULL
)
```

## Arguments

- motion_params_file:

  file containing 6 motion parameters

- col.names:

  names of columns in `motion_params_file`, in order from left to right

- drop_volumes:

  number of volumes to drop from beginning of motion params

- last_volume:

  final volume to include from motion params (if truncated at end). If
  `NULL`, then the end of the time series is not truncated.

- rot_units:

  The units of the rotation parameters. Default is "rad" for radians

- tra_units:

  The units of the translation parameter. Only support "mm" right now
  millimeters
