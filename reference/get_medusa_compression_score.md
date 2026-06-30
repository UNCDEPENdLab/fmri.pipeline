# front-end function for taking a list of windowed time series by mask value, interpolating them onto a time grid, and (optionally) averaging across voxels/units within a value to derive the mean interpolated time series

front-end function for taking a list of windowed time series by mask
value, interpolating them onto a time grid, and (optionally) averaging
across voxels/units within a value to derive the mean interpolated time
series

## Usage

``` r
get_medusa_compression_score(
  fmri_obj,
  event = NULL,
  time_before = -3,
  time_after = 3,
  collide_before = NULL,
  collide_after = NULL,
  group_by = NULL
)
```

## Arguments

- fmri_obj:

  an fmri_ts object containing a single run of data with corresponding
  events

- event:

  the event to which the fmri time series should be aligned (column in
  `fmri_obj$event_data`)

- time_before:

  How many seconds before the `event` do we want data. Default: -3

- time_after:

  How many seconds after the `event` do we want data. Default: 3

- collide_before:

  An optional vector of column names in `trial_df` that set a boundary
  on the earliest time point used in interpolation. This effectively
  truncates the time series for interpolation to a smaller window than
  specified by `time_before`.

- collide_after:

  An optional vector of column names in `trial_df` that set a boundary
  on the latest time point used in interpolation. This effectively
  truncates the time series for interpolation to a smaller window than
  specified by `time_after`.

- group_by:

  return interpolated time series for each combination of group_by
  variables. Default is to provide one interpolated time series per
  trial.

## Author

Michael Hallquist
