# Function to convert dmat (runs x regressor list) to a time-oriented representation. This yields a list of runs where each element is data.frame of volumes x regressors

Function to convert dmat (runs x regressor list) to a time-oriented
representation. This yields a list of runs where each element is
data.frame of volumes x regressors

## Usage

``` r
place_dmat_on_time_grid(
  dmat,
  convolve = TRUE,
  run_timing = NULL,
  bdm_args,
  lg = NULL
)
```

## Arguments

- dmat:

  A runs x regressors 2-d list where each element is a matrix containing
  onset, duration, and value for a signal

- convolve:

  If `TRUE` (default), convolve the time-oriented signals with an HRF

- run_timing:

  A vector of cumulative start times for each run in a multi-run dataset

- bdm_args:

  A list of arguments passed to build_design_matrix, as well as a few
  fields added during the initial argument parsing. See
  build_design_matrix for details. Should contain:

  - convolve_wi_run TRUE/FALSE

  - run_volumes Numeric vector of run length

  - normalizations Character vector of HRF normalizations for each
    regressor. Options are "none", "durmax_1", or "evtmax_1".

  - add_derivs A logical vector (`TRUE/FALSE`) of regressors whose
    temporal derivatives should be included. Temporal derivatives are
    only applied if `convolve` is `TRUE`.

  - convmax_1 A logical vector (`TRUE/FALSE`) denoting whether to
    rescale max height to 1 after convolution

  - high_pass The cutoff frequency (in Hz) used for high-pass filtering.
    If `NULL`, no filtering is applied.

  - tr The repetition time (sometimes called TR) in seconds

  - hrf_parameters The parameters for the double-gamma HRF

- lg:

  An lgr logger object used for logging messages

## Details

Note that any volumes dropped from the beginning of each run should
already be reflected in the timings of regressors in `dmat`. This
prevents us from needing to have a drop_volumes implementation inside
convolve_regressor, which is confusing anyhow. Likewise, run_timing
should reflect the post-drop cumulative volumes.

## Author

Michael Hallquist
