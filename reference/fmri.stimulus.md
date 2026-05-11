# Convolve a regressor with a hemodynamic response function for fMRI analysis.

Extended from `fmri` package to allow for continuous-valued regressor,
which is passed using the values parameter.

## Usage

``` r
fmri.stimulus(
  n_vols = 1,
  onsets = NULL,
  durations = NULL,
  values = NULL,
  units = "time",
  center_values = FALSE,
  rm_zeros = TRUE,
  convolve = TRUE,
  tr = 2,
  ts_multiplier = NULL,
  demean = TRUE,
  convmax_1 = FALSE,
  a1 = 6,
  a2 = 12,
  b1 = 0.9,
  b2 = 0.9,
  cc = 0.35,
  conv_method = "r",
  microtime_resolution = 20
)
```

## Arguments

- n_vols:

  The number of volumes (scans) to be output in the convolved regressor

- onsets:

  A vector of times (in volumes or seconds) specifying event onsets

- durations:

  A vector of durations (in volumes or seconds) for each event

- values:

  A vector of parametric values used as regressor heights prior to
  convolution

- units:

  Specifies the units of measurement for \`onsets\` and \`durations\`.
  If \`"seconds"\`, these are interpreted in terms of time in seconds.
  If \`"volumes"\`, onsets and durations are interpreted in terms of the
  volumes of the scan. Thus, an onset of \`1\` in volumes refers to time
  = 0s, and by extension, \`time = (volume - 1)\*tr\`. Default:
  \`"seconds"\`.

- center_values:

  Whether to demean values vector before convolution

- rm_zeros:

  Whether to remove zeros from events vector prior to convolution.
  Generally a good idea since we typically center values prior to
  convolution, and retaining zeros will lead them to be non-zero after
  mean centering.

- convolve:

  Whether to convolve the regressor with the HRF. If FALSE, a time
  series of events, durations, and heights is returned.

- tr:

  The repetition time (sometimes called TR) in seconds

- ts_multiplier:

  A time series of length `n_vols` that is multiplied against the
  stimulus before convolution.

- demean:

  Whether to demean the regressor after convolution. Default: TRUE

- convmax_1:

  Whether to rescale the convolved regressor to a maximum height of 1.

- a1:

  The a1 parameter of the double gamma

- a2:

  The a2 parameter of the double gamma

- b1:

  The b1 parameter of the double gamma

- b2:

  The b2 parameter of the double gamma

- cc:

  The cc parameter of the double gamma

- conv_method:

  Method for convolving HRF with stimulus. Either "r" or "cpp". The "r"
  method uses an FFT-based internal convolution with convolve(x, y,
  conj=TRUE). The "cpp" method uses an internal C++ function with a
  loop-based convolution over the vectors. Unfortunately, at present,
  the C++ approach is noticeably slower since it does not use FFT to
  obtain the filter.

- microtime_resolution:

  The number of bins between TRs used for calculating regressor values
  in continuous time

## Details

The function also supports mean centering of parametric regressor prior
to convolution to dissociate it from stimulus occurrence (when event
regressor also in model)
