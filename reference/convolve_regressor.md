# This function convolves a regressor with a normalized HRF whose peak amplitude is 1.0

It extends `fmri.stimulus` by allowing for two normalization approaches
(building on AFNI dmUBLOCK): 1) "evtmax_1": pre-convolution HRF max=1.0
normalization of each stimulus regardless of duration: identical to
dmUBLOCK(1) 2) "durmax_1": pre-convolution HRF max=1.0 normalization for
long events (15+ sec) – height of HRF is modulated by duration of event:
identical to dmUBLOCK(0)

## Usage

``` r
convolve_regressor(
  n_vols,
  reg,
  tr = 1,
  normalization = "none",
  rm_zeros = TRUE,
  center_values = TRUE,
  convmax_1 = FALSE,
  demean_convolved = FALSE,
  high_pass = NULL,
  convolve = TRUE,
  ts_multiplier = NULL,
  hrf_parameters = c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35),
  microtime_resolution = 20L,
  lg = NULL
)
```

## Arguments

- n_vols:

  The number of volumes (scans) to be output in the convolved regressor

- reg:

  A matrix containing the trial, onset, duration, and value for each
  event

- tr:

  The repetition time in seconds

- normalization:

  The HRF normalization method used: "none", "durmax_1", or "evtmax_1"

- rm_zeros:

  Whether to remove zeros from events vector prior to convolution.
  Generally a good idea since we typically center values prior to
  convolution, and retaining zeros will lead them to be non-zero after
  mean centering.

- center_values:

  Whether to demean values vector before convolution. Default `TRUE`.

- convmax_1:

  Whether to rescale the convolved regressor to a maximum height of 1.

- demean_convolved:

  Whether to demean the regressor after convolution (default: `TRUE`)

- high_pass:

  The cutoff frequency (in Hz) used for high-pass filtering. If `NULL`,
  no filtering is applied.

- convolve:

  If `TRUE`, the regressor is convolved with the HRF. If `FALSE`, the
  regressor values are simply aligned onto the time grid without
  convolution based on the corresponding onsets, durations, and values.

- ts_multiplier:

  A vector that is n_vols in length that will be multiplied against the
  stimulus vector before convolution.

- hrf_parameters:

  A named vector of parameters passed to `fmri.stimulus` that control
  the shape of the double gamma HRF. Default:
  `c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)`.

- microtime_resolution:

  Number of microtime bins per TR used when evaluating and normalizing
  the convolved regressor.

- lg:

  An lgr object used for logging messages

## Author

Michael Hallquist
