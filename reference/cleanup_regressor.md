# Small utility function to handle any NAs in the regressor specification. Also, optionally, removes values of 0 from the regressor prior to convolution. The latter is usually a good idea to avoid a non-event (zero) becoming an event after mean centering of a parametric regressor. Convolving a zero-height regressor with an HRF will maintain all zeros, whereas a non-zero value (resulting from mean centering) will not.

Small utility function to handle any NAs in the regressor specification.
Also, optionally, removes values of 0 from the regressor prior to
convolution. The latter is usually a good idea to avoid a non-event
(zero) becoming an event after mean centering of a parametric regressor.
Convolving a zero-height regressor with an HRF will maintain all zeros,
whereas a non-zero value (resulting from mean centering) will not.

## Usage

``` r
cleanup_regressor(times, durations, values, rm_zeros = TRUE)
```

## Arguments

- times:

  Vector of times (in seconds) for regressor events

- durations:

  Vector of durations (in seconds) for regressor events

- values:

  Vector of values (parametric heights) for regressor events

- rm_zeros:

  Logical indicating whether to remove zero-valued events from regressor

## Author

Michael Hallquist
