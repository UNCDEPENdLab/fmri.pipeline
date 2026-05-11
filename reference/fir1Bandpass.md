# Apply a FIR-1 bandpass filter a signal. Can low- or high-pass filter by specifying 0 for low or \>= Nyquist for high.

Apply a FIR-1 bandpass filter a signal. Can low- or high-pass filter by
specifying 0 for low or \>= Nyquist for high.

## Usage

``` r
fir1Bandpass(
  x,
  TR = 2,
  low = 0.009,
  high = 0.08,
  n = 250,
  plotFilter = FALSE,
  forward_reverse = TRUE,
  padx = 0,
  detrend = 1
)
```

## Arguments

- x:

  The time series to be filtered

- TR:

  The sampling frequency in seconds

- low:

  The lower filter cutoff in Hz. Fluctuations below this frequency will
  be filtered out

- high:

  The upper filter cutoff in Hz. Fluctuations above this frequency will
  be filtered out

- n:

  The order of the filter coefficients. Should probably leave this alone
  in general
