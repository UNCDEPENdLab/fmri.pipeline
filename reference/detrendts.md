# Detrend a time series up to quadratic trend. Used by fir1Bandpass prior to filtering

Detrend a time series up to quadratic trend. Used by fir1Bandpass prior
to filtering

## Usage

``` r
detrendts(x, order = 0)
```

## Arguments

- x:

  A time series to be detrended

- order:

  The polynomial order used for detrending. 0=demean; 1=linear;
  2=quadratic
