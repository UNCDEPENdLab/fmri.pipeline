# compute a moving average smooth over a time series (here, a vector of RTs)

used to fit smoothed RTs (`clock_model` object). Function is an adapted
version of `runmean` from the `caTools` package.

## Usage

``` r
runmean(x, k = 5)
```
