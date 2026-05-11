# Round and format the mean, min, and max as a single character string

Round and format the mean, min, and max as a single character string

## Usage

``` r
ol_mean_min_max(x, digits = 2, prespace = 0, showNA = "ifany", newline = TRUE)
```

## Arguments

- x:

  the numeric variable to summarize

- digits:

  the number of digits to use in rounding

- prespace:

  An integer indicating the number of spaces to include at the beginning
  of the string. Default: 0

- showNA:

  whether to print information about NAs. As with `table`, the options
  are "no", "ifany", and "always".

- newline:

  if TRUE, include a newline character `\n` at the end of the string

## Value

a formatted character string with mean, min, and max
