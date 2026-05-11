# Print the top n unique values and frequencies of an atomic variable

Print the top n unique values and frequencies of an atomic variable

## Usage

``` r
ol_unique(
  x,
  max_show = 5,
  max_allowed = 25,
  prespace = 0,
  showNA = "ifany",
  newline = TRUE
)
```

## Arguments

- x:

  the variable to tabulate

- max_show:

  the maximum number of unique values to show

- max_allowed:

  the maximum number of allowed unique values before the function
  summarizes as mean, min, max

- prespace:

  An integer indicating the number of spaces to include at the beginning
  of the string. Default: 0

- showNA:

  whether to print information about NAs. As with `table`, the options
  are "no", "ifany", and "always".

- newline:

  if TRUE, include a newline character `\n` at the end of the string

## Value

a formatted character string with frequencies for the top `max_show`
categories
