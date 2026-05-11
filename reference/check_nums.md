# helper function to validate numbers in an input string or vector

helper function to validate numbers in an input string or vector

## Usage

``` r
check_nums(inp, lower = 0, upper = 1e+10, as_string = TRUE)
```

## Arguments

- inp:

  A string, character vector, or numeric vector

- lower:

  The lowest valid value for each number

- upper:

  The highest valid value for each number

- as_string:

  If TRUE, return numbers as a space-separated string. If FALSE, return
  the validated numeric vector itself.
