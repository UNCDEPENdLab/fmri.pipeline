# Write or read a compiled FWE plan

Write or read a compiled FWE plan

## Usage

``` r
write_fwe_plan(x, file, overwrite = FALSE)

read_fwe_plan(file)
```

## Arguments

- x:

  an `fwe_plan`.

- file:

  an `.rds` file path.

- overwrite:

  whether an existing file may be replaced.

## Value

`write_fwe_plan()` invisibly returns the normalized file path;
`read_fwe_plan()` returns a validated `fwe_plan`.
