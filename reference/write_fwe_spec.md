# Write or read an FWE correction specification as YAML

Write or read an FWE correction specification as YAML

## Usage

``` r
write_fwe_spec(x, file, overwrite = FALSE)

read_fwe_spec(file)
```

## Arguments

- x:

  an `fwe_spec` or equivalent list.

- file:

  YAML file path.

- overwrite:

  whether an existing file may be replaced.

## Value

`write_fwe_spec()` invisibly returns the normalized output path;
`read_fwe_spec()` returns an `fwe_spec`.
