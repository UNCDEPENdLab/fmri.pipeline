# Read a collected FWE result or artifact manifest

Read a collected FWE result or artifact manifest

## Usage

``` r
read_fwe_result(path)

read_fwe_artifact_manifest(path)
```

## Arguments

- path:

  an RDS/CSV file or correction output directory.

## Value

`read_fwe_result()` returns an `fwe_result`;
`read_fwe_artifact_manifest()` returns a validated data frame.
