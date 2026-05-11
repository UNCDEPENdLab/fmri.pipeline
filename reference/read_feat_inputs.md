# This is a small helper function that returns the inputs provided in the feat_files field for a set of .fsf files.

This is a small helper function that returns the inputs provided in the
feat_files field for a set of .fsf files.

## Usage

``` r
read_feat_inputs(feat_list, recursive = FALSE)
```

## Arguments

- feat_list:

  a vector of design.fsf filenames, and/or .feat/.gfeat directory names

- recursive:

  a boolean indicating whether to drill down and find lower-level inputs
  for .gfeat/.feat inputs

## Details

One can also pass in .feat or .gfeat directories and the function will
use the design.fsf files within each of these to find the inputs

## Author

Michael Hallquist
