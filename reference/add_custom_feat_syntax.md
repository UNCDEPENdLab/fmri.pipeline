# Helper function to find/replace or insert additional arguments to include in FSF syntax for a FEAT model.

Helper function to find/replace or insert additional arguments to
include in FSF syntax for a FEAT model.

## Usage

``` r
add_custom_feat_syntax(fsf_syntax, feat_args, lg = NULL)
```

## Arguments

- fsf_syntax:

  a character vector containing all fsf syntax

- feat_args:

  a named list

## Value

a modified version of `fsf_syntax` that now includes all custom
arguments from `feat_args`

## Details

If a tag already exists in the fsf syntax, it is replaced by the custom
argument in `feat_args`. If the tag does not exist, it is inserted.
