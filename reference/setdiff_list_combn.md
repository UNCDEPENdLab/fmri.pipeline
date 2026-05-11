# internal function to print bidirectional set differences in a list containing vectors to be compared

internal function to print bidirectional set differences in a list
containing vectors to be compared

## Usage

``` r
setdiff_list_combn(l)
```

## Arguments

- l:

  a named list containing vectors to be compared

## Value

NULL (invisibly)

## Details

All pairwise combinations of vectors in the list will be compared. The
setdiff() operation is run twice for each pair, reflecting the two
directions of comparison (s1
