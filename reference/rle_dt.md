# R6 class for a keyed data.table object that uses run-length encoding to reduce storage size

R6 class for a keyed data.table object that uses run-length encoding to
reduce storage size

R6 class for a keyed data.table object that uses run-length encoding to
reduce storage size

## Methods

### Public methods

- [`rle_dt$new()`](#method-rle_dt-new)

- [`rle_dt$get()`](#method-rle_dt-get)

- [`rle_dt$clone()`](#method-rle_dt-clone)

------------------------------------------------------------------------

### Method `new()`

Create an RLE-encoded copy of the data

#### Usage

    rle_dt$new(data, keys = NULL, optimize_order = (length(keys) <= 4))

#### Arguments

- `data`:

  A data.frame or data.table object containing original data

- `keys`:

  A character vector of keys within `data` that should be RLE-encoded

- `optimize_order`:

  A boolean indicating whether to search for the smallest encoding
  scheme. The default it to search in data with 4 or fewer keys, but not
  to search for 5+. Optionally, if a positive integer, randomly test
  ordering over this number of permutations.

------------------------------------------------------------------------

### Method [`get()`](https://rdrr.io/r/base/get.html)

Simple method to return the data.table with all columns in their
original form.

#### Usage

    rle_dt$get()

#### Details

Note that the data are modified slightly in that the keys columns are
placed first, and the data are ordered in the order of the keys (as
originally provided, left-to-right)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    rle_dt$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
