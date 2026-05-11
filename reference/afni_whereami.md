# wrapper class for AFNI whereami

wrapper class for AFNI whereami

wrapper class for AFNI whereami

## Methods

### Public methods

- [`afni_whereami$new()`](#method-afni_whereami-new)

- [`afni_whereami$get_call()`](#method-afni_whereami-get_call)

- [`afni_whereami$get_omask_call()`](#method-afni_whereami-get_omask_call)

- [`afni_whereami$get_whereami_df()`](#method-afni_whereami-get_whereami_df)

- [`afni_whereami$get_overlap_df()`](#method-afni_whereami-get_overlap_df)

- [`afni_whereami$get_omask_df()`](#method-afni_whereami-get_omask_df)

- [`afni_whereami$get_outputs()`](#method-afni_whereami-get_outputs)

- [`afni_whereami$run()`](#method-afni_whereami-run)

- [`afni_whereami$clone()`](#method-afni_whereami-clone)

------------------------------------------------------------------------

### Method `new()`

create new afni_whereami object

#### Usage

    afni_whereami$new(
      afni_3dclusterize_obj = NULL,
      omask = NULL,
      coord_file = NULL,
      coord_file_columns = NULL,
      coord_vector = NULL,
      atlases = NULL,
      coord_orientation = NULL,
      coord_space = NULL,
      output_file = NULL,
      omask_output_file = NULL,
      afnidir = NULL
    )

#### Arguments

- `afni_3dclusterize_obj`:

  Use an existing 3dClusterize object to build the whereami lookup

- `omask`:

  Use this integer-valued cluster mask to lookup percent overlap with
  different atlas parcels

- `coord_file`:

  The text file containing coordinates for the clusters to lookup

- `coord_file_columns`:

  The columns in `coord_file` containing the x, y, and z coordinates.
  Note that these should be zero-based indices, so the first column is
  0, second column is 1, etc.

- `coord_vector`:

  Rather than lookup multiple locations from a file, `coord_vector`
  requires three numeric values (x, y, and z coordinates) that are used
  to lookup a single region in the atlases.

- `atlases`:

  a character vector of atlases in whereami to use for lookup

- `coord_orientation`:

  LPI or RAI

- `coord_space`:

  The template space for the coordinates. Default: 'MNI'

- `output_file`:

  The file name/path for the output of whereami

- `omask_output_file`:

  The file name/path for the output of whereami with the omask
  (percentage overlap)

- `afnidir`:

  The location of afni

------------------------------------------------------------------------

### Method `get_call()`

return the whereami call for this object

#### Usage

    afni_whereami$get_call()

------------------------------------------------------------------------

### Method `get_omask_call()`

return the whereami call for the omask lookup

#### Usage

    afni_whereami$get_omask_call()

------------------------------------------------------------------------

### Method `get_whereami_df()`

return the data.frame of coordinates and labels from the whereami lookup

#### Usage

    afni_whereami$get_whereami_df()

------------------------------------------------------------------------

### Method `get_overlap_df()`

return the data.frame of percent overlap between each atlas and
clusters/ROIs in the input mask

#### Usage

    afni_whereami$get_overlap_df()

------------------------------------------------------------------------

### Method `get_omask_df()`

return the data.frame of the percentage overlap for each parcel (not
implemented)

#### Usage

    afni_whereami$get_omask_df()

------------------------------------------------------------------------

### Method `get_outputs()`

return a list of the output files from this whereami object

#### Usage

    afni_whereami$get_outputs(exclude_missing = TRUE)

#### Arguments

- `exclude_missing`:

  if TRUE, any expected output file that does not exist will be NA in
  the returned list

#### Returns

a list of all output files

------------------------------------------------------------------------

### Method `run()`

Run whereami for this object

#### Usage

    afni_whereami$run(force = FALSE)

#### Arguments

- `force`:

  if TRUE, re-run a whereami call even if the expected output file
  exists

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    afni_whereami$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
