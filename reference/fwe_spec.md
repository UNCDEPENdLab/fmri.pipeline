# R6 class for an FWE correction method that can apply to one or more models

R6 class for an FWE correction method that can apply to one or more
models

## Methods

### Public methods

- [`fwe_spec$new()`](#method-fwe_spec-initialize)

- [`fwe_spec$submit()`](#method-fwe_spec-submit)

- [`fwe_spec$clone()`](#method-fwe_spec-clone)

------------------------------------------------------------------------

### `fwe_spec$new()`

create a new FWE correction instance

#### Usage

    fwe_spec$new(..., fwe_data = NULL)

#### Arguments

- `...`:

  not used yet

- `fwe_data`:

  not used yet

------------------------------------------------------------------------

### `fwe_spec$submit()`

submit FWE estimation to the cluster

#### Usage

    fwe_spec$submit()

------------------------------------------------------------------------

### `fwe_spec$clone()`

The objects of this class are cloneable with this method.

#### Usage

    fwe_spec$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
