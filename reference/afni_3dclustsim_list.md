# R6 class for a list of 3dClustSim runs

R6 class for a list of 3dClustSim runs

R6 class for a list of 3dClustSim runs

## Methods

### Public methods

- [`afni_3dclustsim_list$new()`](#method-afni_3dclustsim_list-new)

- [`afni_3dclustsim_list$submit()`](#method-afni_3dclustsim_list-submit)

- [`afni_3dclustsim_list$get_objs()`](#method-afni_3dclustsim_list-get_objs)

- [`afni_3dclustsim_list$clone()`](#method-afni_3dclustsim_list-clone)

------------------------------------------------------------------------

### Method `new()`

create a new afni_3dclustsim_list object

#### Usage

    afni_3dclustsim_list$new(obj_list = NULL, ...)

#### Arguments

- `obj_list`:

  A list of afni_3dclustsim objects.

- `...`:

  One or more afni_3dclustsim objects. Used if `obj_list` is NULL.

------------------------------------------------------------------------

### Method `submit()`

submit all jobs in this list

#### Usage

    afni_3dclustsim_list$submit(force = FALSE)

#### Arguments

- `force`:

  If TRUE, pass force = TRUE to all the submit method of all subsidiary
  3dClustSim objects

------------------------------------------------------------------------

### Method `get_objs()`

return the clustsim objects

#### Usage

    afni_3dclustsim_list$get_objs()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    afni_3dclustsim_list$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
