# Description of R_batch_sequence R6 class

Description of R_batch_sequence R6 class

## Public fields

- `sequence_id`:

  Optional identifier for all jobs in this sequence, used for job
  tracking

## Methods

### Public methods

- [`batch_sequence$new()`](#method-batch_sequence-initialize)

- [`batch_sequence$add()`](#method-batch_sequence-add)

- [`batch_sequence$submit()`](#method-batch_sequence-submit)

- [`batch_sequence$generate()`](#method-batch_sequence-generate)

- [`batch_sequence$get_job_ids()`](#method-batch_sequence-get_job_ids)

- [`batch_sequence$clone()`](#method-batch_sequence-clone)

------------------------------------------------------------------------

### `batch_sequence$new()`

create a new R_batch_sequence object

#### Usage

    batch_sequence$new(..., joblist = NULL, sequence_id = NULL)

#### Arguments

- `...`:

  One or more R_batch_job objects to be run in sequence

- `joblist`:

  Optional list of jobs to be used instead of ...

- `sequence_id`:

  Optional identifier for all jobs in this sequence; used for job
  tracking

------------------------------------------------------------------------

### `batch_sequence$add()`

add one or more R_batch_job objects to the sequence

#### Usage

    batch_sequence$add(...)

#### Arguments

- `...`:

  One or more R_batch_job objects to be added to sequence

------------------------------------------------------------------------

### `batch_sequence$submit()`

submit the job sequence to the scheduler or local compute

#### Usage

    batch_sequence$submit()

------------------------------------------------------------------------

### `batch_sequence$generate()`

Calls each job's \$generate() method so that scripts can be examined
without running the sequence

#### Usage

    batch_sequence$generate()

------------------------------------------------------------------------

### `batch_sequence$get_job_ids()`

Return a named vector of the job ids for all jobs in this sequence

#### Usage

    batch_sequence$get_job_ids()

------------------------------------------------------------------------

### `batch_sequence$clone()`

The objects of this class are cloneable with this method.

#### Usage

    batch_sequence$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
