# Description of R_batch_sequence R6 class

Description of R_batch_sequence R6 class

Description of R_batch_sequence R6 class

## Public fields

- `sequence_id`:

  Optional identifier for all jobs in this sequence, used for job
  tracking

## Methods

### Public methods

- [`R_batch_sequence$new()`](#method-batch_sequence-new)

- [`R_batch_sequence$add()`](#method-batch_sequence-add)

- [`R_batch_sequence$submit()`](#method-batch_sequence-submit)

- [`R_batch_sequence$generate()`](#method-batch_sequence-generate)

- [`R_batch_sequence$get_job_ids()`](#method-batch_sequence-get_job_ids)

- [`R_batch_sequence$clone()`](#method-batch_sequence-clone)

------------------------------------------------------------------------

### Method `new()`

create a new R_batch_sequence object

#### Usage

    R_batch_sequence$new(..., joblist = NULL, sequence_id = NULL)

#### Arguments

- `...`:

  One or more R_batch_job objects to be run in sequence

- `joblist`:

  Optional list of jobs to be used instead of ...

- `sequence_id`:

  Optional identifier for all jobs in this sequence; used for job
  tracking

------------------------------------------------------------------------

### Method `add()`

add one or more R_batch_job objects to the sequence

#### Usage

    R_batch_sequence$add(...)

#### Arguments

- `...`:

  One or more R_batch_job objects to be added to sequence

------------------------------------------------------------------------

### Method `submit()`

submit the job sequence to the scheduler or local compute

#### Usage

    R_batch_sequence$submit()

------------------------------------------------------------------------

### Method `generate()`

Calls each job's \$generate() method so that scripts can be examined
without running the sequence

#### Usage

    R_batch_sequence$generate()

------------------------------------------------------------------------

### Method `get_job_ids()`

Return a named vector of the job ids for all jobs in this sequence

#### Usage

    R_batch_sequence$get_job_ids()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    R_batch_sequence$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
