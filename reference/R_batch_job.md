# Description of R_batch_job R6 class

Description of R_batch_job R6 class

Description of R_batch_job R6 class

## Public fields

- `parent_jobs`:

  a vector of parent_jobs that are upstream of this job and may
  influence its execution

- `depends_on_parents`:

  logical or character string indicating whether this job should wait
  until `parent_jobs` complete. If a character string is passed, it
  indicates which named elements of `parent_jobs` must complete before
  this job begins.

- `job_name`:

  a user-defined name for the job used for specifying job dependencies
  and informative job status queries on a job scheduler

- `wait_for_children`:

  if TRUE, code will be inserted to wait for all jobs in a vector called
  `child_job_ids` to finish before the batch exits. It's up to your code
  to use this variable name

- `all_children_successful`:

  if TRUE, all jobs in vector `child_job_ids` must be successful for
  this job to finish (see `wait_for_children` field)

- `wall_time`:

  The amount of time requested on the job scheduler, following
  d-hh:mm:ss format. Defaults to "4:00:00", which is 4 hours.

- `n_nodes`:

  The number of nodes to be requested on the job scheduler

- `n_cpus`:

  The number of cores (aka 'cpus', ignoring hyperthreading) to be
  requested on the job scheduler

- `mem_total`:

  The total amount of memory (RAM) requested by the job

- `mem_per_cpu`:

  The amount of memory (RAM) requested per cpu (total = mem_per_cpu \*
  n_cpus)

- `input_objects`:

  An environment containing all objects to be written to an RData object
  and passed to the batch job at execution

- `input_rdata_file`:

  The name of the environment to be loaded at the beginning of the R
  batch prior to executing

- `output_rdata_file`:

  The name of the environment to be saved at the end of the R batch
  execution, which can then be loaded by subsequent jobs.

- `sqlite_db`:

  File path to tracking SQLite database

- `batch_directory`:

  Location of batch scripts to be written

- `batch_code`:

  Shell code to be included in the batch script prior to the R code to
  be run. This can include module load statements, environment variable
  exports, etc.

- `r_code`:

  The R code to be executed by the scheduler. This can be a character
  vector that includes multiple R statements or an expression object
  containing the R code to be evaluated

- `post_children_r_code`:

  The R code to be executed after child jobs have completed. This can be
  a character vector that includes multiple R statements or an
  expression object containing the R code to be evaluated. Only relevant
  if wait_for_children = TRUE

- `r_packages`:

  The R packages to be loaded into the environment before job execution.
  These are loaded by pacman::p_load, which will install any missing
  packages before attempting to load

- `scheduler`:

  The job scheduler to be used for this batch. Options are: "slurm",
  "torque", or "local".

- `scheduler_options`:

  An optional character vector of scheduler arguments to be included in
  the batch script header that control additional features such as job
  emails or group permissions. These directives are added with \#SBATCH
  or \#PBS headings, depending on the scheduler, and are ignored if the
  scheduler is "local".

- `repolling_interval`:

  The number of seconds to wait between successive checks on whether
  parent jobs have completed. This is mostly relevant to the 'local'
  scheduler.

- `print_session_info`:

  If TRUE, print the \`sessionInfo()\` and \`Sys.info()\` when the job
  starts. Useful for debugging problems with the compute environment or
  R installation.

- `print_environment`:

  If \`TRUE\`, print the \`Sys.getenv()\` when the job starts. This can
  produce a lot of output, but can be useful if certain environment
  variables are not being found when your job runs, leading it to fail.

- `sqlite_cache_obj`:

  If \`TRUE\`, copy this object into the SQLite database at submission
  time. This is mostly useful for detailed debugging and is not
  generally recommended because it can increase the size of the database
  considerably

## Active bindings

- `sequence_id`:

  Used to identify sequence when using R_batch_sequence

## Methods

### Public methods

- [`R_batch_job$new()`](#method-batch_job-new)

- [`R_batch_job$generate()`](#method-batch_job-generate)

- [`R_batch_job$submit()`](#method-batch_job-submit)

- [`R_batch_job$copy()`](#method-batch_job-copy)

- [`R_batch_job$reset_file_names()`](#method-batch_job-reset_file_names)

- [`R_batch_job$get_job_id()`](#method-batch_job-get_job_id)

- [`R_batch_job$get_child_ids()`](#method-batch_job-get_child_ids)

- [`R_batch_job$clone()`](#method-batch_job-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new R_batch_job object for execution on an HPC cluster

#### Usage

    R_batch_job$new(
      batch_directory = NULL,
      parent_jobs = NULL,
      job_name = NULL,
      n_nodes = NULL,
      n_cpus = NULL,
      wall_time = NULL,
      mem_per_cpu = NULL,
      mem_total = NULL,
      r_code = NULL,
      r_script = NULL,
      post_children_r_code = NULL,
      batch_code = NULL,
      r_packages = NULL,
      scheduler = NULL,
      wait_for_children = NULL,
      all_children_successful = NULL,
      input_rdata_file = NULL,
      input_objects = NULL,
      output_rdata_file = NULL,
      sqlite_db = NULL,
      scheduler_options = NULL,
      repolling_interval = NULL,
      print_session_info = TRUE,
      print_environment = FALSE
    )

#### Arguments

- `batch_directory`:

  The location of batch scripts to be generated

- `parent_jobs`:

  A vector of one or more job ids that are parents to this job. This can
  be a named vector, to be used in conjunction with `depends_on_parents`
  to specify which parent jobs must be completed before this job begins.

- `job_name`:

  The name of the job used in dependency specification and job scheduler
  naming

- `n_nodes`:

  The number of compute nodes to be requested on the scheduler

- `n_cpus`:

  The number of cpus to be requested on the scheduler

- `wall_time`:

  The compute time requested on the cluster dd-HH:MM:SS

- `mem_per_cpu`:

  The amount of memory to be requested per cpu

- `mem_total`:

  The total amount of memory to requested by the job

- `r_code`:

  A character vector or expression containing R code to be executed

- `r_script`:

  The path to an R script to be executed by the batch (mutually
  exclusive with `r_code`).

- `post_children_r_code`:

  A character vector of R code to be executed after waiting for child
  jobs finishes

- `batch_code`:

  A character vector of code to be included in the batch script for job
  scheduling

- `r_packages`:

  A character vector of R packages to be loaded when compute script runs

- `scheduler`:

  The scheduler to be used for this compute. Options are 'slurm',
  'torque', or 'local'.

- `wait_for_children`:

  If TRUE, do not end this job until all child jobs have completed

- `all_children_successful`:

  If TRUE, don't count this job as successful unless all child jobs are
  successful

- `input_rdata_file`:

  The name of the environment to be loaded at the beginning of the R
  batch prior to executing code

- `input_objects`:

  A list object in the current execution environment to be cached and
  used as input to the R batch. This is mutually exclusive with
  input_rdata_file at present.

- `output_rdata_file`:

  The name of the environment to be saved at the end of the R batch
  execution

- `sqlite_db`:

  The location of the SQLite database to be used for job tracking. If
  \`NULL\`, job tracking will be disabled.

- `scheduler_options`:

  A character vector of scheduler options to be added to the header of
  the batch script

- `repolling_interval`:

  The number of seconds to wait before rechecking whether parent jobs
  have completed

- `print_session_info`:

  If TRUE (default), print information about the R environment
  \`sessionInfo()\` and compute environment \`Sys.info()\` when the job
  starts.

- `print_environment`:

  If TRUE, print the session environment via \`Sys.getenv()\` when the
  job starts. Default: FALSE.

------------------------------------------------------------------------

### Method `generate()`

Helper function that generates the batch and compute files for a job

#### Usage

    R_batch_job$generate(force = FALSE)

#### Arguments

- `force`:

  if TRUE, the RData, batch, and compute will be regenerated/rewritten

#### Details

this is called by `$submit` when a job is submitted and is provided here
in case the user wants to generate the batch files without executing
them

------------------------------------------------------------------------

### Method `submit()`

Submit job to scheduler or local compute

#### Usage

    R_batch_job$submit()

------------------------------------------------------------------------

### Method `copy()`

Function to create a deep copy of a batch job

#### Usage

    R_batch_job$copy(
      job_name = NULL,
      n_nodes = NULL,
      n_cpus = NULL,
      wall_time = NULL,
      mem_total = NULL,
      r_code = NULL,
      post_children_r_code = NULL
    )

#### Arguments

- `job_name`:

  The name of the job used in dependency specification and job scheduler
  naming

- `n_nodes`:

  The number of compute nodes to be requested on the scheduler

- `n_cpus`:

  The number of cpus to be requested on the scheduler

- `wall_time`:

  The compute time requested on the cluster dd-HH:MM:SS

- `mem_total`:

  The total amount of memory (RAM) requested by the job

- `r_code`:

  A character vector or expression containing R code to be executed

- `post_children_r_code`:

  The R code to be executed after child jobs have completed

#### Details

Note that this also resets the compute_file_name and batch_file_name
fields so that the copied object doesn't create files that collide with
the original. This method also exposes a few named parameters that can
be used to override the copied fields with new values to avoid needing
to change these one-by-one using obj\$\<field\> \<- x syntax

------------------------------------------------------------------------

### Method `reset_file_names()`

helper function to reset names of compute and batch files that will be
generated by this job.

#### Usage

    R_batch_job$reset_file_names()

#### Details

This needs to be exposed as a public method for copied objects to be
able to reset these private fields.

------------------------------------------------------------------------

### Method `get_job_id()`

Return the job id of this job (populated by job submission)

#### Usage

    R_batch_job$get_job_id()

------------------------------------------------------------------------

### Method `get_child_ids()`

Return the ids of all child jobs launched by this job

#### Usage

    R_batch_job$get_child_ids()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    R_batch_job$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
