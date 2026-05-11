# Add parent/child id relationship to tracking database

Add parent/child id relationship to tracking database

## Usage

``` r
add_tracked_job_parent(
  sqlite_db = NULL,
  job_id = NULL,
  parent_job_id = NULL,
  child_level = 1
)
```

## Arguments

- sqlite_db:

  Path to SQLite database used for tracking

- job_id:

  Job id of job for which to add a parent

- parent_job_id:

  Job id of the parent job to job_id

- child_level:

  Level of child; currently supports two levels
