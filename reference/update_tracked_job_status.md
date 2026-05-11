# Update job status in tracking SQLite database

Update job status in tracking SQLite database

## Usage

``` r
update_tracked_job_status(
  sqlite_db = NULL,
  job_id = NULL,
  status,
  cascade = FALSE,
  exclude = NULL
)
```

## Arguments

- sqlite_db:

  Character string specifying the SQLite database used for job tracking

- job_id:

  Character string specifying the job id to update as failed

- status:

  Character string specifying the job status to set. Must be one of:
  "QUEUED", "STARTED", "FAILED", "COMPLETED", "FAILED_BY_EXT"

- exclude:

  Any job ids to ignore when cascading a status
