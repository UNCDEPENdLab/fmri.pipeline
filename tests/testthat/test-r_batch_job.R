# tests for R batch job, including SQLite tracker

# library(DBI)
# library(RSQLite)
# db_file <- "job_tracker.sqlite"
# mydb <- dbConnect(RSQLite::SQLite(), db_file)


# job_spec_sql <- "
# CREATE TABLE batch_job (
#     id INTEGER PRIMARY KEY,
#     parent_id INTEGER,
#     job_id VARCHAR NOT NULL UNIQUE,
#     job_name VARCHAR,
#     batch_directory VARCHAR,
#     n_nodes INTEGER CHECK (n_nodes >= 1),
#     n_cpus INTEGER CHECK (n_cpus >= 1),
#     wall_time VARCHAR,
#     mem_per_cpu VARCHAR,
#     mem_total VARCHAR,
#     r_code TEXT,
#     scheduler VARCHAR,
#     scheduler_options VARCHAR,
#     job_obj BLOB,
#     time_submitted INTEGER,
#     time_started INTEGER,
#     time_ended INTEGER,
#     status VARCHAR(24),
#     FOREIGN KEY (parent_id) REFERENCES batch_job (id)
# );
# "


# job_spec_sql <- "
# CREATE TABLE batch_job (
#   id INTEGER PRIMARY KEY,
#   parent_id INTEGER,
#   job_id varchar(128) NOT NULL UNIQUE,
#   job_name TEXT,
#
#   job_obj BLOB,
#   time_submitted INTEGER,
#   time_started INTEGER,
#   time_ended INTEGER,
#   status VARCHAR(24),
#   FOREIGN KEY (parent_id) References batch_job (id)
# );
# "

# dbSendStatement(mydb, job_spec_sql)

# dbGetQuery(mydb, "SELECT * FROM batch_job;")


# #this file has the SQL syntax to setup (and reset) the database
# reset_sql <- "
# SET foreign_key_checks = 0;
# DROP TABLE IF EXISTS batch_job;
# SET foreign_key_checks = 1;
# "



library(glue)
test_df <- data.frame(x=1:10)
library(fmri.pipeline)
d_batch <- R_batch_job$new(
  job_name = glue("test1"), n_cpus = 1, mem_per_cpu = "4g",
  wall_time = "10:00:00", scheduler = "sbatch",
  # pass relevant vars to the batch
  input_objects = fmri.pipeline:::named_list(test_df),
  r_packages = "fmri.pipeline",
  r_code = c(
    "Sys.sleep(1)"
  )
)

d_batch$submit()


# local batch
d_batch <- R_batch_job$new(
  job_name = glue("test1"), n_cpus = 1, mem_per_cpu = "4g",
  wall_time = "10:00:00", scheduler = "local",
  sqlite_db = "~/job_tracker.sqlite",
  # pass relevant vars to the batch
  input_objects = fmri.pipeline:::named_list(test_df),
  r_packages = "fmri.pipeline",
  r_code = c(
    "Sys.getenv('JOBID')",
    "Sys.sleep(8)"
  )
)

d_batch$submit()
sqlite_db <- "~/job_tracker.sqlite"

library(DBI)
con <- dbConnect(RSQLite::SQLite(), sqlite_db)
job_df <- dbReadTable(con, "batch_job")
if (nrow(job_df) > 0L) job_df$job_obj <- lapply(job_df$job_obj, function(x) if (!is.null(x)) unserialize(x))
dbDisconnect(con)

get_batch_job_status("78929", return_children=FALSE, sqlite_db)




# debug update query
con <- dbConnect(RSQLite::SQLite(), sqlite_db)
sql <- glue::glue("UPDATE batch_job SET STATUS = 'COMPLETED', time_ended = '{Sys.time()}' WHERE job_id == '77932'")
#sql <- glue::glue("UPDATE batch_job SET STATUS = 'COMPLETED' WHERE job_id == '77932'")
dbExecute(con, sql)
dbDisconnect(con)

update_batch_job_status(sqlite_db, "77932", "COMPLETED")

# # SQL schema for job ids
res <- dbSendQuery(con, "SELECT last_insert_rowid()")
dbFetch(res)

res <- dbSendQuery(con, "INSERT INTO batch_job last_insert_rowid()")
dbFetch(res)

res <- dbSendStatement(con, "UPDATE Cellar SET WS_ID = ? WHERE Cellar_ID = ?", param=list(ws_id_added, ws_import_cellarid))
stopifnot(dbGetRowsAffected(res) == 1L) #make sure we have updated a row
dbClearResult(res)

# basic ideas

# 1) always add code to beginning of job to update time_started. Basically a command to execute an update to the SQLITE
# 2) always add code to the end of the job to update time_ended and job_status. If the script runs to the last line, then we know that the script completed successfully
# 3) trap exit with on.exit() for failures
# 4) potentially trap any other problems in the batch script itself, though that would require some sort of command line SQLITE update

# test batch sequence


# library(glue)
# library(fmri.pipeline)
d_batch <- R_batch_job$new(
  job_name = glue("test1"), n_cpus = 1, mem_per_cpu = "4g",
  wall_time = "10:00:00", scheduler = "local",
  # pass relevant vars to the batch
  r_packages = "fmri.pipeline",
  r_code = c(
    "Sys.sleep(10)"
  )
)

d_batch2 <- d_batch$copy(job_name = "test2")
d_batch3 <- d_batch$copy(job_name = "test3")

d_seq <- R_batch_sequence$new(
  d_batch, d_batch2, d_batch3
)

d_seq$submit()
jids <- d_seq$get_job_ids()
