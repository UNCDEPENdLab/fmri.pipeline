library(fmri.pipeline)
# note that HPC settings like nodes and batch_code are ignored for local scheduler
w <- R_batch_job$new(
  batch_directory="~/batch_test",
  job_name="step1",
  n_nodes=1,
  n_cpus=1,
  cpu_time="10:00",
  r_code=c(
    "Sys.sleep(1)",
    "print('hi')",
    "x <- 2+2"
  ),
  r_packages=c("lme4"),
  batch_code = c(
    "module use /proj/mnhallqlab/sw/modules",
    "module load r/4.0.3_depend"
  ),
  scheduler="local",
  repolling_interval=1 #repoll every second for this test (should be slower for real jobs)
)

#individual job file creation and execution
#x$generate() #to generate files but not submit
#x$submit() #both generates files and submits compute

# build jobs from each other (save-as style)
x <- w$copy()
x$job_name <- "step2"
x$depends_on_parents <- "step1"

y <- w$copy()
y$job_name <- "step3"
y$depends_on_parents <- "step2" # depend on step2, which in turn depends on step1

waiter <- w$copy()
waiter$job_name <- "local_wait"
waiter$r_code <- c(
  "Sys.sleep(20)"
)
waiter$generate() #generate script for deferred local submission

z <- w$copy()
z$job_name <- "step4"
z$depends_on_parents <- c("step3") # step4 depends on step1 completion only
z$wait_for_children <- TRUE # wait for all child jobs to complete
z$r_code <- c(
  z$r_code,
  "j1 <- fmri.pipeline::cluster_job_submit('/nas/longleaf/home/mnhallq/batch_test/submit_batch_local_wait.sh', scheduler='local')",
  "j2 <- fmri.pipeline::cluster_job_submit('/nas/longleaf/home/mnhallq/batch_test/submit_batch_local_wait.sh', scheduler='local')",
  "child_job_ids <- c(j1, j2)",
  "cat('child job ids', j1, j2)"
)

# this should wait until the children of step4 complete -- this is accomplished by step4 
# not ending until its children are complete
a <- w$copy()
a$job_name <- "step5"
a$depends_on_parents <- "step4"

batch_seq <- R_batch_sequence$new(w, x, y, z, a)

# submit batch sequence
library(tictoc)
tic()
batch_seq$submit()
toc()
