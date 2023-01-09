library(fmri.pipeline)
# note that HPC settings like nodes and batch_code are ignored for local scheduler
w <- R_batch_job$new(
  batch_directory="~/batch_test",
  job_name="step1",
  n_nodes=1,
  n_cpus=1,
  cpu_time="10:00",
  r_code=c(
    "Sys.sleep(100)",
    "print('hi')",
    "x <- 2+2"
  ),
  r_packages=c("lme4"),
  batch_code = c(
    "module use /proj/mnhallqlab/sw/modules",
    "module load r/4.0.3_depend"
  ),
  scheduler="slurm",
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
y$depends_on_parents <- "step2" # step2 and step3 both depend on step 1, but do not depend on each other

waiter <- w$copy()
waiter$job_name <- "local_wait"
waiter$r_code <- c(
  "Sys.sleep(15)"
)
waiter$generate()

z <- w$copy()
z$job_name <- "step4"
z$depends_on_parents <- c("step1") # step4 depends on step1 completion only
z$wait_for_children <- TRUE
z$r_code <- c(
  z$r_code,
  "cat('Sys.sleep(15)', file='sleep15.R')",
  "j1 <- fmri.pipeline::cluster_job_submit('/nas/longleaf/home/mnhallq/batch_test/submit_batch_local_wait.sh', scheduler='slurm')",
  "j2 <- fmri.pipeline::cluster_job_submit('/nas/longleaf/home/mnhallq/batch_test/submit_batch_local_wait.sh', scheduler='slurm')",
  "child_job_ids <- c(j1, j2)"
)

batch_seq <- R_batch_sequence$new(w, z) # x, y,

# submit batch sequence
system.time(batch_seq$submit())
