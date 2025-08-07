
test_data_base_dir = "/proj/mnhallqlab/projects/fmri.pipeline_test_data"

# Read in trial, run, and subject dataframes
trial_df <- read.csv(file.path(test_data_base_dir, "sample_trial_data.csv.gz"))
run_df <- read.csv(file.path(test_data_base_dir, "sample_run_data.csv"))
subj_df <- read.csv(file.path(test_data_base_dir, "sample_subject_data.csv"))

gpa <- setup_glm_pipeline(
  analysis_name = "gpa_tests",
  scheduler = "slurm",
  trial_data = trial_df,
  run_data = run_df,
  subject_data = subj_df, # output_directory = tempdir(), # vm = c(id = "id", run_number = "run"),
  n_expected_runs = 8,
  tr = 1.0,
  drop_volumes = 2,
  l1_models=NULL, l2_models=NULL, l3_models=NULL,
  confound_settings = list(
    motion_params_file = "motion.par",
    confound_input_file = "nuisance_regressors.txt",
    confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
    l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
    exclude_run = "max(framewise_displacement) > 5 | sum(framewise_displacement > .9)/length(framewise_displacement) > .10",
    exclude_subject =  "n_good_runs < 4",
    truncate_run = "(framewise_displacement > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL
  )
)

gpa <- setup_compute_environment(gpa, preselect_action = 5L)  

# Build L1 models
gpa <- build_l1_models(gpa, from_spec_file = "/proj/mnhallqlab/projects/fmri.pipeline_test_data/sample_2_L1_spec.yaml")

# Build L2 models
gpa <- build_l2_models(gpa) #, from_spec_file = file.path(test_data_base_dir, spec_file))

# Build L3 models
gpa <- build_l3_models(gpa) #, from_spec_file = file.path(test_data_base_dir, spec_file))
