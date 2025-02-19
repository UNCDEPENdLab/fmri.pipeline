# ---- testing with l1, l2 (emotion), l3 (age) models ----
gpa <- create_gpa(
    analysis_name = "gpa_tests_l2l3_18feb2025",
    test_data_base_dir = "/proj/mnhallqlab/projects/fmri.pipeline_test_data",
    l1_spec_file = "sample_2_L1_spec.yaml",
    trial_data_file = "sample_trial_data.csv.gz",
    run_data_file = "sample_run_data.csv",
    subject_data_file = "sample_subject_data.csv",
    gpa_cache_file = "gpa_l2l3_18feb2025.rds",
    cache_file = "gpa_base_l2l3_18feb2025.rds") 

gpa$parallel$finalize_time <- "10:00:00"
gpa$parallel$fsl$l2_feat_time <- "4:00:00"
gpa$parallel$fsl$l3_feat_time <- "50:00:00" #extend l3_feat_time

run_glm_pipeline(gpa)

# gzipping the nifti input files

# find files in all subfolders in /proj/mnhallqlab/projects/fmri.pipeline_test_data/ which have nfaswuktm_clock...ni
# nifti_files <- list.files(
#     path = "/proj/mnhallqlab/projects/fmri.pipeline_test_data/",
#     pattern = "nfaswuktm_clock.*\\.nii",
#     recursive = TRUE,
#     full.names = TRUE
# )
# # apply this system("gzip -f /proj/mnhallqlab/projects/fmri.pipeline_test_data/derivatives/gpa_tests_l2l3/first_level/subject_1/run_1/cope1.nii") to all nifti_files
# sapply(nifti_files, function(x) {
#     dest_file <- paste0(x, ".gz")
#     system(glue("gzip -c {x} > {dest_file}"))
# })

# change nifti filepath to the current file .nii.qz
# test_data_base_dir <- "/proj/mnhallqlab/projects/fmri.pipeline_test_data/"
# run_data_file <- "sample_run_data.csv"
# run_df <- read.csv(file.path(test_data_base_dir, run_data_file))
# run_df <- run_df %>% mutate(run_nifti = gsub(".nii.gz", ".nii", run_nifti))
# write.csv(run_df, file.path(test_data_base_dir, run_data_file), row.names = FALSE)
