# fmri.pipeline
R package for running full fMRI GLM analyses in FSL and other software.

## overview

This folder contains all scripts for running a model-based fMRI analysis using FSL FILM GLS for level 1 analyses, OLS
for level 2 analyses (combining runs), and and FLAME 1+2 for level 3 analyses. The pipeline is structured to run on a
high-performance computing cluster using the TORQUE scheduler. The pipeline will spawn one job submission for each first-level
(run) model, then setup job dependencies for the level 2 and level 3 analyses so that the pipeline completes in a funnel-like
fashion (many first-level runs funneling into a few group models).

For the pipeline to work, processed MRI data should be structured as follows:

`<fmri_dir>/<subject_id>/<expectdir>/<paradigm_name><run_number>/<expectfile>`

For descriptions of these field names, see documentation below for `run_fsl_pipeline.R`.

In general, most scripts in the pipeline expect at least two environment variables to be set (passed using qsub -v in the cluster environment). The first, `fsl_pipeline_file` is the full path to the analysis specifications defined in `run_fsl_pipeline.R` and saved to an RData object in the configuration_files subfolder. The second, `run_model_index` specific element within `l1_model_variants` that is the current object of analysis. That is parallelization is achieved in part by spawning many jobs that are split in terms of which level 1 model is being executed. These environment variables should be passed by the `push_pipeline` function called in run_fsl_pipeline.R, but knowing about their existence is helpful for setting up or debugging subsidiary scripts (e.g., the beta extraction scripts).

## core scripts

### `run_fsl_pipeline.R`
  This is the primary script for running the entire voxelwise GLM pipeline. It constructs a full model specification in a list called fsl_model_arguments, which is saved as an RData object and passed to subsidiary scripts. The core ingredients of that object are described below.

  Once the fsl_model_arguments object is created, it is saved, then `push_pipeline` executes the full analysis pipeline based on the specification.

### `extract_sceptic_map_betas.R`
  This script looks for whole-brain-significant clusters in each covariate/predictor included in a level 3 analysis. It uses `zthresh` and `clustsize` settings in the `fsl_model_arguments` list to find the clusters, employing AFNI's 3dclust program specifically.

  Importantly, for a three-level analysis in FSL, each level 1 contrast becomes the dependent variable in a level 3 model such that there are we need a third-level FLAME analysis for each level 1 contrast, yielding one .gfeat folder per level 1 contrast. This script anticipates this structure and walks over all level 1 and level 3 regressors.

  This is accomplished by a for loop over the level 1 contrasts (`for (l1 in 1:n_l1_copes) {`), level 2 contrasts (`for (l2 in 1:n_l2_contrasts) {`), and level 3 contrats (`for (l3 in 1:n_l3_copes) {`). For each combination, 3dclust is run to identiy whole-brain significant clusters. The coordinates and labels are identified using AFNI whereami, then the cluster mask (integer-valued, one per cluster) is read into R. The mean regression coefficient for each subject and/or run (if `calculate_l1_betas` is `TRUE`) is extracted for each cluster. These are compiled into a long-form tidy data.frame (ending in 'betas.csv'). The script generates on cluster file per level 1 contrast, with columns in that file denoting level 2 and level 3 contrasts.

  For each level 1 contrast, two tidy CSV are produced, one labeled with 'metadata' and the other with 'betas'.

### `extract_sceptic_map_betas_atlas.R`
  This script extracts voxelwise regression coefficients for each voxel within an a priori set of mask files, `atlas_files`. Whereas `extract_sceptic_map_betas.R` extracts one coefficient within each significant cluster, `extract_sceptic_map_betas_atlas.R` provides regression coefficients for each voxel in a mask. It produces csv files labeled with `*atlas_betas.csv.gz` for the regression coefficients or `*atlas_zstats.csv.gz` for the corresponding voxelwise z-statistics (from FSL level 1 analysis). The structure of these files is very close to the voxelwise deconvolution outputs from `extract_sceptic_atlas_deconvolved.R`, but those contain BOLD/neural estimates at each time point, while this script is just one coefficient for a single subject (from the second-level GLM). The files produced by this script formed the basis of voxelwise analyses along the hippocampal long axis depicted in Figure 2 of Dombrovski, Luna, and Hallquist.

  - numid: numeric ID of subject (ascending); match with design.fsf
  - vnum: voxel number in the atlas file
  - beta: regression coefficient at this voxel
  - atlas_value: numeric value stored at this voxel within the atlas (e.g., position along the hippocampal long axis on the 0--1 scale)
  - x: L/R coordinate in MNI space
  - y: A/P coordinate in MNI space
  - z: I/S coordinate in MNI space
  - atlas_name: name of NIfTI atlas for extraction
  - ID: subject id
  - fsldir: second-level (subject) GLM directory used for extraction
  - l1_contrast: level 1 regression contrast (for beta)
  - l2_contrast: level 2 regression contrast (for beta)

### `extract_sceptic_atlas_deconvolved.R`

  This script walks over a set of atlas NIfTI files, `atlas_files`, runs a voxelwise deconvolution algorithm on each BOLD acquisition in the study (i.e., multiple runs per subject, multiple subjects). The algorithm used in the Dombrovski, Luna, & Hallquist paper is from Bush and Cisler (2013) and is called Bu12 in that paper. We have largely adapted an optimized C++ codebase provided by Keith Bush for speeding up the voxelwise calculation. This program, called `deconvolvefilter`, is committed in compiled form to Github, but is an ELF binary that will only run on Linux-based machines. The fallback position would be to run the `deconvolve_nlreg` function in the `dependlab` R package, which is the same algorithm, but a mixture of C++ and R (and slower).

  The program produces two long-form data.frame for each run of data, one containing the deconvolved time series (labeled 'deconvolved') and the other containing the input timeseries (i.e., the preprocessed fMRI data without deconvolution; labeled 'original'). The structure of these data.frames is:

  - subid: subject ID
  - run_num: run number
  - contingency: reinforcement contingency
  - emotion: central stimulus
  - time: time since acquisition start, in seconds
  - atlas_name: name of NIfTI atlas for extraction
  - atlas_value: numeric value stored at this voxel within the atlas (e.g., position along the hippocampal long axis on the 0--1 scale)
  - vnum: voxel number in the atlast file
  - x: L/R coordinate in MNI space
  - y: A/P coordinate in MNI space
  - z: I/S coordinate in MNI space
  - decon: deconvolved time series value
  - BOLD_z: (only in 'original files') normalized BOLD time series value


## fsl_model_arguments object

The fsl_model_arguments file gets saved by run_fsl_pipeline.R and is a critical input to most subsidiary scripts in the pipeline and in additional analysis scripts. Its structure is as follows:

- analysis name: This is a string used to name the resulting group analysis folder. Set it to something that describes the overall dataset, preprocessing approach, and model-based input.

- trial_statistics: this is a long-form data.frame containing trial statistics for each subject to be analyzed. Trials and subjects should be stacked on the rows, leading to a data.frame with subjects x trials total rows. This data.frame is used to calculate model-based fMRI regressors and setup the GLM design matrices for each subject.

- subject_covariates: A data.frame containing any covariates that could be used in the group analysis. There should be one row per subject.

- id_col: a string specifying the column name of the subject ID in all data provided as input (esp. subject_covariates and trial_statistics). This is used to merge data together as needed.

- fmri_dir: The root path to preprocessed data for the entire project

- expectdir: The subfolder name within `fmri_dir` containing preprocessed fMRI data for a given subject.

- expectfile: A regular expression denoting the expected filename for the preprocessed fMRI data to be entered into first-level GLMs.

- usepreconvolve: a boolean indicating whether regressors should be convolved with the HRF prior to entry into FSL FEAT. Convolution is achieved using the `build_design_matrix` function in the [dependlab R package](https://github.com/PennStateDEPENdLab/dependlab). I recommend TRUE, which is better tested. If you use FALSE, the pipeline will output so-called three-column timing files (onset, duration, parametric value), which are entered into FEAT, and convolution occurs within FEAT.

- ncpus: controls how many CPUs are used when running FEAT level 2 analysis (OLS combination of runs). I should probably rename this for clarity.

- drop_volumes: controls the number of initial volumes that are dropped from each run of fMRI data before GLM analysis. Useful for eliminating volumes before steady state magnetization.

- tr: the repetition time (in seconds) for the fMRI scan sequence. This is important for various temporal filtering steps to be applied correctly.

- spikeregressors: a boolean controlling whether 1/0 regressors for each large head movement are added as confound regressors in the GLM. Handled downstream in model_clock_fmri.

- l1_model_variants: a list where each element is a vector of signals that should be included simultaneously in the level 1 analysis. For each element of the list, the level 1 data are analyzed with the requested combination of regressors, and group models for this combination are also conducted. Thus, the total number of models is the product of l1_model_variants (level 1) and group_model_variants (level 3).

- group_model_variants: a list where each element is a vector of covariates that should be included in a group analysis. The names refer to columns in the `subject_covariates` data.frame.

- l1_contrasts: a named list whose elements are lists specifying a set of level 1 contrasts that should be included. The name of the vector is used to cross-reference the corresponding name in `l1_model_variants`. For example, if a model is named 'm1' in `l1_model_variants`, you can add level 1 contrasts to this model using 'm1' as the name of the list in `l1_contrasts`.

- execute_feat: a boolean indicating whether to execute FEAT for the level 1 analysis after creating the FSF file. This is somewhat vestigial from a pipeline variant that was run on a desktop computer. For cluster-based execution, this should be FALSE, which creates the FSF file, but hands off execution to the scheduler.

- model_suffix: an optional string that is appended to the folder names for each output the level 1 analyses specified in `l1_model_variants`. This is useful if you are running two or more models with the same level 1 regressors, but the model is nevertheless somewhat different (e.g., different contrasts).

- root_workdir: location for temporary job scripts generated for each run and model. These are not cleaned up after the pipeline completes (since it is distribted), but they can be deleted manually if you don't need them.

- n_cluster_beta_cpus: number of CPUs used simultaneously to extract regression coefficients (betas) from single runs/subjects based on significant group-level clusters.

- badids: an optional vector of IDs that should be excluded from analysis. These are dorpped before the pipeline begins.

- zthresh: the minimum z-stastistic threshold to be used in identifying whole-brain-significant clusters in extraction scripts. If not passed, defaults to |z| > 3.09.

- clustsize: the minimum number of voxels above `zthresh` that must be present for the cluster to be retained for extraction. This is the second component of a cluster-based thresholding approach for whole-brain significance (e.g., using 3dClustSim). If not specified, defaults to 34.
