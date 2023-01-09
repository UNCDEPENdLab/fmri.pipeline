#!/bin/bash

#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_136_2.3mm.nii.gz"  sbatch_medusa.bash
#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/vmpfc_ah/vmpfc_ba_max_2.3mm_filled.nii.gz"  sbatch_medusa.bash
#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/pauli_da_midbrain_2.3mm.nii.gz" sbatch_medusa.bash
#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/bilateral_striatum_tight_7Networks_2.3mm.nii.gz" sbatch_medusa.bash
#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/pauli_combined_integermask_2.3mm.nii.gz" sbatch_medusa.bash


#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_Cont_2.3mm.nii.gz"  sbatch_medusa.bash
#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_Default_2.3mm.nii.gz"  sbatch_medusa.bash
sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_DorsAttn_2.3mm.nii.gz"  sbatch_medusa.bash
#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_Limbic_2.3mm.nii.gz"  sbatch_medusa.bash
sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_SalVentAttn_2.3mm.nii.gz"  sbatch_medusa.bash
#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_SomMot_2.3mm.nii.gz"  sbatch_medusa.bash
#sbatch --export=atlas="/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_Vis_2.3mm.nii.gz"  sbatch_medusa.bash

