#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=72g
#SBATCH -n 18
#SBATCH -t 2-

module use /proj/mnhallqlab/sw/modules
module load r/4.0.1
module load fsl/6.0.4
module load afni/20.1.05

export atlas #pass through env variable
suffix=$( basename $atlas )
suffix=${suffix/.nii.gz/}

R CMD BATCH --no-save --no-restore extract_sceptic_atlas_deconvolved.R "extract_sceptic_atlas_deconvolved_${suffix}.Rout"
