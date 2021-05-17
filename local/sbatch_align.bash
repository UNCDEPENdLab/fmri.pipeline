#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=96g
#SBATCH -n 18
#SBATCH -t 1-

module use /proj/mnhallqlab/sw/modules
module load r/4.0.1
module load fsl/6.0.4
module load afni/20.1.05

#export atlas #pass through env variable
#suffix=$( basename $atlas )
#suffix=${suffix/.nii.gz/}

export ncores=18
R CMD BATCH --no-save --no-restore decon_alignment.R decon_alignment.Rout
