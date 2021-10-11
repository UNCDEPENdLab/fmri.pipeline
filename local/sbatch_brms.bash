#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=12g
#SBATCH -n 1
#SBATCH -t 10:00:00

module use /proj/mnhallqlab/sw/modules
module load r/4.0.3_depend

R CMD BATCH --no-save --no-restore rt_prediction_meta_analysis.R rt_prediction_meta_analysis.Rout
