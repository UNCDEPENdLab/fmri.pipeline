#!/usr/bin/env sh

#PBS -A mnh5174_c_g_sc_default
#PBS -l nodes=1:ppn=20
#PBS -l pmem=8gb
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m n
#PBS -W group_list=mnh5174_collab

cd $PBS_O_WORKDIR

export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

#module load r/3.5.0
#module load fsl/5.0.11
#module load afni/18.1.15
#module load openblas/0.2.20

module load r/3.6.0
module load fsl/6.0.1
module load afni/19.0.26
module load gsl/2.5

ni_tools="$G/lab_resources"

#location of MRI template directory for preprocessing
MRI_STDDIR="${ni_tools}/standard"

#add preprocessing scripts that may be called in this pipeline
PATH="${ni_tools}/c3d-1.1.0-Linux-x86_64/bin:${ni_tools}/fmri_processing_scripts:${ni_tools}/fmri_processing_scripts/autopreproc:${ni_tools}/bin:${PATH}"

export PATH MRI_STDDIR

#the fsl_pipeline_file environment variable must be passed in through qsub, which is picked up by the R script
#run_model_index is also passed in to determine which model variant to run within sceptic_run_variant

if [ -z "${R_SCRIPT}" ]; then
    echo "R_SCRIPT variable not set"
    exit 1
fi

if [ ! -f "${R_SCRIPT}" ]; then
    echo "R_SCRIPT file does not exist"
    exit 1
fi

if [ -d "$( dirname $R_SCRIPT )/outputs" ]; then
    outfile="outputs/${R_SCRIPT/.R/.Rout}"
else
    outfile="${R_SCRIPT/.R/.Rout}"
fi

if [ -n "${run_model_index}" ]; then
    outfile="${outfile/.Rout/_${run_model_index}.Rout}" #separate .Rout files by pipeline/index
fi

R CMD BATCH --no-save --no-restore "${R_SCRIPT}" "${outfile}"
