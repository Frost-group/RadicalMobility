#!/bin/bash
#PBS -N dftb_xtb_annealing
#PBS -l select=1:ncpus=32:mem=256gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# Load required modules (Modify based on your HPC setup)
module load miniforge/3  
eval "$(~/miniforge3/bin/conda shell.bash hook)"   

# Activate your Conda environment
source /rds/general/user/jh121/home/miniforge3/bin/activate xtb_env

# Navigate to the working directory
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export OMP_STACKSIZE=2G


python annealing_xtb.py
