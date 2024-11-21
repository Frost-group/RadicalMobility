#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=4:mem=4gb
#PBS -N main_python_job
#PBS -o main_output.log
#PBS -e main_error.log

TARGET_DIR="/rds/general/user/jh121/home" 

if [ "$PWD" != "$TARGET_DIR" ]; then
    echo "Switching to the target directory: $TARGET_DIR"
    cd "$TARGET_DIR" || { echo "Failed to change directory to $TARGET_DIR"; exit 1; }
fi

# Load the necessary module and activate the Python environment
module load anaconda3/personal
source /rds/general/user/jh121/home/miniforge3/etc/profile.d/conda.sh
conda activate xtb_env

# Run the Python script
python3 main.py