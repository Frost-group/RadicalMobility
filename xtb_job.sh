#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=4:mem=4gb
#PBS -N xtb_job
#PBS -o xtb_output.log
#PBS -e xtb_error.log

TARGET_DIR="/rds/general/user/jh121/home" 

if [ "$PWD" != "$TARGET_DIR" ]; then
    echo "Switching to the target directory: $TARGET_DIR"
    cd "$TARGET_DIR" || { echo "Failed to change directory to $TARGET_DIR"; exit 1; }
fi

if [ ! -f ".UHFfrag" ]; then
    echo "Error: .UHFfrag file not found in $TARGET_DIR"
    exit 1
fi

module load anaconda3/personal
source /rds/general/user/jh121/home/miniforge3/etc/profile.d/conda.sh
conda activate xtb_env

export PATH="/rds/general/user/jh121/home/xtb_install/bin:$PATH"

/rds/general/user/jh121/home/xtb_install/bin/xtb optdimerradicaldimerBase3SF5.mol --dipro --chrg 0 --uhf 2