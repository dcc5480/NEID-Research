#!/bin/bash
#SBATCH --partition=open
##time in h:mm:ss
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --output=/storage/home/dcc5480/work/logs/multi_jobarray_%A_%a.log

#SBATCH --mail-user=dcc5480@psu.edu
#SBATCH --mail-type=ALL

echo "Starting job $SLURM_JOB_NAME"
echo "Job id: $SLURM_JOB_ID"
echo "Submit Directory: $SLURM_SUBMIT_DIR"
date

echo "running Python sciprt..."
source /storage/home/dcc5480/work/venv1/bin/activate

/storage/home/dcc5480/work/bin/python3 multi.py
date
echo "Finished running Python script"
