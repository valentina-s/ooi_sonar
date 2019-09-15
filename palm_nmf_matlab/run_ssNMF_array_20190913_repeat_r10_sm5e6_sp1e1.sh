#!/bin/bash
# This script runs multiple single core Matlab jobs using
# a slurm job array 

#SBATCH --array=1-50
#SBATCH --time=20:00:00
#SBATCH --mem=6GB
#SBATCH --ntasks-per-node=1
#SBATCH --output=log/slurm_%A-%a.out
#SBATCH --error=log/slurm_%A-%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=leewujung@gmail.com

# Load Matlab environment
module load matlab/2018b

# Create a temporary directory on scratch for any Job related files
slurmArrayID="${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}"
export slurmArrayID
mkdir -p /scratch/ch153/wjl/slurmJobs/$slurmArrayID

PARAM_FILE="'/scratch/ch153/wjl/params_20190913_r10_sm5e6_sp1e1.txt'"
matlab -nodisplay -nosplash -singleCompThread -r "ssNMF_runner_repeat($PARAM_FILE, $SLURM_ARRAY_TASK_ID);exit;"

# remove workspace
rm -rf /scratch/ch153/wjl/slurmJobs/$slurmArrayID
