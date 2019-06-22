#!/bin/bash
# The slurm script file runParallelMultiple.slurm runs 
# multiple single core Matlab jobs using s Slurm job array 

#SBATCH --array=0-9
#SBATCH --time=16:00:00
#SBATCH --mem=6GB
#SBATCH --ntasks-per-node=1
#SBATCH --output=log/slurm_%A-%a.out
#SBATCH --error=log/slurm_%A-%a.err

# Load Matlab environment
module load matlab/2018b

# Create a temporary directory on scratch for any Job related files
slurmArrayID="${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}"
export slurmArrayID
mkdir -p /scratch/ch153/wjl/slurmJobs/$slurmArrayID

SM_ARRAY=(1e5 2e5 5e5 1e6 2e6 5e6 1e7 2e7 5e7 1e8)

matlab -nodisplay -nosplash -singleCompThread -r "smoothness_20190621_sweep_prince_func(3, ${SM_ARRAY[$SLURM_ARRAY_TASK_ID]}, 20, 0.1, 2e4);exit;"

# remove workspace
rm -rf /scratch/ch153/wjl/slurmJobs/$slurmArrayID
