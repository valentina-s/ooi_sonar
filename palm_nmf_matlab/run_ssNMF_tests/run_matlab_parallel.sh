#!/bin/bash
#
#SBATCH --job-name=Matlab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=15GB
#SBATCH --time=12:00:00
 
module purge
module load matlab/2018b
 
cd /scratch/ch153/wjl/matlab-prefdir
 
if [ "$SLURM_JOBTMP" == "" ]; then
    export SLURM_JOBTMP=/state/partition1/$USER/$$
    mkdir -p $SLURM_JOBTMP
fi
 
export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)
 
echo
echo "Hostname: $(hostname)"
echo

cd /home/ch153/wjl/ooi_sonar/palm_nmf_matlab
 
cat<<EOF | srun matlab -nodisplay
parpool('local', 10)
smoothness_20190621_sweep_prince
exit
EOF
 
rm -rf $SLURM_JOBTMP/*
