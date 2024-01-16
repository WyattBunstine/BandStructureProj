#!/bin/bash -l
#SBATCH --job-name=elk
#SBATCH --time=72:00:00
# Can use either shared or parallel queue since we use a whole node
#SBATCH --ntasks-per-node=48
#SBATCH --nodes=2
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=wbunsti1@jhu.edu
# MARCC defaults to gcc, must load intel module
module load intel/2022.2
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=24
export OMP_NESTED=false
export OMP_DYNAMIC=false
export OMP_STACKSIZE=1G
ulimit -Ss unlimited
/home/wbunsti1/elk/elk-8.8.26/src/elk >& elk.log
