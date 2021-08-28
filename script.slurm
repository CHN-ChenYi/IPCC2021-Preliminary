#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -p amd_256
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun --bind-to none ./SLIC $CASE
