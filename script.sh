#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -p amd_256
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1

mpirun ./SLIC
