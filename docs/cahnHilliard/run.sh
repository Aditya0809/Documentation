#!/bin/bash


#SBATCH -t 04:00:00
#SBATCH -n 128
#SBATCH -o "%x.o%j"
#SBATCH -e "%x.e%j"
#SBATCH --job-name="ch2d_explicit"
#SBATCH --mem-per-cpu=1G
#SBATCH -p RM


# Run the main program
mpirun -np 128 ./ch2d_explicit > "${SLURM_JOB_NAME}.o$id"
