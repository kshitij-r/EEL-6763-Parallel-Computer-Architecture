#!/bin/bash
#SBATCH --job-name=mat_mul
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kshitijraj@ufl.edu
#SBATCH --account=eel6763
#SBATCH --qos=eel6763
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100mb
#SBATCH -t 00:05:00
#SBATCH -o mat_mul.txt
#SBATCH -e mat_mul.err

module load intel
mpicc -n 2  hybrid_mat_mult.c -o hyb -fopenmp 
srun --mpi=pmix_v2 ./hyb 4