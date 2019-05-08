#!/bin/bash
#SBATCH --job-name=hw4_a1_ruben-vazquez_$1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ruben.vazquez@ufl.edu
#SBATCH --account=eel6763
#SBATCH --qos=eel6763
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100mb
#SBATCH -t 00:05:00
#SBATCH -o out_hw4_a1/hw4_a1_ruben-vazquez_$1.out
#SBATCH -e err_hw4_a1/hw4_a1_ruben-vazquez_$1.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tesla:1

module purge
module load ufrc intel/2018 cuda/9.2.88

srun nvprof ./cuda $1
