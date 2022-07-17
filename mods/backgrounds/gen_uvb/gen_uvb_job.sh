#!/bin/bash --login

# A template SLURM script for running CIAOLoop jobs in parallel

#SBATCH -N 4
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 128
#SBATCH --mem=200G
#SBATCH -t 02:00:00
#SBATCH --mail-user=tairaeli@msu.edu
#SBATCH --mail-type=ALL

# Do any job preamble stuff here; change directories, load modules, etc
cd /mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/gen_uvb

# Run CIAOLoop in parallel
srun python ./full_cloudy.py