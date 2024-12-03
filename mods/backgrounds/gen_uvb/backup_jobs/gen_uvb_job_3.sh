#!/bin/bash --login

# A template SLURM script for running CIAOLoop jobs in parallel

#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 128
#SBATCH --mem=200G
#SBATCH -t 30:00:00
#SBATCH -A galaxies
#SBATCH --mail-user=tairaeli@msu.edu
#SBATCH --mail-type=ALL

# Do any job preamble stuff here; change directories, load modules, etc
cd /mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/

# Run CIAOLoop in parallel
./CIAOLoop -mp 3 4 -np 128 /mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/gen_uvb/backup_jobs/uvb_params.par