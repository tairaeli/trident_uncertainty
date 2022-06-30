#!/bin/bash --login

# A template SLURM script for running CIAOLoop jobs in parallel

#SBATCH -N 2
#SBATCH --ntasks-per-node 100
#SBATCH --mem=100G
#SBATCH -t 04:00:00
#SBATCH --mail-user=tairaeli@msu.edu
#SBATCH --mail-type=ALL

# Do any job preamble stuff here; change directories, load modules, etc
cd /mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/parallel_tools
# Create the machine list by parsing $SLURM_NODELIST
# Default filename is "machines.dat"
/mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/parallel_tools/make_machine_list_slurm.pl

# Run CIAOLoop in parallel
/mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/CIAOLoop -m /mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/parallel_tools/machines.dat /mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/gen_uvb/uvb_params.par
