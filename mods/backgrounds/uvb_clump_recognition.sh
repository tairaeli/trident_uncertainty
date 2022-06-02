#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=tairaeli@msu.edu
#SBATCH --mail-type=ALL
#SBATCH -o /mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/err.txt

########## Command Lines for Job Running ##########

DS=$1
NRAYS=$2
UV_BACK=$3
CONFIG_LOC=$4
ABUN=$5

cd /mnt/home/tairaeli/trident_uncertainty/mods/backgrounds

module load Python/3.6.4

python make_rays.py --ds $DS --nrays $NRAYS #-- halo $HALO

python define_uvb.py --background_file $UV_BACK --config_loc $CONFIG_LOC

python make_tables.py --ds $DS --nrays $NRAYS --abun $ABUN #--halo $HALO