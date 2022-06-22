#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=2G
#SBATCH --account galaxies
#SBATCH --job-name abun_job_submission1

########## Command Lines for Job Running ##########

cd /mnt/home/f0104093/trident_uncertainty/mods/abundances/scripts
module load python/3.6.4
srun -n 128 python sal_the_snek.py --ds /mnt/scratch/f0104093/cgm_abundance_variance/  --nrays 4 --abun /mnt/home/f0104093/new_abundances/cgm_abundances.txt
pyhton id_clumps_test.py
python clumps_stats.py

#scontrol show job $SLURM_JOB_ID
#js -j $SLURM_JOB_ID
