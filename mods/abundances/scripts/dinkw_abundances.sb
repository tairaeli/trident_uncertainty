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
module unload Python
which python
srun -n 128 python sal_the_snek.py --ds /mnt/scratch/f0104093/condensed_pipeline_tests/  --nrays 4 --abun /mnt/home/f0104093/new_abundances/new_abundances/cgm_abundances.txt

#scontrol show job $SLURM_JOB_ID
#js -j $SLURM_JOB_ID
