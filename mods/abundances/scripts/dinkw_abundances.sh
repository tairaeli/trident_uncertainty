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
srun -n 128 python sal_the_snek.py --ds /mnt/scratch/f0104093/condensed_pipeline_tests/  --nrays 4 --abun /mnt/home/f0104093/new_abundances/cgm_abundances.txt

for ION in C_II C_IV O_VI
do
    ls /mnt/scratch/f0104093/condensed_pipeline_tests/data/*${ION}.txt > temp_filelist.txt
    python plot_clumps.py --ds /mnt/scratch/f0104093/condensed_pipeline_tests/visuals/ --fn temp_filelist.txt
    cat /dev/null > temp_filelist.txt
done
echo "All done :)"

pyhton id_clumps_test.py
python make_hist_new.py

#scontrol show job $SLURM_JOB_ID
#js -j $SLURM_JOB_ID
