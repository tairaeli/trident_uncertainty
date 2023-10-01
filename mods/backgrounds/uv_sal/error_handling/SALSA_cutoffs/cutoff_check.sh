#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=2G
#SBATCH --account galaxies
#SBATCH --job-name abun_job_submission1

########## Command Lines for Job Running ##########

start_cutoff=10

pcw_2019=/mnt/home/tairaeli/trident_uncertainty/mods/abundances/data_bin/par_test.h5
pcw_name=PCW_2019

hm_2012=/mnt/home/tairaeli/trident_uncertainty/mods/abundances/data_bin/hm2012_ss_hr.h5
hm_name=HM_2012

for i in $(seq $start_cutoff 0.05 13);
do
    echo "current iteration:"
    echo $i

    python SALSA_dens_cutoffs.py -min_dens $i -uvb_path $pcw_2019 -uvb_name $pcw_name

    python SALSA_dens_cutoffs.py -min_dens $i -uvb_path $hm_2012 -uvb_name $hm_name

    python SALSA_dens_cutoffs.py -min_dens $i -make_plot True

done