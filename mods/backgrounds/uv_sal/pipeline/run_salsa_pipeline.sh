#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=2G
#SBATCH --account galaxies
#SBATCH --job-name abun_job_submission1

########## Command Lines for Job Running ##########

module load Conda/3

conda activate astro_env

path1=/mnt/scratch/tairaeli/uvb_dat/fg2009_ss_hr.h5
name1=FG_2009

path2=/mnt/scratch/tairaeli/trident_inputs/fg_test.h5
name2=FG_2020

path3=/mnt/scratch/tairaeli/uvb_dat/hm2012_ss_hr.h5
name3=HM_2012

path4=/mnt/scratch/tairaeli/trident_inputs/pcw_test.h5
name4=PCW_2019

# running salsa
# python sal_the_super_uvb.py -uvb_path $path1 -uvb_name $name1
# python sal_the_super_uvb.py -uvb_path $path2 -uvb_name $name2
# python sal_the_super_uvb.py -uvb_path $path3 -uvb_name $name3
# python sal_the_super_uvb.py -uvb_path $path4 -uvb_name $name4

# perfoming analysis
python uvb_abun_pairwise_compare.py -uvb_path1 $path1 -uvb_name1 $name1 -uvb_path2 $path2 -uvb_name2 $name2
python uvb_abun_pairwise_compare.py -uvb_path1 $path1 -uvb_name1 $name1 -uvb_path2 $path3 -uvb_name2 $name3
python uvb_abun_pairwise_compare.py -uvb_path1 $path2 -uvb_name1 $name2 -uvb_path2 $path4 -uvb_name2 $name4
python uvb_abun_pairwise_compare.py -uvb_path1 $path3 -uvb_name1 $name3 -uvb_path2 $path4 -uvb_name2 $name4

# generating plots 
python pair_comp_plot_paper.py
