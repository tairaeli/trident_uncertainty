#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=5G
#SBATCH --account galaxies
#SBATCH --job-name abun_job_submission1

########## Command Lines for Job Running ##########

module load Conda/3

conda activate astro_env

path1=/mnt/scratch/tairaeli/trident_inputs/fg2009_ss_hr.h5
name1=FG_2009

path2=/mnt/scratch/tairaeli/trident_inputs/fg_test.h5
name2=FG_2020

path3=/mnt/scratch/tairaeli/trident_inputs/hm2012_ss_hr.h5
name3=HM_2012

path4=/mnt/scratch/tairaeli/trident_inputs/pcw_test.h5
name4=PCW_2019

# running salsa
python absorber_extractor.py -uvb_path $path1 -uvb_name $name1
python absorber_extractor_H_I.py -uvb_path $path1 -uvb_name $name1

python absorber_extractor.py -uvb_path $path2 -uvb_name $name2
python absorber_extractor_H_I.py -uvb_path $path2 -uvb_name $name2

python absorber_extractor.py -uvb_path $path3 -uvb_name $name3
python absorber_extractor_H_I.py -uvb_path $path3 -uvb_name $name3

python absorber_extractor.py -uvb_path $path4 -uvb_name $name4
python absorber_extractor_H_I.py -uvb_path $path4 -uvb_name $name4

# perfoming analysis
python uvb_abun_pairwise_compare.py -uvb_path1 $path1 -uvb_name1 $name1 -uvb_path2 $path2 -uvb_name2 $name2
python uvb_abun_pairwise_compare.py -uvb_path1 $path1 -uvb_name1 $name1 -uvb_path2 $path3 -uvb_name2 $name3
python uvb_abun_pairwise_compare.py -uvb_path1 $path4 -uvb_name1 $name4 -uvb_path2 $path2 -uvb_name2 $name2
python uvb_abun_pairwise_compare.py -uvb_path1 $path3 -uvb_name1 $name3 -uvb_path2 $path4 -uvb_name2 $name4

# making some cross comparisons
python uvb_abun_pairwise_compare.py -uvb_path1 $path1 -uvb_name1 $name1 -uvb_path2 $path4 -uvb_name2 $name4
python uvb_abun_pairwise_compare.py -uvb_path1 $path3 -uvb_name1 $name3 -uvb_path2 $path2 -uvb_name2 $name2

# generating plots 
python plot_total_column.py
python plot_pair_comp.py
