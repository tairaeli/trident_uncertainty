#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=2G
#SBATCH -A galaxies
#SBATCH --job-name cutoff_checking

########## Command Lines for Job Running ##########

conda activate astro_env

nrays=4
ions=ion_list.txt

arg_1=cutoff_frac
start_1=0.6
end_1=0.9
stepsize_1=0.05

arg_2=min_dens
start_2=11
end_2=15.5
stepsize_2=0.5

path_1=/mnt/scratch/tairaeli/uvb_dat/fg2009_ss_hr.h5
name_1=FG_2009

path_2=/mnt/scratch/tairaeli/trident_inputs/fg_test.h5
name_2=FG_2020

out_dir=/mnt/scratch/tairaeli/cutoff_bin_fg2009_fg2020

for i in $(seq $start_1 $stepsize_1 $end_1);
do
    echo "current iteration:"
    echo $i

    python SALSA_dens_cutoffs.py -$arg_1 $i -var_arg $arg_1 -ion_list $ions -nrays $nrays -uvb_path $path_1 -uvb_name $name_1 -out_dir $out_dir

    python SALSA_dens_cutoffs.py -$arg_1 $i -var_arg $arg_1 -ion_list $ions -nrays $nrays -uvb_path $path_2 -uvb_name $name_2 -out_dir $out_dir

    python SALSA_dens_cutoffs.py -$arg_1 $i -var_arg $arg_1 -ion_list $ions -nrays $nrays -make_plot True -out_dir $out_dir -uvb_name $name_1 -uvb_name2 $name_2

done

for j in $(seq $start_2 $stepsize_2 $end_2);
do
    echo "current iteration:"
    echo $j

    python SALSA_dens_cutoffs.py -$arg_2 $j -var_arg $arg_2 -ion_list $ions -nrays $nrays -uvb_path $path_1 -uvb_name $name_1 -out_dir $out_dir

    python SALSA_dens_cutoffs.py -$arg_2 $j -var_arg $arg_2 -ion_list $ions -nrays $nrays -uvb_path $path_2 -uvb_name $name_2 -out_dir $out_dir

    python SALSA_dens_cutoffs.py -$arg_2 $j -var_arg $arg_2 -ion_list $ions -nrays $nrays -make_plot True -out_dir $out_dir

done

python count_lonely_clumps.py -ion_list $ions -nrays $nrays