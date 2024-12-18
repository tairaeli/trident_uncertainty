#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=2G
#SBATCH -A galaxies
#SBATCH --job-name cutoff_checking

########## Command Lines for Job Running ##########

module load Conda/3
conda activate astro_env

nrays=100
ions=ion_list.txt

arg_1=cutoff_frac
start_1=0.6
end_1=0.9
stepsize_1=0.05

iter_vals_1=($start_1,$end_1,$stepsize_1)

arg_2=min_dens
start_2=11.0
end_2=15.5
stepsize_2=0.5

iter_vals_2=($start_2,$end_2,$stepsize_2)

path=/mnt/scratch/tairaeli/trident_inputs/fg2009_ss_hr.h5
name=FG_2009

out_dir=/mnt/scratch/tairaeli/cutoff_bin_fg2009

for i in $(seq $start_1 $stepsize_1 $end_1);
do
    echo "current iteration:"
    echo $i

    python SALSA_cutoff_finder.py -$arg_1 $i -var_arg $arg_1 -ion_list $ions -nrays $nrays -uvb_path $path -uvb_name $name -out_dir $out_dir

done

# plotting cutoff fraction stuff
python SALSA_cutoff_plotter.py -var_arg $arg_1 -iter $iter_vals_1 -ion_list $ions -nrays $nrays -uvb_name $name -out_dir $out_dir

for j in $(seq $start_2 $stepsize_2 $end_2);
do
    echo "current iteration:"
    echo $j

    python SALSA_cutoff_finder.py -$arg_2 $j -var_arg $arg_2 -ion_list $ions -nrays $nrays -uvb_path $path -uvb_name $name -out_dir $out_dir

done

# plotting min density stuff
python SALSA_cutoff_plotter.py -var_arg $arg_2 -iter $iter_vals_2 -ion_list $ions -nrays $nrays -uvb_name $name -out_dir $out_dir