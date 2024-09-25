import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import pickle
import configparser
import os
import h5py
import pandas as pd
import roman

def get_true_rs(val): ##define how to get actual rshift numbers
    if val == 20:
        true_rs = '2.0'
    elif val == 18:
        true_rs = '2.5'
    return true_rs

# reading in arguments
sal_args = configparser.ConfigParser()
sal_args.read("/mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/uv_sal/pipeline/sal_params.par")

# set desired halo pattern
halo = sal_args["galaxy_settings"]["gal_pattern"]
rs = sal_args["galaxy_settings"]["redshift"]
nrays = int(sal_args["galaxy_settings"]["nrays"])

# gets the true rs needed
true_rs = get_true_rs(int(rs))

# identifying the path argument as a variable
out_file = sal_args["base_settings"]["output_file"]
path = os.path.expandvars(os.path.expanduser(out_file))
halo_path = path+'/halo'+f'{halo}'
rs_path = halo_path + '/redshift'+f'{true_rs}'

# initializing UVB data
uvb_names = sal_args["uvb_analysis"]["uvb_names"].split(" ")
uvb_filenames = sal_args["uvb_analysis"]["uvb_filenames"].split(" ")
sal_dat = None

rs_ids = [11,4,11,4]

old_gen = [uvb_names[0], uvb_names[2], uvb_names[1]]
old_fname = [uvb_filenames[0], uvb_filenames[2], uvb_filenames[1]]
old_rs = [rs_ids[0], rs_ids[2], rs_ids[1]]

new_gen = [uvb_names[1], uvb_names[3], uvb_names[3]]
new_fname = [uvb_filenames[1], uvb_filenames[3], uvb_filenames[3]]
new_rs = [rs_ids[1], rs_ids[3], rs_ids[3]]

ion_list = sal_args["galaxy_settings"]["ions"].split(" ")
num_ion = len(ion_list)

ax_lab_size = 15
title_size = 15

lb_cutoff = -27
ub_cutoff = -0.00001

# Parameter1 = column density: column 2
# Parameter2 = redshift: column 3
# Temperature = column 4
comp_paths = {}
for i in range(len(old_gen)):
    print(f"Comparing {new_gen[i]} and {old_gen[i]}") 
    for ion in ion_list:
        print("Ion: "+ion)
        fig, ax = plt.subplots(2,2, figsize=(10,10))
        # colors = sns.color_palette("rocket",3)
        w_size = 10

        atom, istate = ion.split("_")
        istate_num = roman.fromRoman(istate)
        
        old_f = h5py.File(old_fname[i],'r')
        old_col_dens = old_f[atom].attrs["Parameter1"]
        old_temp = old_f[atom].attrs["Temperature"]
        old_f_filter = old_f[atom][istate_num,:,old_rs[i],:]
        old_f_filter[old_f_filter < lb_cutoff] = np.nan
        old_f_filter[old_f_filter > ub_cutoff] = np.nan
        old_f.close()
        
        new_f = h5py.File(new_fname[i],'r')
        new_col_dens = new_f[atom].attrs["Parameter1"]
        new_temp = new_f[atom].attrs["Temperature"]
        new_f_filter = new_f[atom][istate_num,:,new_rs[i],:]
        new_f_filter[new_f_filter < lb_cutoff] = np.nan
        new_f_filter[new_f_filter > ub_cutoff] = np.nan
        new_f.close()

        # print(old_f_filter.attrs)
        comp_paths[f"{old_gen[i]}_{new_gen[i]}"] = rs_path+f"/{old_gen[i]}_{new_gen[i]}_comp"
        if os.path.exists(comp_paths[f"{old_gen[i]}_{new_gen[i]}"]) == False:
            os.mkdir(comp_paths[f"{old_gen[i]}_{new_gen[i]}"])

        ion_path = comp_paths[f"{old_gen[i]}_{new_gen[i]}"]+"/"+ion
        if os.path.exists(ion_path) == False:
            os.mkdir(ion_path)

        if new_f_filter.shape != old_f_filter.shape:
            new_f_filter = new_f_filter[0:old_f_filter.shape[0]]
            new_col_dens = new_col_dens[0:len(old_col_dens)]

        f1 = ax[0,0].pcolormesh(old_col_dens, old_temp, old_f_filter.T)
        ax[0,0].set_ylabel("$log_{10}$ Temperature (K)", fontsize=ax_lab_size)
        ax[0,0].set_xlabel(r"$log_{10}$(n) ($cm^{-3}$)", fontsize=ax_lab_size)
        ax[0,0].set_title(f"{ion} Fraction {old_gen[i]}", fontsize=title_size)
        fig.colorbar(f1,ax = ax[0,0])

        f2 = ax[1,1].pcolormesh(new_col_dens, new_temp, new_f_filter.T)
        ax[1,1].set_ylabel("$log_{10}$ Temperature (K)", fontsize=ax_lab_size)
        ax[1,1].set_xlabel(r"$log_{10}$(n) ($cm^{-3}$)", fontsize=ax_lab_size)
        ax[1,1].set_title(f"{ion} Fraction {new_gen[i]}", fontsize=title_size)
        fig.colorbar(f2,ax = ax[1,1])

        f3 = ax[1,0].pcolormesh(old_col_dens, old_temp, (new_f_filter - old_f_filter).T)
        ax[1,0].set_ylabel("$log_{10}$ Temperature (K)", fontsize=ax_lab_size)
        ax[1,0].set_xlabel(r"$log_{10}$(n) ($cm^{-3}$)", fontsize=ax_lab_size)
        ax[1,0].set_title(f"{ion} Fraction {new_gen[i]}/{old_gen[i]}", fontsize=title_size)
        fig.colorbar(f3,ax = ax[1,0])
        
        ax[0,1].axis("off")

        plt.tight_layout()
        plt.savefig(rs_path+f"/{old_gen[i]}_{new_gen[i]}_comp/{ion}/ion_frac_{old_gen[i]}_{new_gen[i]}_{ion}.pdf")
        plt.clf()



plt.savefig(ion_path+f"/phys_quantities_{old_gen[i]}_{new_gen[i]}_{ion}.pdf")

