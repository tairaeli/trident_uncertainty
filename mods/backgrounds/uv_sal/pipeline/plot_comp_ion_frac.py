import numpy as np
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import configparser
import os
import h5py
import roman

def get_true_rs(val): ##define how to get actual rshift numbers
    if val == 20:
        true_rs = '2.0'
    elif val == 18:
        true_rs = '2.5'
    return true_rs

# reformatting ion labels
ion_name_dict = {"H_I":r"H $\mathrm{\i}$",
                 "Si_II":r"Si $\mathrm{\i\i}$",
                 "Si_III":r"Si $\mathrm{\i\i\i}$",
                 "C_III":r"C $\mathrm{\i\i\i}$",
                 "Si_IV":r"Si $\mathrm{\i v}$",
                 "C_IV":r"C $\mathrm{\i v}$",
                 "N_V":r"N $\mathrm{v}$",
                 "O_VI":r"O $\mathrm{v\i}$"}

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

# contains shortened UVB names
short_uvb_names = {"FG_2009":"FG09", "FG_2020":"FG20",
                   "HM_2012":"HM12", "PCW_2019":"PW19"}

# defining ion list and nunber of ions
ion_list = ["H_I", "O_VI"]
# can also take in list of arguments from parameter file
# ion_list = sal_args["galaxy_settings"]["ions"].split(" ")
num_ion = len(ion_list)

# defining size of axis labels
ax_lab_size = 20
title_size = 20

# defining bounds of colorbar (in logscale)
lb_cutoff = -6 # 10**-30 # 10**-3
ub_cutoff = 0 #-0.00001

# Parameter1 = gas density: column 2
# Parameter2 = redshift: column 3
# Temperature = column 4

comp_paths = {}

# iterating through each parwise comparison
for i in range(len(old_gen)):
    print(f"Comparing Ion Fractions {new_gen[i]} and {old_gen[i]}") 
    # iterating through each ion
    for ion in ion_list:
        print("Ion: "+ion)

        # creating figure
        fig, ax = plt.subplots(1,3, figsize=(15,5))
        cmap = cm.get_cmap('viridis')
        im = cm.ScalarMappable()
        w_size = 10

        atom, istate = ion.split("_")
        istate_num = roman.fromRoman(istate)
        
        # loading in "old" data and masking out data based on set bounds
        old_f = h5py.File(old_fname[i],'r')
        old_col_dens = old_f[atom].attrs["Parameter1"]
        old_temp = old_f[atom].attrs["Temperature"]
        old_f_filter = old_f[atom][istate_num-1,:,old_rs[i],:]
        old_f_filter[old_f_filter < lb_cutoff] = np.nan
        old_f_filter[old_f_filter > ub_cutoff] = np.nan
        old_f.close()
        
        # repeating process with "new" data
        new_f = h5py.File(new_fname[i],'r')
        new_col_dens = new_f[atom].attrs["Parameter1"]
        new_temp = new_f[atom].attrs["Temperature"]
        new_f_filter = new_f[atom][istate_num-1,:,new_rs[i],:]
        new_f_filter[new_f_filter < lb_cutoff] = np.nan
        new_f_filter[new_f_filter > ub_cutoff] = np.nan
        new_f.close()
        
        # loading comparison paths. creating dictionaries if they don't exist
        comp_paths[f"{old_gen[i]}_{new_gen[i]}"] = rs_path+f"/{old_gen[i]}_{new_gen[i]}_comp"
        if os.path.exists(comp_paths[f"{old_gen[i]}_{new_gen[i]}"]) == False:
            os.mkdir(comp_paths[f"{old_gen[i]}_{new_gen[i]}"])

        ion_path = comp_paths[f"{old_gen[i]}_{new_gen[i]}"]+"/"+ion
        if os.path.exists(ion_path) == False:
            os.mkdir(ion_path)

        # adjusting shapes of dictionaries if they do not match 
        if new_f_filter.shape != old_f_filter.shape:
            new_f_filter = new_f_filter[0:old_f_filter.shape[0]]
            new_col_dens = new_col_dens[0:len(old_col_dens)]

        # setting up the 'old' 2d histogram
        f1 = ax[0].pcolormesh(10**old_col_dens, 10**old_temp, old_f_filter.T, 
                                cmap = cmap)
        # setting up contours 
        ax[0].contour(10**old_col_dens, 10**old_temp, old_f_filter.T, 
                        colors = "black")

        ax[0].set_ylabel("$T$ [K]", fontsize=ax_lab_size)
        ax[0].set_title("(a)")
        ax[0].set_yscale("log")
        ax[0].set_xscale("log")

        ax[0].text(0.1,0.82,short_uvb_names[old_gen[i]], 
                fontsize=ax_lab_size+4, transform=ax[0].transAxes)

        # setting up 'new' 2d historgram
        f2 = ax[1].pcolormesh(10**new_col_dens, 10**new_temp, new_f_filter.T, 
                                cmap = cmap)
        # adding contours
        ax[1].contour(10**new_col_dens, 10**new_temp, new_f_filter.T, 
                        colors = "black")

        ax[1].set_title("(b)")
        ax[1].set_yscale("log")
        ax[1].set_xscale("log")
        cb1 = fig.colorbar(f1,ax = ax[1])

        ax[1].text(1.23, 0.2, 
                f'log({ion_name_dict[ion]}) Ion Fraction', fontsize=ax_lab_size,
                transform=ax[1].transAxes, rotation=270)
        
        ax[1].text(0.1,0.82,short_uvb_names[new_gen[i]], 
                fontsize=ax_lab_size+4, transform=ax[1].transAxes)
        
        # setting up ratio comparison 2d histogram
        f3 = ax[2].pcolormesh(10**old_col_dens, 10**old_temp, (new_f_filter - old_f_filter).T)
        ax[2].set_title("(c)")
        ax[2].set_yscale("log")
        ax[2].set_xscale("log")
        ax[2].set_ylabel("$T$ [K]", fontsize=ax_lab_size)
        cb2 = fig.colorbar(f3,ax = ax[2])

        ax[2].text(1.28, 0.1, 
                f'log({ion_name_dict[ion]}) Fraction Ratio', fontsize=ax_lab_size,
                transform=ax[2].transAxes, rotation=270)
        
        ax[2].text(0.1,0.82,short_uvb_names[new_gen[i]]+"/"+short_uvb_names[old_gen[i]],
                fontsize=ax_lab_size+4, transform=ax[2].transAxes)

        fig.supxlabel(r"n [$cm^{-3}$]", fontsize=ax_lab_size)

        plt.tight_layout()
        # saving figure
        plt.savefig(rs_path+f"/{old_gen[i]}_{new_gen[i]}_comp/{ion}/ion_frac_{old_gen[i]}_{new_gen[i]}_{ion}.png",
                    dpi=400,bbox_inches='tight')
        plt.clf()


