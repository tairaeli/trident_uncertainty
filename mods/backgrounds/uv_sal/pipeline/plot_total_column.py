"""
Makes plot of total column densities of selected ions for 
different ultraviolet backgrounds
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import seaborn as sns
import pandas as pd
import pickle
import configparser
import os
import yt

def get_true_rs(val): ##define how to get actual rshift numbers
    if val == 20:
        true_rs = '2.0'
    elif val == 18:
        true_rs = '2.5'
    return true_rs

def lonely_hunter(uvb1,uvb2):
    """
    Isolates instances where there exist lonely clumps of gas

    args:

        uvb1 (dataframe) - data on the cloud for a given ray and ion
        uvb2 (dataframe) - see above variable
    
    returns:

        lone1 (dataframe) - clump information on all lonely gas clumps in array
        lone2 (dataframe) - see above variable

        red_uvb1 (dataframe) - clump information on all non-lonely gas clumps in array
        red_uvb2 (dataframe) - see above variable
    """

    uvb1_mask = np.where(np.array(uvb1["col_dens"])==0)
    # preparing data containers
    red_uvb1 = {}
    lone1 = {}
    red_uvb2 = {}
    lone2 = {}

    # removing lonely clumps from each key in dict
    for key in uvb1.keys():
        # removing 0 vals in uvb1
        red_uvb1[key] = np.delete(np.array(uvb1[key]),uvb1_mask)

        # isolating lonely clumps in uvb2
        lone2[key] = np.array(uvb2[key])[uvb1_mask]

        # removing lonely clumps from uvb2
        red_uvb2[key] = np.delete(np.array(uvb2[key]), uvb1_mask)

    assert len(uvb1["col_dens"]) >= len(red_uvb1["col_dens"]), "first masking failed (uvb1)"
    assert len(uvb2["col_dens"]) >= len(red_uvb2["col_dens"]), "first masking failed (uvb2)"

    # making second mask with previously masked data
    uvb2_mask = np.where(np.array(red_uvb2["col_dens"])==0)

    # removing lonely clumps from each key in dict
    for key in uvb1.keys():
        # removing 0 vals in uvb2
        red_uvb2[key] = np.delete(red_uvb2[key], uvb2_mask)

        # isolating lonely clumps in uvb1
        lone1[key] = red_uvb1[key][uvb2_mask]

        # removing lonely clumps from uvb1
        red_uvb1[key] = np.delete(red_uvb1[key], uvb2_mask)

    assert len(red_uvb1["col_dens"]) == len(red_uvb2["col_dens"]), "second masking failed"

    return lone1, lone2, red_uvb1, red_uvb2

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

sal_dat = None

old_gen = [uvb_names[0], uvb_names[2], uvb_names[1]]#, uvb_names[0], uvb_names[2]]
new_gen = [uvb_names[1], uvb_names[3], uvb_names[3]]#, uvb_names[3], uvb_names[1]]

short_uvb_names = {"FG_2009":"FG09", "FG_2020":"FG20",
                   "HM_2012":"HM12", "PCW_2019":"PW19"}

for i,name in enumerate(uvb_names):
    dat_path = rs_path +f'/{name}/data'

    if i == 0:
        with open(f'{dat_path}/salsa_out_dict.pickle',"rb") as dat:
            sal_dat = pickle.load(dat)

    else:
        with open(f'{dat_path}/salsa_out_dict.pickle',"rb") as dat:
            uvb_dict = pickle.load(dat)

            for ion in uvb_dict.keys():
        
                sal_dat[ion][name] = uvb_dict[ion][name]

# loading in data
ion_list = sal_args["galaxy_settings"]["ions"].split(" ")

# creating dictionaries to store our data if they don't already exist
uvb_dist_path = rs_path+"/uvb_dists"
if os.path.exists(uvb_dist_path) == False:
    os.mkdir(uvb_dist_path)

# creating column density plots for each UVB
num_ion = len(ion_list)

# reformatting ion labels
ion_name_dict = {"H_I":r"H $\mathrm{\i}$",
                 "Si_II":r"Si $\mathrm{\i\i}$",
                 "Si_III":r"Si $\mathrm{\i\i\i}$",
                 "C_III":r"C $\mathrm{\i\i\i}$",
                 "Si_IV":r"Si $\mathrm{\i v}$",
                 "C_IV":r"C $\mathrm{\i v}$",
                 "N_V":r"N $\mathrm{v}$",
                 "O_VI":r"O $\mathrm{v\i}$"}

# setting up total column density figure
fig, ax = plt.subplots(2, (num_ion//2)+(num_ion%2), figsize=(25,10),
                       constrained_layout = True)
                       
colors = sns.color_palette("rocket",len(old_gen))
w_size = 30
ion_txt_size = 40
j=0
k=0

for ion in ion_list:    
    tot_ray_dens = {}
    n = len(str(nrays))

    for i in range(len(uvb_names)):
        tot_ray_dens[uvb_names[i]] = np.zeros(nrays)

    for ray in range(nrays):
        # loading in ray data
        n = len(str(nrays))

        if ion == "H_I":
            ray_dat = yt.load(rs_path+f"/rays_HI/ray{ray:0{n}d}.h5")
            ray_len = np.max(ray_dat.r[("gas","dl")].to("cm"))
        
        else:
            ray_dat = yt.load(rs_path+f"/rays/ray{ray:0{n}d}.h5")
            ray_len = np.max(ray_dat.r[("gas","dl")].to("cm"))

        for i in range(len(uvb_names)):
            if ion == "H_I":
                with open(rs_path +f'/{uvb_names[i]}/ray_dat/{ion}/ray_{ray:0{n}d}_dens_HI.pickle', "rb") as dens_dat:
                    num_dens = pickle.load(dens_dat)
                tot_ray_dens[uvb_names[i]][ray] = np.sum(num_dens*ray_len)
            else:
                with open(rs_path +f'/{uvb_names[i]}/ray_dat/{ion}/ray_{ray:0{n}d}_dens.pickle', "rb") as dens_dat:
                    num_dens = pickle.load(dens_dat)
                tot_ray_dens[uvb_names[i]][ray] = np.sum(num_dens*ray_len)

    for i in range(len(old_gen)):
        avg_dens = np.log10((tot_ray_dens[old_gen[i]] + tot_ray_dens[new_gen[i]])/2)
        dens_diff = np.log10(tot_ray_dens[old_gen[i]]) - np.log10(tot_ray_dens[new_gen[i]])
        
        moving_avg = np.convolve(dens_diff, np.ones(w_size) / w_size, mode='valid')
        x_vals = avg_dens[w_size-1:]
        
        x_vals = np.sort(x_vals)
        moving_avg = moving_avg[x_vals.argsort()]
        
        # adding ion label
        if i == 0: 
            
            element,istate = ion.split("_")
            
            if (ion =="N_V") or (ion=="O_VI"):
                ax[j,k].text(0.2, 0.2, 
                      s=ion_name_dict[ion], fontsize=ion_txt_size,
                      transform=ax[j,k].transAxes)
            else:
                ax[j,k].text(0.7, 0.2, 
                      s=ion_name_dict[ion], fontsize=ion_txt_size,
                      transform=ax[j,k].transAxes)
            
        ax[j,k].scatter(x = avg_dens, y = dens_diff, 
                        label=f"{short_uvb_names[old_gen[i]]}/{short_uvb_names[new_gen[i]]}",
                        color = colors[i])
        # ax[j,k].plot(x_vals, moving_avg, color = colors[i])
        if (ion == "H_I") or (ion == "Si_II"):
            ax[j,k].set_yticks([1,0,-1,-2,-3]) 
        elif (ion == "C_IV"):
            ax[j,k].set_yticks([0.1,0,-0.1,-0.2,-0.3,-0.4])
        elif (ion == "N_V"):
            ax[j,k].set_yticks([0.2,0,-0.2,-0.4,-0.6])
        elif (ion == "O_VI"):
            ax[j,k].set_yticks([0.2,0,-0.2,-0.4,-0.6,-0.8,-1.0])
        
    ax[j,k].grid()
    ax[j,k].tick_params(axis='both', labelsize=20)

    if (k>=(num_ion//2)+(num_ion%2)-1):
        j+=1
        k=0
    else:
        k+=1

ax[0,0].legend(loc='lower left', prop={'size': 15})
fig.supxlabel(r"log($\overline{N_{LOS}}$) [$cm^{-2}$]", fontsize=35)
fig.supylabel(r"log($N_{LOS,1}/N_{LOS,2}$) [$cm^{-2}$]", fontsize=35)
plt.savefig(uvb_dist_path+f"/tot_col_dens_comp.png",
                dpi=400,bbox_inches='tight')
plt.clf()