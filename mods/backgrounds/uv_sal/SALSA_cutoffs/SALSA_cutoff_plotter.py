# importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sn
import scipy as sc
import pickle
import yt
from salsa.utils import check_rays
import argparse
import os

parser = argparse.ArgumentParser(description = "Select cutoff and UVB for analysis")

parser.add_argument("-var_arg", action='store',
                    required=True, dest="var_arg",
                    help="argument that is varied between each iteration.\
                    Can either be 'cutoff_frac' or 'min_dens'.")

parser.add_argument("-iter", action='store',
                    required=True, dest="iter_list", type=str,
                    help="shows the start, stop and stepsize of arg\
                    iteration respectively")

parser.add_argument('-ion_list', action='store', 
                    required=False, dest='ion_list', 
                    help='Ions to be analyzed')

parser.add_argument('-nrays', action='store', 
                    required=False, dest='nrays', type = int,
                    help='Number of rays to be used for analysis')

parser.add_argument('-out_dir', action='store', 
                    required=False, dest='out_dir', 
                    help='Directory where files are output')

parser.add_argument('-uvb_path', action='store', 
                    required=False, dest='uvb', 
                    help='Path to UVB file')

parser.add_argument('-uvb_name', action='store', 
                    required=False, dest='uvb_name', 
                    help='Label to assign to uvb.')

args = parser.parse_args()
dic_args = vars(args)

# initializing ion list
ion_list = pd.read_csv(args.ion_list, sep = " ", header=0)

# preparing iteration values
iter_vals_cf = np.arange(0.6,0.9,0.05)
iter_vals_md = np.arange(11,15.5,0.5)

name1="FG_2009"

plot_df_cf = {}
plot_df_md = {}

for ion in ion_list:
    plot_df_cf[ion] = {}
    plot_df_md[ion] = {}
    for cut in iter_vals_cf:
        plot_df_cf[ion][cut] = np.array([])
    
    for cut in iter_vals_md:
        plot_df_md[ion][cut] = np.array([])

n = len(str(args.nrays))
for ray in range(args.nrays):

    out_path_ray = args.out_dir+f"/ray_{ray:0{n}d}/"

    ray_filename = f'{args.out_dir}/rays/ray{ray:0{n}d}.h5'

    ray_dat = yt.load(ray_filename)

    for ion in ion_list:
        
        nion = ion.replace("_"," ")

        out_path_ion = args.out_dir+"/"+str(ion)+"/"
        
        # loading in cutoff fraction data
        for i,cut in enumerate(iter_vals_cf):

            try:
                with open(out_path_ion+f"density/{name1}/{name1}_dens_cutoff_frac_"+str(np.round(cut,2))+".pickle", "rb") as dens_dat:
                    uvb1_dens = pickle.load(dens_dat)

                with open(out_path_ion+f"clump_data/{name1}/{name1}_clump_dat_cutoff_frac_"+str(np.round(cut,2))+".pickle", "rb") as salsa_dat:
                    uvb1_clump_dat = pickle.load(salsa_dat)

            except AttributeError:
                print("One or more settings found no absorbers. Skipping")
                pass

            ray_pos = ray_dat.r[("gas","l")].to("kpc")

            try:
                uvb1_clumps = uvb1_clump_dat[nion][name1][uvb1_clump_dat[nion][name1]["lightray_index"] == f"{ray:0{n}d}"]
            
            except TypeError:
                print("One or more settings found no absorbers. Skipping")
                pass

            plot_df_cf[ion][cut] = np.concatenate((plot_df_cf[ion][cut],uvb1_clumps["col_dens"]), axis = None)
        
        # loading in minimum density data   
        for i,cut in enumerate(iter_vals_md):

            try:
                with open(out_path_ion+f"density/{name1}/{name1}_dens_min_dens_"+str(np.round(cut,2))+".pickle", "rb") as dens_dat:
                    uvb1_dens = pickle.load(dens_dat)

                with open(out_path_ion+f"clump_data/{name1}/{name1}_clump_dat_min_dens_"+str(np.round(cut,2))+".pickle", "rb") as salsa_dat:
                    uvb1_clump_dat = pickle.load(salsa_dat)

            except AttributeError:
                print("One or more settings found no absorbers. Skipping")
                pass

            ray_pos = ray_dat.r[("gas","l")].to("kpc")

            try:
                uvb1_clumps = uvb1_clump_dat[nion][name1][uvb1_clump_dat[nion][name1]["lightray_index"] == f"{ray:0{n}d}"]
            
            except TypeError:
                print("One or more settings found no absorbers. Skipping")
                pass

            plot_df_md[ion][cut] = np.concatenate((plot_df_md[ion][cut],uvb1_clumps["col_dens"]), axis = None)


print("Data loaded, generating plots...")

n_iter_cf = len(iter_vals_cf)
colors_cf = plt.cm.Set2(np.linspace(0,1,int(n_iter_cf)+1))
n_iter_md = len(iter_vals_md)
colors_md = plt.cm.Set2(np.linspace(0,1,int(n_iter_md)+1))

# for ion in ion_list:
ion = 'H_I'
nion = ion.replace("_"," ")

plot_path = args.out_dir+"/"+str(ion)+"/"
if os.path.exists(plot_path) == False:
    os.mkdir(plot_path)



# determining resolution of figure
res = 100

hist_dat_cf = []
max_val_cf = -np.inf
for i,cut in enumerate(iter_vals_cf):
    temp_hist = np.histogram(plot_df_cf[ion][cut], bins=20)
    hist_dat_cf.append(temp_hist)

    if np.max(temp_hist[0]) > max_val_cf:
        max_val_cf = np.max(temp_hist[0])

hist_dat_md = []
max_val_md = -np.inf
for i,cut in enumerate(iter_vals_md):
    temp_hist = np.histogram(plot_df_md[ion][cut], bins=20)
    hist_dat_md.append(temp_hist)

    if np.max(temp_hist[0]) > max_val_md:
        max_val_md = np.max(temp_hist[0])

fig, ax = plt.subplots(1,2, figsize=(15,8), sharex=True, sharey=True,
                       constrained_layout=True)

col_list = []
for i,cut in enumerate(iter_vals_cf):
    hist = hist_dat_cf[i]
    bheight = np.hstack((hist[0][0], np.array(hist[0])))
    sn.lineplot(x=hist[1], y=bheight/max_val_cf, drawstyle='steps-pre', color=colors_cf[i],
                linewidth=3, ax = ax[0])
    col_list.append(mlines.Line2D([], [], color=colors_cf[i],
                    markersize=10, label="cutoff frac"+" = "+str(np.round(cut,2))))

ax[0].legend(handles=col_list, prop={'size': 15})
ax[0].set_title("(a)", fontsize = 15)
ax[0].grid()

col_list = []
for i,cut in enumerate(iter_vals_md):
    hist = hist_dat_md[i]
    bheight = np.hstack((hist[0][0], np.array(hist[0])))
    sn.lineplot(x=hist[1], y=bheight/max_val_md, drawstyle='steps-pre', color=colors_md[i],
                linewidth=3, ax = ax[1])
    col_list.append(mlines.Line2D([], [], color=colors_md[i],
                    markersize=10, label="min dens"+" = "+str(np.round(cut,2))))

ax[1].legend(handles=col_list, prop={'size': 15})
ax[1].set_title("(b)", fontsize = 15)
ax[1].grid()
fig.supxlabel(r"$\log\sigma_{FG09}$ ($cm^{-2}$)", fontsize = 20)
fig.supylabel(r"% max($\sigma_{FG09}$)", fontsize = 20)
plt.savefig(plot_path+"UVB_dens_compare_"+args.var_arg+"_test.pdf")
plt.show()

print("Run Complete!")