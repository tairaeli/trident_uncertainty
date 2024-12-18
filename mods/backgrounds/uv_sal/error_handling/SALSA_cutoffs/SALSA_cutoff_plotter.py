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

name1 = args.uvb_name

# initializing ion list
ion_list = pd.read_csv(args.ion_list, sep = " ", header=0)

# preparing iteration values
iter_vals = args.iter_list.split(',')
iter_vals = [float(i) for i in iter_vals]

plot_df = {}

for ion in ion_list:
    plot_df[ion] = {}
    for cut in np.arange(iter_vals[0],iter_vals[1],iter_vals[2]):
        plot_df[ion][cut] = np.array([])

n = len(str(args.nrays))
for ray in range(args.nrays):

    out_path_ray = args.out_dir+f"/ray_{ray:0{n}d}/"

    ray_filename = f'{args.out_dir}/rays/ray{ray:0{n}d}.h5'

    ray_dat = yt.load(ray_filename)

    for ion in ion_list:
        
        nion = ion.replace("_"," ")

        out_path_ion = args.out_dir+"/"+str(ion)+"/"

        plot_path = out_path_ion+args.var_arg+"_plots/"

        for i,cut in enumerate(np.arange(iter_vals[0],iter_vals[1],iter_vals[2])):

            try:
                with open(out_path_ion+f"density/{name1}/{name1}_dens_"+args.var_arg+"_"+str(np.round(cut,2))+".pickle", "rb") as dens_dat:
                    uvb1_dens = pickle.load(dens_dat)

                with open(out_path_ion+f"clump_data/{name1}/{name1}_clump_dat_"+args.var_arg+"_"+str(np.round(cut,2))+".pickle", "rb") as salsa_dat:
                    uvb1_clump_dat = pickle.load(salsa_dat)

            except AttributeError:
                print("One or more settings found no absorbers. Skipping")
                pass

            ray_pos = ray_dat.r[("gas","l")].to("kpc")

            try:
                # print(uvb1_clump_dat[nion][name1])
                uvb1_clumps = uvb1_clump_dat[nion][name1][uvb1_clump_dat[nion][name1]["lightray_index"] == f"{ray:0{n}d}"]
            
            except TypeError:
                print("One or more settings found no absorbers. Skipping")
                pass

            plot_df[ion][cut] = np.concatenate((plot_df[ion][cut],uvb1_clumps["col_dens"]), axis = None)


print("Data loaded, generating plots...")

n_iter = (iter_vals[1]-iter_vals[0])//iter_vals[2]
colors = plt.cm.Set1(np.linspace(0,1,int(n_iter)+1))
# for ion in ion_list:
ion = 'H_I'
nion = ion.replace("_"," ")

plot_path = args.out_dir+"/"+str(ion)+"/"
if os.path.exists(plot_path) == False:
    os.mkdir(plot_path)

col_list = []

# determining resolution of figure
res = 100

hist_dat = []
max_val = -np.inf
for i,cut in enumerate(np.arange(iter_vals[0],iter_vals[1],iter_vals[2])):
    temp_hist = np.histogram(plot_df[ion][cut], bins=20)
    hist_dat.append(temp_hist)

    if np.max(temp_hist[0]) > max_val:
        max_val = np.max(temp_hist[0])

for i,cut in enumerate(np.arange(iter_vals[0],iter_vals[1],iter_vals[2])):
    hist = hist_dat[i]
    len(hist[0])
    bheight = np.hstack((hist[0][0], np.array(hist[0])))
    len(bheight)
    sn.lineplot(x=hist[1], y=bheight/max_val, drawstyle='steps-pre', color=colors[i])
    col_list.append(mlines.Line2D([], [], color=colors[i],
                    markersize=10, label=args.var_arg+" = "+str(np.round(cut,2))))

plt.legend(handles=col_list)
plt.grid()
plt.xlabel(r"$\sigma_{FG09}$ ($cm^{-2}$)", fontsize = 15)
plt.title(f"UVB Number Density Comparison "+nion, fontsize = 20)
plt.savefig(plot_path+"UVB_dens_compare_"+args.var_arg+"_test.pdf")
plt.show()

print("Run Complete!")