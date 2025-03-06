"""
Generates a violin plot of the distributions of column density 
differences as a function of ion (in order of increasing 
ionization energy)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import pickle
import configparser
import os
import pandas as pd

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

# dict containing redshift values
get_true_rs = {20:'2.0',18:'2.5'}

# gets the true rs needed
true_rs = get_true_rs[int(rs)]

# identifying the path argument as a variable
out_file = sal_args["base_settings"]["output_file"]
path = os.path.expandvars(os.path.expanduser(out_file))
halo_path = path+'/halo'+f'{halo}'
rs_path = halo_path + '/redshift'+f'{true_rs}'

# initializing UVB data
uvb_names = sal_args["uvb_analysis"]["uvb_names"].split(" ")

sal_dat = None

# list of uvb data for comparison
old_gen = [uvb_names[0], uvb_names[2], uvb_names[3]]
new_gen = [uvb_names[1], uvb_names[3], uvb_names[1]]

# contains shortened UVB names
short_uvb_names = {"FG_2009":"FG09", "FG_2020":"FG20",
                   "HM_2012":"HM12", "PCW_2019":"PW19"}

# reformatting ion labels
ion_name_dict = {"H_I":r"H $\mathrm{\i}$",
                 "Si_II":r"Si $\mathrm{\i\i}$",
                 "Si_III":r"Si $\mathrm{\i\i\i}$",
                 "C_III":r"C $\mathrm{\i\i\i}$",
                 "Si_IV":r"Si $\mathrm{\i v}$",
                 "C_IV":r"C $\mathrm{\i v}$",
                 "N_V":r"N $\mathrm{v}$",
                 "O_VI":r"O $\mathrm{v\i}$"}

# creating a comparison dictionary, then assigning UVB data to dictionary
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

# loading in ion data
ion_list = sal_args["galaxy_settings"]["ions"].split(" ")
ion_dat = pd.read_csv("./ion_dat.txt", delimiter = "  ", index_col = "ion")
ion_dat = ion_dat.sort_values(by=["ionization energy (eV)"])

# creating dictionaries to store our data if they don't already exist
uvb_dist_path = rs_path+"/uvb_dists"
if os.path.exists(uvb_dist_path) == False:
    os.mkdir(uvb_dist_path)

num_ion = len(ion_list)
comp_paths = {}

uvb_colors = plt.cm.tab10(np.linspace(0,1,3))

data_list = []

for j, ion in enumerate(ion_list):
    nion = ion.replace("_"," ")

    for i in range(len(old_gen)):
        try:
            with open(rs_path+f"/uvb_clump_labels_{old_gen[i]}_{new_gen[i]}.pickle","rb") as file:
                    clump_categories = pickle.load(file)

            with open(rs_path+f"/uvb_compare_{old_gen[i]}_{new_gen[i]}.pickle","rb") as file:
                comp_dict = pickle.load(file)
        except:
            try:
                with open(rs_path+f"/uvb_compare_{new_gen[i]}_{old_gen[i]}.pickle","rb") as file:
                    comp_dict = pickle.load(file)
                
                with open(rs_path+f"/uvb_compare_{new_gen[i]}_{old_gen[i]}.pickle","rb") as file:
                    comp_dict = pickle.load(file)

            except:
                print(f"Comparison between {new_gen[i]} and {old_gen[i]} Does not exist")
                continue
        
        print(f"Comparing {new_gen[i]} and {old_gen[i]}")

        # setting up data bins to store data in
        uvb_dens_diff = np.array([])
        old_gen_dat = np.array([])
        new_gen_dat = np.array([])

        lonely_new_tot = 0
        lonely_old_tot = 0
        total_new = 0
        total_old = 0

        for ray in comp_dict[old_gen[i]][nion].keys():
            old_lone = []
            new_lone = []
            
            # creating colorbar for clump categories
            match, shorter, longer, overlap, split, merge, lonely_old, lonely_new = clump_categories[nion][ray]
            
            clumps_containted_1 = 0
            clumps_containted_2 = 0

            # some weird edge case causing me to ad a +2
            for j in range(len(comp_dict[old_gen[i]][nion][ray]["col_dens"])+2):
                if (j) in split:
                    clumps_containted_2 += split[j]
                
            for j in range(len(comp_dict[new_gen[i]][nion][ray]["col_dens"])+2): 
                if (j) in merge:
                    clumps_containted_1 += merge[j]

            for j in range(len(comp_dict[old_gen[i]][nion][ray]["col_dens"])+clumps_containted_1):

                if (j) in lonely_old:
                    old_lone.append(j)

            for j in range(len(comp_dict[new_gen[i]][nion][ray]["col_dens"])+clumps_containted_2): 
                if (j) in lonely_new:
                    new_lone.append(j)
                elif (j) in match:
                    continue

            # removing lonely clumps from comparison
            lonely_new, lonely_old, reduced_uvb_new, reduced_uvb_old = lonely_hunter(comp_dict[new_gen[i]][nion][ray],
                                                                                     comp_dict[old_gen[i]][nion][ray])
            # calculating fractions of old and new uvb lonely absorbers          
            assert len(lonely_new["col_dens"]) == len(new_lone), "NEW FAIL "+str(len(lonely_new))+":"+str(len(new_lone))
            assert len(lonely_old["col_dens"]) == len(old_lone), "OLD FAIL "+str(len(lonely_old))+":"+str(len(old_lone))
            assert len(reduced_uvb_old["col_dens"]) == len(reduced_uvb_new["col_dens"]), "Arrays are different sizes. \n old = "+str(len(reduced_uvb_old["col_dens"]))+ "\n new = "+str(len(reduced_uvb_new["col_dens"]))

            # keeping track of total number of lonely clumps
            lonely_new_tot += len(new_lone)
            lonely_old_tot += len(old_lone)

            # keeping track of total clump number
            total_new += len(comp_dict[new_gen[i]][nion][ray]["col_dens"])
            total_old += len(comp_dict[old_gen[i]][nion][ray]["col_dens"])

            # calculating density difference and concatenating it to list
            dens_diff =  reduced_uvb_new["col_dens"] - reduced_uvb_old["col_dens"]
            uvb_dens_diff = np.concatenate((uvb_dens_diff,dens_diff))

            # storing old/new generation
            old_gen_dat = np.concatenate((old_gen_dat, reduced_uvb_old["col_dens"]))
            new_gen_dat = np.concatenate((new_gen_dat, reduced_uvb_new["col_dens"]))


        # setting x-positions for each ion or each UVB comparison


        data_list.append(uvb_dens_diff)

outlier_list = []

for i,dat in enumerate(data_list):
    
    outlier_low, median, outlier_high = np.percentile(data_list[i], [0.5,50,99.5])
    
    outlier_ids = np.where((dat < outlier_low) | (dat > outlier_high))
    # not_outliers = np.where(np.abs(dat[outlier_ids] - median) < 1 )
    # outlier_ids = np.where(np.abs(dat[outlier_ids] - median) > 1 )
    rest_ids = np.where((dat > outlier_low) & (dat < outlier_high))
    # rest_ids = np.concatenate((rest_ids[0],not_outliers[0]))
    dat_temp = dat[outlier_ids]
    data_list[i] = dat[rest_ids] 
    outlier_list.append(dat_temp)

# creating figure
fig = plt.figure(figsize=(20,12))
gs = GridSpec(1,2, width_ratios=[1/len(ion_list),1],
                  wspace=0.15)

print("PLOTTING HI FIGURE")
plot_HI = fig.add_subplot(gs[0])
violin = plot_HI.violinplot(data_list[0:3], [0,1,2],
                        showmeans=False, showmedians=False, showextrema=False)

for i,pc in enumerate(violin['bodies']):
    pc.set_facecolor(uvb_colors[i%3])
    pc.set_edgecolor('black')
    pc.set_alpha(1)

for i in range(3):
    q1, medians, q3 = np.percentile(data_list[i], [25, 50, 75])
    plot_HI.scatter([i], medians, marker='.', color='white', s=30, zorder=3)
    plot_HI.vlines([i], q1, q3, color='k', linestyle='-', lw=5)
    plot_HI.scatter([i]*len(outlier_list[i]),outlier_list[i], color = uvb_colors[i%3])

plot_HI.set_ylabel(r"$log$($N_{1}}$/$N_{2})$", size = 40)
plot_HI.set_yticks([8,6,4,2,0,-2,-4,-6,-8])
plot_HI.set_ylim(-8,8)
plot_HI.set_xticks([1], 
                   labels=[ion_name_dict["H_I"]+"\n"], 
                   size=40)
plot_HI.text(0.25, -0.1, 
                s=str(ion_dat["ionization energy (eV)"][0]), fontsize=20,
                transform=plot_HI.transAxes)

plot_HI.tick_params(axis='y', which='major', labelsize=30)
plot_HI.grid(axis='y')

print("PLOTTING FIGURE")
pos = [4*(k//3)+k%3 for k in range(21)]
# print(len(data_list[3:]), len(pos))
plot = fig.add_subplot(gs[1])
violin = plot.violinplot(data_list[3:], pos,
                        showmeans=False, showmedians=False, 
                        showextrema=False)

for i in range(len(pos)):
    q1, medians, q3 = np.percentile(data_list[3+i], [25, 50, 75])
    plot.scatter(pos[i], medians, marker='.', color='white', s=30, zorder=3)
    plot.vlines(pos[i], q1, q3, color='k', linestyle='-', lw=5)
    plot.scatter([pos[i]]*len(outlier_list[3+i]),outlier_list[3+i], color = uvb_colors[i%3])

for i,pc in enumerate(violin['bodies']):
    pc.set_facecolor(uvb_colors[i%3])
    pc.set_edgecolor('black')
    pc.set_alpha(1)

nion_list = [" "]*len(ion_list)
for i,ion in enumerate(ion_list):
    nion_list[i] = f'{ion_name_dict[ion]} \n'
    if i==0:
        continue
    else:
        plot.text(0.055+(i-1)*0.138, -0.1, 
                s=str(ion_dat["ionization energy (eV)"][i]), fontsize=20,
                transform=plot.transAxes)

legend_labs = []
for i in range(len(old_gen)):
    labels = short_uvb_names[new_gen[i]]+"/"+short_uvb_names[old_gen[i]]
    legend_labs.append(Line2D([0], [0], color=uvb_colors[i], lw=4, label=labels))

plot.set_xticks(np.arange(1, len(pos) + 8,4), 
                labels=nion_list[1:], 
                size=40)
plot.set_yticks([4,3,2,1,0,-1,-2,-3,-4])
plot.set_ylim(-4,4)
plot.tick_params(axis='y', which='major', labelsize=30)
plot.grid(axis='y')
plot.legend(loc='upper right', fontsize=30, handles=legend_labs)

plot.text(0.2,-0.16,
          s="Ionization Energy (eV)", size = 30,
          transform=plot.transAxes)

plt.savefig(rs_path+"/uvb_dists"+"/summary_plot.png",
            dpi=400,bbox_inches='tight')
print("PLOTTING COMPLETE")