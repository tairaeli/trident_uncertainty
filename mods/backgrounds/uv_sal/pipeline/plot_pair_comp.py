"""
Plots pairwise comparisons of column densities along with gas densities
and temperatures
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import pickle
import configparser
import os
from sklearn.linear_model import LinearRegression
import unyt as u
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

ion_list = sal_args["galaxy_settings"]["ions"].split(" ")

# loading in ion data
ion_dat = pd.read_csv("./ion_dat.txt", delimiter = "  ", index_col = "ion")
ion_dat = ion_dat.sort_values(by=["ionization energy (eV)"])

# creating dictionaries to store our data if they don't already exist
uvb_dist_path = rs_path+"/uvb_dists"
if os.path.exists(uvb_dist_path) == False:
    os.mkdir(uvb_dist_path)

num_ion = len(ion_list)
comp_paths = {}

# setting labels for clump categories
clump_cat_labels = ["match","diff_size","overlap","merge"]
prop_list = ["density","temperature"]
plot_titles = ["N","T"]
prop_unit = ["[$cm^{-3}$]","[K]"]
palt = plt.cm.tab10(np.linspace(0,1,4))

# defining the super grid parameters
category_list = ["col_dens","dens","t"]
wrats = [1]*len(category_list)
wrats.append(0.5)
ax_lab_size = 25
title_size = 15
ion_txt_size=40

plot_axis = {"FG_2020:FG_2009":{"col_dens":[(11.9,16.5),(-7.5,7)],
                                "avg_ray_dens":[(-6,1)],
                                "avg_ray_temp":[(2.5,6)]},
             "PCW_2019:HM_2012":{"col_dens":[(11.4,22),(-5,6)],
                                "avg_ray_dens":[(-6,1)],
                                "avg_ray_temp":[(2.5,6)]},
             "FG_2020:PCW_2019":{"col_dens":[(11.4,22),(-6.1,4.5)],
                                "avg_ray_dens":[(-6,1)],
                                "avg_ray_temp":[(2.5,6)]}}

for i in range(len(old_gen)):
    # setting up grid
    fig = plt.figure(figsize=(20,25))
    gs = GridSpec(len(ion_list),len(category_list)+1, width_ratios=wrats,
                  wspace=0.15)
    
    # making the super plot
    for m,ion in enumerate(ion_list):
        nion = ion.replace("_"," ")
        print(f"Creating {nion}")

        # making sure that the right comp file is read in
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
        print("fraction of missing rays:"+str(comp_dict["bad_ray_perc"]))

        # creating dictionaries to store our data if they don't already exist
        comp_paths[f"{old_gen[i]}_{new_gen[i]}"] = rs_path+f"/{old_gen[i]}_{new_gen[i]}_comp"
        if os.path.exists(comp_paths[f"{old_gen[i]}_{new_gen[i]}"]) == False:
            os.mkdir(comp_paths[f"{old_gen[i]}_{new_gen[i]}"])
        
        ion_path = comp_paths[f"{old_gen[i]}_{new_gen[i]}"]+"/"+ion
        if os.path.exists(ion_path) == False:
            os.mkdir(ion_path)

        # setting up data bins to store data in
        uvb_dens_diff = np.array([])
        old_gen_dat = np.array([])
        new_gen_dat = np.array([])

        lonely_new_tot = 0
        lonely_old_tot = 0
        total_new = 0
        total_old = 0

        # this doesn't feel like a good way to store our data
        phys_quant = {}
        phys_quant["avg_ray_dens"] = {"mean":np.array([]),
                                      "upper":np.array([]),
                                      "lower":np.array([])}
        phys_quant["avg_ray_temp"] = {"mean":np.array([]),
                                      "upper":np.array([]),
                                      "lower":np.array([])}

        color_arr = []

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
                if (j) in match:
                    color_arr.append("match")

                elif (j) in split:
                    color_arr.append("merge")

                elif (j) in lonely_old:
                    old_lone.append(j)
                
                elif (j) in longer:
                    color_arr.append("diff_size")
                
                elif (j) in shorter:
                    color_arr.append("diff_size")

                elif (j) in overlap:
                    color_arr.append("overlap")

            for j in range(len(comp_dict[new_gen[i]][nion][ray]["col_dens"])+clumps_containted_2): 
                if (j) in lonely_new:
                    new_lone.append(j)
                elif (j) in merge:
                    color_arr.append("merge")
                elif (j) in match:
                    continue
                # else:
                    # print("Error: New", ray,j)

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

            assert len(uvb_dens_diff) == len(color_arr), "Arrays diff sizes at ray "+str(ray)+". uvb_dens_diff:"+str(len(uvb_dens_diff))+" color_arr:"+str(len(color_arr))+" "+old_gen[i]+" "+new_gen[i]

            # adding upper and lower bounds
            new_old_col_dens = np.array([reduced_uvb_old["col_dens"],
                                         reduced_uvb_new["col_dens"]])

            for k, prop in enumerate(prop_list):
                # finding upper and lower bounds of data
                ub = np.zeros_like(reduced_uvb_old["col_dens"])
                lb = np.zeros_like(reduced_uvb_old["col_dens"])
                
                old = reduced_uvb_old[prop]
                new = reduced_uvb_new[prop]
                key = list(phys_quant.keys())[k]

                # creating data masks for 
                ub_mask = np.where(old>new)
                lb_mask = np.where(old<new)

                # all indices where the phys quantity is larger goes in ub
                ub[ub_mask] = old[ub_mask]
                ub[lb_mask] = new[lb_mask]
                
                # all indices where the phys quantity is smaller goes in lb
                lb[lb_mask] = old[lb_mask]
                lb[ub_mask] = new[ub_mask]

                old_new_arr = np.array([new,old])

                # assigning data to dictionary
                phys_quant[key]["mean"] = np.concatenate((phys_quant[key]["mean"],
                                                          np.average(old_new_arr, axis = 0,
                                                                     weights = new_old_col_dens)))
                phys_quant[key]["lower"] = np.concatenate((phys_quant[key]["lower"], lb))
                phys_quant[key]["upper"] = np.concatenate((phys_quant[key]["upper"], ub))
        
        color_arr = np.array(color_arr)
        ion_mass_g = ion_dat["mass (amu)"][ion]*1.6605E-24

        new_lone_frac = lonely_new_tot/total_new
        old_lone_frac = lonely_old_tot/total_old
        
        print("fraction of new lonely absorbers:"+str(new_lone_frac))
        print("fraction of old lonely absorbers:"+str(old_lone_frac))

        # 1:1 comp
        # setting legend
        legend_labs = []
        
        # line where both UVB quantities match
        match_line = np.linspace(min(np.min(new_gen_dat),np.min(old_gen_dat)),max(np.max(new_gen_dat),np.max(old_gen_dat)))
        
        # creating subplot object
        plot = fig.add_subplot(gs[m,0])
        plot.axhline(0, color = "black", linestyle = '--')
        
        plot.set_xlim(plot_axis[f"{new_gen[i]}:{old_gen[i]}"]["col_dens"][0])
        plot.set_ylim(plot_axis[f"{new_gen[i]}:{old_gen[i]}"]["col_dens"][1])

        # adding ion label
        plot.text(0.68, 0.15, 
                  s=nion, fontsize=ion_txt_size,
                  transform=plot.transAxes)
        
        # plotting data
        for k,cat in enumerate(clump_cat_labels):
            cat_mask = np.where(color_arr == cat)
            plot.plot(old_gen_dat[cat_mask], uvb_dens_diff[cat_mask],
                      linestyle = "None", marker = ".", color = palt[k], label = cat)
            legend_labs.append(Line2D([0], [0], color=palt[k], lw=4, label=cat))
        
        # making a linear fit to the data
        lr = LinearRegression()
        lr.fit(old_gen_dat.reshape(-1,1), uvb_dens_diff)

        # plotting fit line
        plot.plot(match_line, lr.predict(match_line.reshape(-1,1)), 
                  color = "green", linestyle = ":", linewidth=4)

        # calculating RMSE
        d = np.abs(old_gen_dat*lr.coef_[0] - new_gen_dat + lr.intercept_)/np.sqrt(lr.coef_[0]**2+1)
        rmse = np.sqrt(np.mean(d))

        # other plot settings
        if m == len(ion_list)-1:
            if i==0:
                 plot.set_xlabel(r"$log(\sigma_{FG09}$) [$cm^{-2}$]", fontsize=ax_lab_size)
            elif i==1:
                plot.set_xlabel(r"$log(\sigma_{HM12}$) [$cm^{-2}$]", fontsize=ax_lab_size)
            elif i==2:
                plot.set_xlabel(r"$log(\sigma_{PW19}$) [$cm^{-2}$]", fontsize=ax_lab_size)

            # plot.set_xlabel(rf"$log(\sigma_{short_uvb_names[old_gen[i]]}$) [$cm^{-2}$]", fontsize=ax_lab_size)
        
        plot.grid()
        if m == 0:
            plot.legend(handles=legend_labs, loc='upper right')

        # converting gas density to number density
        phys_quant["avg_ray_dens"]["mean"] = phys_quant["avg_ray_dens"]["mean"]/ion_mass_g
        phys_quant["avg_ray_dens"]["lower"] = phys_quant["avg_ray_dens"]["lower"]/ion_mass_g
        phys_quant["avg_ray_dens"]["upper"] = phys_quant["avg_ray_dens"]["upper"]/ion_mass_g

        # creating comparison plots
        for j, quant in enumerate(phys_quant.keys()):
            
            # creating a linear fit to calculate variance of data
            lr = LinearRegression()
            phys_quant_mean = phys_quant[quant]["mean"][~np.isnan(phys_quant[quant]["mean"])]
            uvb_dens_diff_rmna = uvb_dens_diff[~np.isnan(phys_quant[quant]["mean"])]
            
            lr.fit(phys_quant_mean.reshape(-1,1), uvb_dens_diff_rmna)

            match_line = np.linspace(np.min(phys_quant_mean),
                                     np.max(phys_quant_mean),
                                     100)
            rss = np.sum((uvb_dens_diff - lr.intercept_ - lr.coef_*phys_quant[quant]["mean"])**2)
            var = rss/(len(uvb_dens_diff)-2)
            print("Variance of",quant+":", var)

            # plotting data
            plot = fig.add_subplot(gs[m,1+j])
            plot.set_xlim(plot_axis[f"{new_gen[i]}:{old_gen[i]}"][quant][0])
            plot.set_ylim(plot_axis[f"{new_gen[i]}:{old_gen[i]}"]["col_dens"][1])

            for k,cat in enumerate(clump_cat_labels):
                cat_mask = np.where(color_arr == cat)
                plot.plot(np.log10(phys_quant[quant]["mean"][cat_mask]), uvb_dens_diff[cat_mask],
                        linestyle = "None", marker = ".", color = palt[k], label = cat)
                plot.hlines(y=uvb_dens_diff[cat_mask],xmin=np.log10(phys_quant[quant]["lower"][cat_mask]), 
                    xmax=np.log10(phys_quant[quant]["upper"][cat_mask]),
                    color = palt[k])
            
            ylab_xpos = 0.08
            if m == (len(ion_list)-1):
                # plot.set_xlabel(names_one_to_one[old_gen[i]], fontsize=35)
                # hardcoding some lines in because latex and variable strings do not play well together
                if i==0:
                    # plot.set_ylabel(r"$log_{10}$($\sigma_{FG20}}$/$\sigma_{FG09}$) [$cm^{-2}$]", fontsize=ax_lab_size)
                    fig.text(ylab_xpos, 0.5, r"$log$($\sigma_{FG20}}$/$\sigma_{FG09}$) [$cm^{-2}$]", 
                             ha='center', va='center', rotation='vertical', fontsize=35)
                elif i==1:
                    # plot.set_ylabel(r"$log_{10}$($\sigma_{PW19}$/$\sigma_{HM12}$) [$cm^{-2}$]", fontsize=ax_lab_size)
                    fig.text(ylab_xpos, 0.5, r"$log$($\sigma_{PW19}$/$\sigma_{HM12}$) [$cm^{-2}$]", 
                             ha='center', va='center', rotation='vertical', fontsize=35)
                elif i==2:
                    # plot.set_ylabel(r"$log_{10}$($\sigma_{FG20}$/$\sigma_{PW19}$) [$cm^{-2}$]", fontsize=ax_lab_size)
                    fig.text(ylab_xpos, 0.5, r"$log$($\sigma_{FG20}$/$\sigma_{PW19}$) [$cm^{-2}$]", 
                             ha='center', va='center', rotation='vertical', fontsize=35)
                # fig.text(0.083+0.285, 0.5, names_one_to_one[new_gen[i]], ha='center', va='center', 
                # rotation='vertical', fontsize=35)

                if j==0:
                    plot.set_xlabel("$log$"+r"$(\overline{n})$"+f" {prop_unit[j]}", fontsize=ax_lab_size)
                else:
                    plot.set_xlabel("$log$"+r"$(\overline{T})$"+f" {prop_unit[j]}", fontsize=ax_lab_size)
            
            # plot.set_yticks([])
            plot.grid(True)

        # creating a histogram of column density differences 
        hist = fig.add_subplot(gs[m,len(prop_list)+1])
        for k,cat in enumerate(clump_cat_labels):
            cat_mask = np.where(color_arr == cat)
            hist.hist(uvb_dens_diff[cat_mask], orientation="horizontal", 
                      color = palt[k], density=True, bins=20, alpha = 0.7,
                      ec = palt[k], lw=1)
        # hist.set_yticks([])
        hist.grid(True)
    
    # figure settings
    plt.tight_layout()
    plt.savefig(comp_paths[f"{old_gen[i]}_{new_gen[i]}"]+f"/super_plot_{old_gen[i]}_{new_gen[i]}.png",
                dpi=400,bbox_inches='tight')
    plt.show()
    plt.clf()

# creating 1:1 column density comparison figures
wrats = [1,1,1]
fig = plt.figure(figsize=(20,25))
gs = GridSpec(len(ion_list),len(old_gen), width_ratios=wrats,
              wspace=0.4)

map_i = {0:0, 1:2, 2:4}

# contains shortened UVB names
names_one_to_one = {"FG_2009":r"$\sigma_{FG09}$", "FG_2020":r"$\sigma_{FG20}$",
                   "HM_2012":r"$\sigma_{HM12}$", "PCW_2019":r"$\sigma_{PW19}$"}

# puts each 1:1 comparison in above figure
for i in range(len(old_gen)):
    for m,ion in enumerate(ion_list):
        nion = ion.replace("_"," ")
        print(f"Creating {nion}")

        # making sure that the right comp file is read in
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
        print("fraction of missing rays:"+str(comp_dict["bad_ray_perc"]))
        
        # setting up data bins to store data in
        uvb_dens_diff = np.array([])
        old_gen_dat = np.array([])
        new_gen_dat = np.array([])
        color_arr = []


        for ray in comp_dict[old_gen[i]][nion].keys():
            old_lone = []
            new_lone = []
            
            # creating colorbar for clump categories
            match, shorter, longer, overlap, split, merge, lonely_old, lonely_new = clump_categories[nion][ray]
            
            clumps_containted_1 = 0
            clumps_containted_2 = 0

            # adding data to colorbar categories
            # some weird edge case causing me to ad a +2
            for j in range(len(comp_dict[old_gen[i]][nion][ray]["col_dens"])+2):
                if (j) in split:
                    clumps_containted_2 += split[j]
                
            for j in range(len(comp_dict[new_gen[i]][nion][ray]["col_dens"])+2): 
                if (j) in merge:
                    clumps_containted_1 += merge[j]

            for j in range(len(comp_dict[old_gen[i]][nion][ray]["col_dens"])+clumps_containted_1):
                if (j) in match:
                    color_arr.append("match")

                elif (j) in split:
                    color_arr.append("merge")

                elif (j) in lonely_old:
                    old_lone.append(j)
                
                elif (j) in longer:
                    color_arr.append("diff_size")
                
                elif (j) in shorter:
                    color_arr.append("diff_size")

                elif (j) in overlap:
                    color_arr.append("overlap")

            for j in range(len(comp_dict[new_gen[i]][nion][ray]["col_dens"])+clumps_containted_2): 
                if (j) in lonely_new:
                    new_lone.append(j)
                elif (j) in merge:
                    color_arr.append("merge")
                elif (j) in match:
                    continue

            # removing lonely clumps from comparison
            lonely_new, lonely_old, reduced_uvb_new, reduced_uvb_old = lonely_hunter(comp_dict[new_gen[i]][nion][ray],
                                                                                     comp_dict[old_gen[i]][nion][ray])
            
            # checking the data 
            assert len(lonely_new["col_dens"]) == len(new_lone), "NEW FAIL "+str(len(lonely_new))+":"+str(len(new_lone))
            assert len(lonely_old["col_dens"]) == len(old_lone), "OLD FAIL "+str(len(lonely_old))+":"+str(len(old_lone))
            assert len(reduced_uvb_old["col_dens"]) == len(reduced_uvb_new["col_dens"]), "Arrays are different sizes. \n old = "+str(len(reduced_uvb_old["col_dens"]))+ "\n new = "+str(len(reduced_uvb_new["col_dens"]))

            dens_diff =  reduced_uvb_new["col_dens"] - reduced_uvb_old["col_dens"]

            uvb_dens_diff = np.concatenate((uvb_dens_diff,dens_diff))

            # storing old/new generation
            old_gen_dat = np.concatenate((old_gen_dat, reduced_uvb_old["col_dens"]))
            new_gen_dat = np.concatenate((new_gen_dat, reduced_uvb_new["col_dens"]))

            assert len(uvb_dens_diff) == len(color_arr), "Arrays diff sizes at ray "+str(ray)+". uvb_dens_diff:"+str(len(uvb_dens_diff))+" color_arr:"+str(len(color_arr))+" "+old_gen[i]+" "+new_gen[i]

        color_arr = np.array(color_arr)

        # creating data for a line of equivalency between the data
        # match_line = np.linspace(min(np.min(new_gen_dat),np.min(old_gen_dat)),max(np.max(new_gen_dat),np.max(old_gen_dat)))
        match_line = np.linspace(11,20,100)
        
        # creating subplot object
        # plot_i = map_i[i]
        splot = fig.add_subplot(gs[m,i])
        splot.set_xlim(11,20)
        splot.set_ylim(11,20)
        splot.plot(match_line, match_line, color = "black", linestyle = '--')
        
        # adding ion label
        if i == 0:
            # splot.text(np.min(old_gen_dat)+np.ptp(old_gen_dat)*0.7, np.min(old_gen_dat)+np.ptp(old_gen_dat)*0.2, 
            #         s=nion, fontsize=ion_txt_size)
            splot.text(17, 13, 
                    s=nion, fontsize=ion_txt_size)
        
        legend_labs = []

        # plotting data
        for k,cat in enumerate(clump_cat_labels):
            cat_mask = np.where(color_arr == cat)
            splot.plot(old_gen_dat[cat_mask], new_gen_dat[cat_mask],
                      linestyle = "None", marker = ".", color = palt[k], label = cat)
            legend_labs.append(Line2D([0], [0], color=palt[k], lw=4, label=cat))
        
        # ax[m,i].set_ylabel(new_gen[i])
        splot.grid()
        if m == len(ion_list)-1:
            splot.set_xlabel(names_one_to_one[old_gen[i]], fontsize=35)
            fig.text(0.083+0.285*i, 0.5, names_one_to_one[new_gen[i]], ha='center', va='center', 
                    rotation='vertical', fontsize=35)
        
        if (m==0) and (i==0):
            splot.legend(loc='upper left', fontsize=13)

# plotting figure
# fig.supylabel("FG20",fontsize = 25)
plt.tight_layout()
plt.savefig(rs_path+"/uvb_dists"+f"/super_one_to_one.png",
                dpi=400,bbox_inches='tight')
plt.show()
plt.clf()