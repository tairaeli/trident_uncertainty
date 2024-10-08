import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import pickle
import configparser
import os
from sklearn.linear_model import LinearRegression

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

old_gen = [uvb_names[0], uvb_names[2], uvb_names[1]]
new_gen = [uvb_names[1], uvb_names[3], uvb_names[3]]


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

# # creating dictionaries to store our data if they don't already exist
uvb_dist_path = rs_path+"/uvb_dists"
if os.path.exists(uvb_dist_path) == False:
    os.mkdir(uvb_dist_path)

num_ion = len(ion_list)
comp_paths = {}

# setting labels for clump categories
clump_cat_labels = ["match","diff_size","overlap","merge"]
prop_list = ["density","temperature"]
plot_titles = ["n","T"]
prop_unit = ["($cm^{-3}$)","(K)"]
palt = plt.cm.tab10(np.linspace(0,1,4))

# making the super grid
category_list = ["col_dens","dens","t"]
wrats = [1]*len(category_list)
wrats.append(0.5)
fig = plt.figure(figsize=(20,16))
gs = GridSpec(len(ion_list),len(category_list)+1, width_ratios=wrats)
ax_lab_size = 15
title_size = 15
ion_txt_size=30

for i in range(len(old_gen)):
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
        # phys_quant["avg_ray_met"] = {"mean":np.array([]),
        #                               "upper":np.array([]),
        #                               "lower":np.array([])}
        # temporary: adding color array
        # color_arr = np.array([])
        c = 0
        
        color_arr = []
        # old_lone = []
        # new_lone = []
        for ray in comp_dict[old_gen[i]][nion].keys():
            old_lone = []
            new_lone = []
            
            # creating colorbar for clump categories
            match, shorter, longer, overlap, split, merge, lonely_old, lonely_new = clump_categories[nion][ray]
            # print(ray, clump_categories[nion][ray], old_gen[i], new_gen[i], len(comp_dict[old_gen[i]][nion][ray]["col_dens"]))
            
            clumps_containted_1 = 0
            clumps_containted_2 = 0
            # some weird edge case causing me to ad a +2
            # print(comp_dict.keys())
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
                
                # else:
                #     print("Error: Old",ray,j)

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
            
            # if len(reduced_uvb_old["col_dens"])+len(lonely_old["col_dens"]) != 0:
                
            #     frac_lone_old = len(lonely_old["col_dens"])/(len(reduced_uvb_old["col_dens"])+len(lonely_old["col_dens"]))
            #     frac_lone_new = len(lonely_new["col_dens"])/(len(reduced_uvb_new["col_dens"])+len(lonely_new["col_dens"]))
            # else:
            #     frac_lone_old = "INVALID"
            #     frac_lone_new = "INVALID"
            
            assert len(lonely_new["col_dens"]) == len(new_lone), "NEW FAIL "+str(len(lonely_new))+":"+str(len(new_lone))
            assert len(lonely_old["col_dens"]) == len(old_lone), "OLD FAIL "+str(len(lonely_old))+":"+str(len(old_lone))
            assert len(reduced_uvb_old["col_dens"]) == len(reduced_uvb_new["col_dens"]), "Arrays are different sizes. \n old = "+str(len(reduced_uvb_old["col_dens"]))+ "\n new = "+str(len(reduced_uvb_new["col_dens"]))

            # keeping track of total number of lonely clumps
            lonely_new_tot += len(new_lone)
            lonely_old_tot += len(old_lone)

            # keeping track of total clump number
            total_new += len(comp_dict[new_gen[i]][nion][ray]["col_dens"])
            total_old += len(comp_dict[old_gen[i]][nion][ray]["col_dens"])

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
                
                ub = np.zeros_like(reduced_uvb_old["col_dens"])
                lb = np.zeros_like(reduced_uvb_old["col_dens"])
                
                old = reduced_uvb_old[prop]
                new = reduced_uvb_new[prop]
                key = list(phys_quant.keys())[k]

                ub_mask = np.where(old>new)
                lb_mask = np.where(old<new)

                ub[ub_mask] = old[ub_mask]
                ub[lb_mask] = new[lb_mask]

                lb[lb_mask] = old[lb_mask]
                lb[ub_mask] = new[ub_mask]

                old_new_arr = np.array([new,old])
                phys_quant[key]["mean"] = np.concatenate((phys_quant[key]["mean"],
                                                          np.average(old_new_arr, axis = 0,
                                                                     weights = new_old_col_dens)))
                phys_quant[key]["lower"] = np.concatenate((phys_quant[key]["lower"], lb))
                phys_quant[key]["upper"] = np.concatenate((phys_quant[key]["upper"], ub))
        
        color_arr = np.array(color_arr)
        # 1:1 comp
        # setting legend
        new_lone_frac = lonely_new_tot/total_new
        old_lone_frac = lonely_old_tot/total_old
        legend_labs = [Line2D([0], [0], color='w', lw=4, label=f"{old_gen[i]}:{np.round(old_lone_frac,2)}\n {new_gen[i]}:{np.round(new_lone_frac,2)}")]
        
        # line where both UVB quantities match
        match_line = np.linspace(min(np.min(new_gen_dat),np.min(old_gen_dat)),max(np.max(new_gen_dat),np.max(old_gen_dat)))
        
        # creating subplot object
        plot = fig.add_subplot(gs[m,0])
        plot.plot(match_line, match_line, color = "black", linestyle = '--')
        
        # adding ion label
        plot.text(np.min(old_gen_dat)+np.ptp(old_gen_dat)*0.7, np.min(old_gen_dat)+np.ptp(old_gen_dat)*0.2, 
                  s=nion, fontsize=ion_txt_size)
        
        # plotting data
        for k,cat in enumerate(clump_cat_labels):
            cat_mask = np.where(color_arr == cat)
            plot.plot(old_gen_dat[cat_mask], new_gen_dat[cat_mask],
                      linestyle = "None", marker = ".", color = palt[k], label = cat)
            legend_labs.append(Line2D([0], [0], color=palt[k], lw=4, label=cat))
        
        # making a linear fit to the data
        lr = LinearRegression()
        lr.fit(old_gen_dat.reshape(-1,1), new_gen_dat)

        # plotting fit line
        plot.plot(match_line, lr.predict(match_line.reshape(-1,1)), 
                  color = "green", linestyle = ":", linewidth=4)
        legend_labs.append(Line2D([0], [0], color="green", lw=4, 
                                  label=f"Fit line", linestyle = ":"))

        # calculating RMSE
        d = np.abs(old_gen_dat*lr.coef_[0] - new_gen_dat + lr.intercept_)/np.sqrt(lr.coef_[0]**2+1)
        rmse = np.sqrt(np.mean(d))
        legend_labs.append(Line2D([0], [0], color="w", lw=4, label=f"Fit RMSE:{np.round(rmse,3)}"))

        # other plot settings
        plot.set_xlabel(f"{old_gen[i]} $log_{10}$($\sigma$) ($cm^{-2}$)", fontsize=ax_lab_size)
        plot.set_ylabel(f"{new_gen[i]} $log_{10}$($\sigma$) ($cm^{-2}$)", fontsize=ax_lab_size)
        plot.set_title(f"{ion} Absorber $\sigma$ Compare", fontsize=title_size)
        plot.grid()
        plot.legend(handles=legend_labs)

        # comparing other phys quantities
        for j, quant in enumerate(phys_quant.keys()):
            plot = fig.add_subplot(gs[m,1+j])
            for k,cat in enumerate(clump_cat_labels):
                cat_mask = np.where(color_arr == cat)
                plot.plot(np.log10(phys_quant[quant]["mean"][cat_mask]), uvb_dens_diff[cat_mask],
                        linestyle = "None", marker = ".", color = palt[k], label = cat)
                plot.hlines(y=uvb_dens_diff[cat_mask],xmin=np.log10(phys_quant[quant]["lower"][cat_mask]), 
                    xmax=np.log10(phys_quant[quant]["upper"][cat_mask]),
                    color = palt[k])
            
            plot.set_xlabel("$log_{10}$ Mean "+plot_titles[j]+f" {prop_unit[j]}", fontsize=ax_lab_size)
            plot.set_title(f"{ion} Mean "+plot_titles[j]+f" {old_gen[i]} v {new_gen[i]}", fontsize=title_size)
            plot.grid()
            
            # hardcoding some lines in because latex and variable strings do not play well together
            if (i==0) and (j==0):
                plot.set_ylabel(r"$log_{10}$($\frac{\sigma_{FG20}}{\sigma_{FG09}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
            elif (i==1) and (j==0):
                plot.set_ylabel(r"$log_{10}$($\frac{\sigma_{PCW19}}{\sigma_{HM12}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
            elif (i==2) and (j==0):
                plot.set_ylabel(r"$log_{10}$($\frac{\sigma_{PCW19}}{\sigma_{FG20}})$) ($cm^{-2}$)", fontsize=ax_lab_size)

        hist = fig.add_subplot(gs[m,len(prop_list)+1])
        # might use these
        # weights = 10**uvb_dens_diff / sum(10**uvb_dens_diff)

        for k,cat in enumerate(clump_cat_labels):
            cat_mask = np.where(color_arr == cat)
            hist.hist(uvb_dens_diff[cat_mask], orientation="horizontal", 
                      color = palt[k], density=True, bins=20, alpha = 0.7,
                      ec = palt[k], lw=1)
        hist.grid()
    
    # figure settings
    plt.tight_layout()
    plt.savefig(comp_paths[f"{old_gen[i]}_{new_gen[i]}"]+f"/super_plot_{old_gen[i]}_{new_gen[i]}.pdf")
    plt.show()
    plt.clf()