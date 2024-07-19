
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
    uvb1_omask = np.where(np.array(uvb1["col_dens"])!=0)

    # preparing data containers
    red_uvb1 = {}
    lone1 = {}
    red_uvb2 = {}
    lone2 = {}

    # removing lonely clumps from each key in dict
    for key in uvb1.keys():
        red_uvb1[key] = np.array(uvb1[key])
        red_uvb1[key] = red_uvb1[key][uvb1_omask]
        lone1[key] = np.array(uvb1[key])
        lone1[key] = lone1[key][uvb1_mask]

        red_uvb2[key] = np.array(uvb2[key])
        red_uvb2[key] = red_uvb2[key][uvb1_omask]

    assert len(uvb1["col_dens"]) >= len(red_uvb1["col_dens"]), "first masking failed"

    # making second mask with previously masked data
    uvb2_mask = np.where(np.array(red_uvb2["col_dens"])==0)
    uvb2_omask = np.where(np.array(red_uvb2["col_dens"])!=0)

    # removing lonely clumps from each key in dict
    for key in uvb1.keys():
        red_uvb1[key] = red_uvb1[key][uvb2_omask]

        lone2[key] = np.array(red_uvb2[key])
        lone2[key] = lone2[key][uvb2_mask]

        red_uvb2[key] = red_uvb2[key][uvb2_omask]

    assert len(red_uvb1["col_dens"]) >= len(red_uvb1["col_dens"]), "second masking failed"

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

old_gen = [uvb_names[0], uvb_names[2]]
new_gen = [uvb_names[1], uvb_names[3]]

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

# creating dictionaries to store our data if they don't already exist
uvb_dist_path = rs_path+"/uvb_dists"
if os.path.exists(uvb_dist_path) == False:
    os.mkdir(uvb_dist_path)

# creating column density plots for each UVB

for ion in sal_dat.keys():

    fig, ax = plt.subplots(2,2, figsize=[8,8])
    i = 0
    j = 0
    for name in uvb_names:
        sns.kdeplot(sal_dat[ion][name]["col_dens"], ax = ax[i,j])
        ax[i,j].set_xlabel(r"log10(Column Density) $cm^{-2}$")
        ax[i,j].set_ylabel(r"Density")
        ax[i,j].set_title(name)

        if (i>0):
            j+=1
            i=0
        else:
            i+=1
    
    plt.tight_layout()
    plt.savefig(uvb_dist_path+f"/dist_plot_{ion}.pdf")

plt.clf()

print("Sorting Data")

comp_paths = {}

# making col dens diff histograms for each ion
for ion in ion_list:

    nion = ion.replace("_"," ")
    print(f"Creating {nion}")

    niter = len(uvb_names)
    colors = plt.cm.Set1(np.linspace(0,1,int(niter)+1))
    k = 0
    # creating pairwise comp plots for each UVB
    for i in range(len(uvb_names)-1):
        for j in range(i+1,len(uvb_names)):
            try:                
                with open(rs_path+f"/uvb_compare_{uvb_names[i]}_{uvb_names[j]}.pickle","rb") as file:
                    comp_dict = pickle.load(file)
            except:
                print(f"Comparison between {uvb_names[i]} and {uvb_names[j]} Does not exist")
                continue

            print(f"Comparing {uvb_names[i]} and {uvb_names[j]}")

            uvb_dens_diff = np.array([])

            # iterating over each ray
            lonely1_tot = 0
            lonely2_tot = 0

            for ray in range(nrays):
                # removing lonely clumps from comparison
                lonely1, lonely2, reduced_uvb1, reduced_uvb2 = lonely_hunter(comp_dict[uvb_names[i]][nion][ray],
                                                                             comp_dict[uvb_names[j]][nion][ray])
                lonely1_tot += len(lonely1)
                lonely2_tot += len(lonely2)

                dens_diff =  reduced_uvb1["col_dens"] - reduced_uvb2["col_dens"]

                uvb_dens_diff = np.concatenate((uvb_dens_diff,dens_diff))
            
            # creating kde
            sns.kdeplot(uvb_dens_diff, 
                        label=f"{uvb_names[i]}/{uvb_names[j]}:{lonely1_tot}:{lonely2_tot}",
                        warn_singular=False)
    
    # creating dictionaries to store our data if they don't already exist
    # ion_path = uvb_dist_path+"/"+ion
    # if os.path.exists(ion_path) == False:
    #     os.mkdir(ion_path)

    # Line2D()
    plt.xlabel(r"log10(Column Density Difference) ($cm^{-2}$)", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title(f"{ion} KDE of Column Density Difference for {nrays} Rays", fontsize=15)
    plt.legend()
    plt.grid()
    plt.savefig(uvb_dist_path+f"/{ion}_dens_diff_hist.pdf")
    plt.clf()

    # fig 4: scatter of old uvbs vs diff btw them and the new iterations
    # fig 5-7: density diff vs physical quantities
    
    for i in range(len(old_gen)):
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

            except:
                print(f"Comparison between {new_gen[i]} and {old_gen[i]} Does not exist")
                continue

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

        # this doesn't feel like a good way to store our data
        phys_quant = {}
        phys_quant["avg_ray_dens"] = {"mean":np.array([]),
                                      "upper":np.array([]),
                                      "lower":np.array([])}
        phys_quant["avg_ray_temp"] = {"mean":np.array([]),
                                      "upper":np.array([]),
                                      "lower":np.array([])}
        phys_quant["avg_ray_met"] = {"mean":np.array([]),
                                      "upper":np.array([]),
                                      "lower":np.array([])}
        # temporary: adding color array
        # color_arr = np.array([])
        c = 0
        
        color_arr = []
        # old_lone = []
        # new_lone = []
        for ray in range(nrays):
            old_lone = []
            new_lone = []

            # creating colorbar for clump categories
            match, shorter, longer, split_new, split_old, lonely_old, lonely_new, overlap = clump_categories[nion][ray]
            
            # color_map = plt.cm.Set1(np.arange(0,7))
            clumps_containted = 0
            for j in range(len(comp_dict[old_gen[i]][nion][ray]["col_dens"])):
                if (j+clumps_containted) in match:
                    color_arr.append(0.0)
                elif (j+clumps_containted) in split_old:
                    color_arr.append(0.5)
                    clumps_containted += split_old[j]
                elif (j+clumps_containted) in lonely_old:
                    old_lone.append(j)
            
            clumps_containted = 0
            for j in range(len(comp_dict[new_gen[i]][nion][ray]["col_dens"])): 
                if (j+clumps_containted) in split_new:
                    color_arr.append(1)
                    clumps_containted += split_new[j]
                elif (j+clumps_containted) in match:
                    continue
                elif (j+clumps_containted) in longer:
                    color_arr.append(0.7)
                elif (j+clumps_containted) in shorter:
                    color_arr.append(0.2)
                elif (j+clumps_containted) in overlap:
                    color_arr.append(0.9)
                elif (j+clumps_containted) in lonely_new:
                    new_lone.append(j)

            # removing lonely clumps from comparison
            # lonely_new, lonely_old, reduced_uvb_new, reduced_uvb_old = lonely_hunter(comp_dict[new_gen[i]][nion][ray],
            #                                                                          comp_dict[old_gen[i]][nion][ray])
            reduced_uvb_new = {}
            reduced_uvb_old = {}
            for key in comp_dict[old_gen[i]][nion][ray].keys():
                reduced_uvb_old[key] = np.delete(np.array(comp_dict[old_gen[i]][nion][ray][key]), old_lone)
                reduced_uvb_new[key] = np.delete(np.array(comp_dict[new_gen[i]][nion][ray][key]), new_lone)

            # print(new_lone, old_lone)
            assert len(reduced_uvb_old["col_dens"]) == len(reduced_uvb_new["col_dens"]), f"Arrays are different sizes. \n old = {len(reduced_uvb_old)} \n new = {len(reduced_uvb_new)}"

            # keeping track of total number of lonely clumps
            lonely_new_tot += len(new_lone)
            lonely_old_tot += len(old_lone)

            dens_diff =  reduced_uvb_new["col_dens"] - reduced_uvb_old["col_dens"]

            uvb_dens_diff = np.concatenate((uvb_dens_diff,dens_diff))

            # storing old/new generation
            old_gen_dat = np.concatenate((old_gen_dat, reduced_uvb_old["col_dens"]))
            new_gen_dat = np.concatenate((new_gen_dat, reduced_uvb_new["col_dens"]))

            # temporary: adding color to each ray
            # ray_color = np.zeros_like(dens_diff) + c
            # color_arr = np.concatenate((color_arr,ray_color))
            # c=ray/nrays

            # adding upper and lower bounds
            for trip in [(reduced_uvb_old["density"],reduced_uvb_new["density"],"avg_ray_dens"),
                         (reduced_uvb_old["temperature"],reduced_uvb_new["temperature"],"avg_ray_temp"),
                         (reduced_uvb_old["metallicity"],reduced_uvb_new["metallicity"],"avg_ray_met")]:
                
                ub = np.zeros_like(reduced_uvb_old["col_dens"])
                lb = np.zeros_like(reduced_uvb_old["col_dens"])
                
                old = trip[0]
                new = trip[1]
                key = trip[2]

                ub_mask = np.where(old>new)
                lb_mask = np.where(old<new)

                ub[ub_mask] = old[ub_mask]
                ub[lb_mask] = new[lb_mask]

                lb[lb_mask] = old[lb_mask]
                lb[ub_mask] = new[ub_mask]

                phys_quant[key]["mean"] = np.concatenate((phys_quant[key]["mean"], (new+old)/2))
                phys_quant[key]["lower"] = np.concatenate((phys_quant[key]["lower"], lb))
                phys_quant[key]["upper"] = np.concatenate((phys_quant[key]["upper"], ub))

        # direct column density comparison
        sns.scatterplot(x=old_gen_dat, y=new_gen_dat,
                        label=f"{old_gen[i]}:{lonely1_tot}\n {new_gen[i]}:{lonely2_tot}")
        match_line = np.linspace(min(np.min(new_gen_dat),np.min(old_gen_dat)),max(np.max(new_gen_dat),np.max(old_gen_dat)))
        sns.lineplot(x=match_line, y=match_line, color = "black", linestyle = '--')
        plt.xlabel(fr"{old_gen[i]} log10(Column Density) ($cm^{-2}$)", fontsize=12)
        plt.ylabel(fr"{new_gen[i]} log10(Column Density) ($cm^{-2}$)", fontsize=12)
        plt.title(f"{ion} Clump By Clump Column Density Compare: {nrays} Rays", fontsize=15)
        plt.grid()
        plt.savefig(ion_path+f"/compare_plot_{old_gen[i]}_{new_gen[i]}_{ion}.pdf")
        plt.clf()

        # old col dens vs col dens diff
        palt = sns.color_palette("viridis", as_cmap=True)
        sns.scatterplot(x=old_gen_dat, y=uvb_dens_diff,
                        label=f"{old_gen[i]}:{lonely1_tot}\n {new_gen[i]}:{lonely2_tot}")
        plt.xlabel(fr"{old_gen[i]} log10(Column Density) ($cm^{-2}$)", fontsize=12)
        plt.ylabel("log10(Column Density Difference)", fontsize=12)
        plt.title(f"{ion} Clump By Clump Density Differences of {old_gen[i]} for {nrays} Rays", fontsize=15)
        plt.grid()
        plt.savefig(ion_path+f"/dens_diff_{old_gen[i]}_{ion}.pdf")
        plt.clf()

        # gas density vs col dens diff
        sns.scatterplot(y=uvb_dens_diff, x=np.log10(phys_quant["avg_ray_dens"]["mean"]), 
                        label=f"{old_gen[i]}:{lonely1_tot}\n {new_gen[i]}:{lonely2_tot}",
                        c = color_arr)
        plt.hlines(y=uvb_dens_diff,xmin=np.log10(phys_quant["avg_ray_dens"]["lower"]), 
                   xmax=np.log10(phys_quant["avg_ray_dens"]["upper"]))
        plt.ylabel("log Column Density Difference ($cm^{-2}$)", fontsize=12)
        plt.xlabel("log Mean Gas Density ($cm^{-3}$)", fontsize=12)
        plt.title(f"{ion} Average Mean Gas Densities of {old_gen[i]} and {new_gen[i]} for {nrays} Rays", fontsize=15)
        plt.grid()
        plt.savefig(ion_path+f"/gas_dens_{old_gen[i]}_{new_gen[i]}_{ion}.pdf")
        plt.clf()

        # temp vs col dens diff
        sns.scatterplot(y=uvb_dens_diff, x=np.log10(phys_quant["avg_ray_temp"]["mean"]), 
                        label=f"{old_gen[i]}:{lonely1_tot}\n {new_gen[i]}:{lonely2_tot}",
                        hue_norm = color_arr, palette=palt)
        plt.hlines(y=uvb_dens_diff,xmin=np.log10(phys_quant["avg_ray_temp"]["lower"]), 
                   xmax=np.log10(phys_quant["avg_ray_temp"]["upper"]))
        plt.ylabel("log Column Density Difference ($cm^{-2}$)", fontsize=12)
        plt.xlabel("Mean Temperature (K)", fontsize=12)
        plt.title(f"{ion} Average Temperature of {old_gen[i]} and {new_gen[i]} for {nrays} Rays", fontsize=15)
        plt.grid()
        plt.savefig(ion_path+f"/temperature_{old_gen[i]}_{new_gen[i]}_{ion}.pdf")
        plt.clf()

        # metallically vs col dens diff
        sns.scatterplot(y=uvb_dens_diff, x=np.log10(phys_quant["avg_ray_met"]["mean"]), 
                        label=f"{old_gen[i]}:{lonely1_tot}\n {new_gen[i]}:{lonely2_tot}",
                        hue_norm = color_arr, palette=palt)
        plt.hlines(y=uvb_dens_diff,xmin=np.log10(phys_quant["avg_ray_met"]["lower"]), 
                   xmax=np.log10(phys_quant["avg_ray_met"]["upper"]))
        plt.ylabel("log10(Column Density Difference) ($cm^{-2}$)", fontsize=12)
        plt.xlabel("Mean Metallicity (K)", fontsize=12)
        plt.title(f"{ion} Average Metallicity of {old_gen[i]} and {new_gen[i]} for {nrays} Rays", fontsize=15)
        plt.grid()
        plt.savefig(ion_path+f"/metallicity_{old_gen[i]}_{new_gen[i]}_{ion}.pdf")
        plt.clf()

# making comparison plots for each ray
colors = plt.cm.cool([0.3,0.4,0.9,1])

for ion in ion_list:
    nion = ion.replace("_"," ")
    for ray in range(nrays):
        # loading in ray data
        n = len(str(nrays))
        ray_dat = yt.load(rs_path+f"/rays/ray{ray:0{n}d}.h5")
        ray_pos = ray_dat.r[("gas","l")].to("kpc")

        for i in range(len(old_gen)):

            plot_path = rs_path+f"/{old_gen[i]}_{new_gen[i]}_comp/{ion}/ray_dat"
            if os.path.exists(plot_path) == False:
                os.mkdir(plot_path)

            old = sal_dat[nion][old_gen[i]][sal_dat[nion][old_gen[i]]["lightray_index"] == f"{ray:0{n}d}"].reset_index()
            new = sal_dat[nion][new_gen[i]][sal_dat[nion][new_gen[i]]["lightray_index"] == f"{ray:0{n}d}"].reset_index()

            with open(rs_path +f'/{old_gen[i]}/ray_dat/{ion}/ray_{ray:0{n}d}_dens.pickle', "rb") as dens_dat:
                    old_num_dens = pickle.load(dens_dat)
            
            with open(rs_path +f'/{new_gen[i]}/ray_dat/{ion}/ray_{ray:0{n}d}_dens.pickle', "rb") as dens_dat:
                    new_num_dens = pickle.load(dens_dat)
            
            plt.figure(figsize = [15,8], dpi = 500, facecolor = "white")
            plt.semilogy(ray_pos, old_num_dens,label = old_gen[i], color = colors[3])
            plt.semilogy(ray_pos, new_num_dens,label = new_gen[i], color = colors[1])

            for j in range(len(old["interval_start"])):
                
                lb = old["interval_start"][j]
                hb = old["interval_end"][j]
                rng = [lb,hb]
                
                yb = old_num_dens[slice(*rng)]
                xb = ray_pos[lb:hb]
                
                plt.fill_between(xb,yb, color = colors[2], alpha = 0.3)
            
            for j in range(len(new["interval_start"])):
                lb = new["interval_start"][j]
                hb = new["interval_end"][j]
                rng = [lb,hb]
                
                yb = new_num_dens[slice(*rng)]
                xb = ray_pos[lb:hb]
                
                plt.fill_between(xb,yb, color = colors[0], alpha = 0.3)
            
            plt.grid()
            plt.legend(fontsize = 20)
            # plt.ylim(10**(-15),10**(-6.5))
            # plt.xlim(500,900)
            plt.xlabel("Position Along Ray (kpc)", fontsize = 25)
            plt.ylabel(r"Density ($cm^{-3}$)", fontsize = 25)
            # plt.ylim(1e-12, 2e-6)
            plt.title(f"UVB Number Density Comparison "+ion, fontsize = 30)
            plt.savefig(plot_path+f"/UVB_dens_compare_ray_{ray:0{n}d}.pdf")
            plt.clf()

# fig 8: making comparisons between species
# C_lis = ["C_II","C_III","C_IV"]
# O_lis = ["O_III","O_VI"]
# N_lis = ["N_II","N_IV","N_V"]
# Si_lis = ["Si_II", "Si_III", "Si_IV"]
# S_lis = ["S_II", "S_III", "S_IV"]
# Mg_lis = ["Mg_II","Mg_X"]
# Al_lis = ["Al_III"]
# Ne_lis = ["Ne_VII","Ne_VIII"]
# H = ["H_I"]

# ion_list = [C_lis,O_lis,N_lis,Si_lis,
#             S_lis,Mg_lis,Al_lis,Ne_lis,H]

# spec_list = ["C","O","N","Si","S","Mg","Al","Ne","H"]

# for i in range(len(old_gen)):
#     try:
#         with open(rs_path+f"/uvb_compare_{old_gen[i]}_{new_gen[i]}.pickle","rb") as file:
#             comp_dict = pickle.load(file)
#     except:
#         try:
#             with open(rs_path+f"/uvb_compare_{new_gen[i]}_{old_gen[i]}.pickle","rb") as file:
#                 comp_dict = pickle.load(file)

#         except:
#             print(f"Comparison between {new_gen[i]} and {old_gen[i]} Does not exist")
#             continue
#     # iterating through each species group in the ion list
#     for k, species in enumerate(ion_list):

#         spec_dict = {}

#         for ion in species:

#             spec_dict[ion] = {"dens_diff":None, "dens_old":None}

#             nion = ion.replace("_"," ")

#             uvb_dens_diff = np.array([])
#             old_gen_dat = np.array([])

#             lonely_old_tot = 0
#             lonely_new_tot = 0
#             for ray in range(nrays):
#                 lonely_new, lonely_old, reduced_uvb_new, reduced_uvb_old = lonely_hunter(comp_dict[new_gen[i]][nion][ray],
#                                                                                          comp_dict[old_gen[i]][nion][ray])

#                 # calculating generation differences
#                 dens_diff = np.array(reduced_uvb_new["col_dens"]) - np.array(reduced_uvb_old["col_dens"])
                
#                 lonely_new_tot += len(lonely_new)
#                 lonely_old_tot += len(lonely_old)

#                 uvb_dens_diff = np.concatenate((uvb_dens_diff,dens_diff))

#                 # storing old generation
#                 old_gen_dat = np.concatenate((old_gen_dat, reduced_uvb_old["col_dens"]))

#             # each ionization state is plotted separately
#             sns.scatterplot(x=uvb_dens_diff, y=old_gen_dat, 
#                             label = f"{nion}| {old_gen[i]}:{lonely1_tot}, {new_gen[i]}:{lonely2_tot}")

#         plt.xlabel("log10(Column Density Difference) ($cm^{-2}$)", fontsize=12)
#         plt.ylabel(f"{old_gen[i]} Column Density", fontsize=12)
#         plt.title(f"{spec_list[k]} Densities of {old_gen[i]} and {new_gen[i]} for {nrays} Rays", fontsize=15)
#         plt.grid()
#         plt.legend()
#         plt.savefig(comp_paths[f"{old_gen[i]}_{new_gen[i]}"]+f"/dens_{old_gen[i]}_{new_gen[i]}_{spec_list[k]}.pdf")
#         plt.clf()

