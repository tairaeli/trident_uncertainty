import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
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

# creating column density plots for each UVB
num_ion = len(ion_list)
fig, ax = plt.subplots((num_ion//2)+(num_ion%2),2, figsize=[8,8])
i = 0
j = 0
for ion in sal_dat.keys():
    for name in uvb_names:
        sns.kdeplot(sal_dat[ion][name]["col_dens"], ax = ax[i,j], label = name)
    ax[i,j].set_xlabel(r"$log_{10}$ $\sigma$ $cm^{-2}$")
    ax[i,j].set_ylabel(r"Density")
    ax[i,j].set_title(ion)
    if (i == 0) and (j==0):
        ax[i,j].legend()

    if (j>0):
        i+=1
        j=0
    else:
        j+=1
    
plt.tight_layout()
plt.grid()
plt.savefig(uvb_dist_path+f"/dist_plot_paper.pdf")
plt.clf()

# making kde for density across entire ray
# fig, ax = plt.subplots((num_ion//2)+(num_ion%2),2, figsize=(10,10))
# colors = sns.color_palette("rocket",3)
# w_size = 10
# j=0
# k=0
# for ion in ion_list:    
#     tot_ray_dens = {}
#     n = len(str(nrays))
#     for i in range(len(uvb_names)):
#         tot_ray_dens[uvb_names[i]] = np.zeros(nrays)

#     for ray in range(nrays):
#         # loading in ray data
#         n = len(str(nrays))
#         ray_dat = yt.load(rs_path+f"/rays/ray{ray:0{n}d}.h5")
#         ray_len = np.max(ray_dat.r[("gas","dl")].to("cm"))

#         for i in range(len(uvb_names)):
#             with open(rs_path +f'/{uvb_names[i]}/ray_dat/{ion}/ray_{ray:0{n}d}_dens.pickle', "rb") as dens_dat:
#                 num_dens = pickle.load(dens_dat)
#             tot_ray_dens[uvb_names[i]][ray] = np.sum(num_dens*ray_len)

#     for i in range(len(old_gen)):
#         avg_dens = np.log10((tot_ray_dens[old_gen[i]] + tot_ray_dens[new_gen[i]])/2)
#         dens_diff = np.log10(tot_ray_dens[new_gen[i]]) - np.log10(tot_ray_dens[old_gen[i]])
        
#         moving_avg = np.convolve(dens_diff, np.ones(w_size) / w_size, mode='valid')
#         x_vals = avg_dens[w_size-1:]

#         sns.scatterplot(x = avg_dens, y = dens_diff, 
#                         label=f"{old_gen[i]}/{new_gen[i]}", 
#                         ax=ax[j,k], color = colors[i])
#         sns.lineplot(x = x_vals, y = moving_avg, ax=ax[j,k], color = colors[i])
#     ax[j,k].set_xlabel(r"log10(Avg Total Column Density)", fontsize=12)
#     ax[j,k].set_ylabel(r"log10($\frac{\sigma_{old}}{\sigma_{new}}$) ($cm^{-2}$)", fontsize=12)
#     ax[j,k].set_title(f"{ion} Total Column Density of {nrays} Rays", fontsize=15)
#     ax[j,k].grid()

#     if (k>0):
#         ax[j,k].legend()
#         j+=1
#         k=0
#     else:
#         k+=1

# plt.tight_layout()
# plt.savefig(uvb_dist_path+f"/tot_col_dens_comp.pdf", bbox_inches="tight")
# plt.clf()

# print("Sorting Data")

# comp_paths = {}

# # making col dens diff histograms for each ion
# for ion in ion_list:

#     nion = ion.replace("_"," ")
#     print(f"Creating {nion}")

#     niter = len(uvb_names)
#     colors = plt.cm.Set1(np.linspace(0,1,int(niter)+1))
#     k = 0
#     # creating pairwise comp plots for each UVB
#     for i in range(len(uvb_names)-1):
#         for j in range(i+1,len(uvb_names)):
#             try:                
#                 with open(rs_path+f"/uvb_compare_{uvb_names[i]}_{uvb_names[j]}.pickle","rb") as file:
#                     comp_dict = pickle.load(file)
#                 with open(rs_path+f"/uvb_compare_{uvb_names[i]}_{uvb_names[j]}.pickle","rb") as file:
#                     ray_perc = pickle.load(file)
#             except:
#                 print(f"Comparison between {uvb_names[i]} and {uvb_names[j]} Does not exist")
#                 continue

#             print(f"Comparing {uvb_names[i]} and {uvb_names[j]}")

#             uvb_dens_diff = np.array([])

#             # iterating over each ray
#             lonely1_tot = 0
#             lonely2_tot = 0

#             for ray in comp_dict[uvb_names[i]][nion].keys():

#                 # removing lonely clumps from comparison
#                 lonely1, lonely2, reduced_uvb1, reduced_uvb2 = lonely_hunter(comp_dict[uvb_names[i]][nion][ray],
#                                                                              comp_dict[uvb_names[j]][nion][ray])
#                 lonely1_tot += len(lonely1["col_dens"])
#                 lonely2_tot += len(lonely2["col_dens"])

#                 dens_diff =  reduced_uvb1["col_dens"] - reduced_uvb2["col_dens"]

#                 uvb_dens_diff = np.concatenate((uvb_dens_diff,dens_diff))
            
#             # creating kde
#             sns.kdeplot(uvb_dens_diff, 
#                         label=f"{uvb_names[i]}/{uvb_names[j]}:{lonely1_tot}:{lonely2_tot}",
#                         warn_singular=False)
    
#     # creating dictionaries to store our data if they don't already exist
#     plt.xlabel(r"log10($\frac{\sigma_{old}}{\sigma_{new}}$) ($cm^{-2}$)", fontsize=12)
#     plt.ylabel("Density", fontsize=12)
#     plt.title(f"{ion} KDE of Column Density Difference for {nrays} Rays", fontsize=15)
#     plt.legend()
#     plt.grid()
#     plt.savefig(uvb_dist_path+f"/{ion}_dens_diff_hist.pdf")
#     plt.clf()

#     # fig 4: scatter of old uvbs vs diff btw them and the new iterations
#     # fig 5-7: density diff vs physical quantities
#     for i in range(len(old_gen)):
#         # making sure that the right comp file is read in
#         try:
#             with open(rs_path+f"/uvb_clump_labels_{old_gen[i]}_{new_gen[i]}.pickle","rb") as file:
#                     clump_categories = pickle.load(file)

#             with open(rs_path+f"/uvb_compare_{old_gen[i]}_{new_gen[i]}.pickle","rb") as file:
#                 comp_dict = pickle.load(file)
#         except:
#             try:
#                 with open(rs_path+f"/uvb_compare_{new_gen[i]}_{old_gen[i]}.pickle","rb") as file:
#                     comp_dict = pickle.load(file)
                
#                 with open(rs_path+f"/uvb_compare_{new_gen[i]}_{old_gen[i]}.pickle","rb") as file:
#                     comp_dict = pickle.load(file)

#             except:
#                 print(f"Comparison between {new_gen[i]} and {old_gen[i]} Does not exist")
#                 continue
        
#         print(f"Comparing {new_gen[i]} and {old_gen[i]}") 

#         # creating dictionaries to store our data if they don't already exist
#         comp_paths[f"{old_gen[i]}_{new_gen[i]}"] = rs_path+f"/{old_gen[i]}_{new_gen[i]}_comp"
#         if os.path.exists(comp_paths[f"{old_gen[i]}_{new_gen[i]}"]) == False:
#             os.mkdir(comp_paths[f"{old_gen[i]}_{new_gen[i]}"])
        
#         ion_path = comp_paths[f"{old_gen[i]}_{new_gen[i]}"]+"/"+ion
#         if os.path.exists(ion_path) == False:
#             os.mkdir(ion_path)

#         # setting up data bins to store data in
#         uvb_dens_diff = np.array([])
#         old_gen_dat = np.array([])
#         new_gen_dat = np.array([])

#         lonely_new_tot = 0
#         lonely_old_tot = 0

#         # this doesn't feel like a good way to store our data
#         phys_quant = {}
#         phys_quant["avg_ray_dens"] = {"mean":np.array([]),
#                                       "upper":np.array([]),
#                                       "lower":np.array([])}
#         phys_quant["avg_ray_temp"] = {"mean":np.array([]),
#                                       "upper":np.array([]),
#                                       "lower":np.array([])}
#         phys_quant["avg_ray_met"] = {"mean":np.array([]),
#                                       "upper":np.array([]),
#                                       "lower":np.array([])}
#         # temporary: adding color array
#         # color_arr = np.array([])
#         c = 0
        
#         color_arr = []
#         # old_lone = []
#         # new_lone = []
#         for ray in comp_dict[old_gen[i]][nion].keys():
#             old_lone = []
#             new_lone = []
            
#             # creating colorbar for clump categories
#             match, shorter, longer, overlap, split, merge, lonely_old, lonely_new = clump_categories[nion][ray]
#             print(ray, clump_categories[nion][ray], old_gen[i], new_gen[i], len(comp_dict[old_gen[i]][nion][ray]["col_dens"]))
            
#             clumps_containted_1 = 0
#             clumps_containted_2 = 0
#             # some weird edge case causing me to ad a +2
#             # print(comp_dict.keys())
#             for j in range(len(comp_dict[old_gen[i]][nion][ray]["col_dens"])+2):
#                 if (j) in split:
#                     clumps_containted_2 += split[j]
                
#             for j in range(len(comp_dict[new_gen[i]][nion][ray]["col_dens"])+2): 
#                 if (j) in merge:
#                     clumps_containted_1 += merge[j]

#             for j in range(len(comp_dict[old_gen[i]][nion][ray]["col_dens"])+clumps_containted_1):
#                 if (j) in match:
#                     color_arr.append("match")

#                 elif (j) in split:
#                     color_arr.append("merge")

#                 elif (j) in lonely_old:
#                     old_lone.append(j)
                
#                 elif (j) in longer:
#                     color_arr.append("diff_size")
                
#                 elif (j) in shorter:
#                     color_arr.append("diff_size")

#                 elif (j) in overlap:
#                     color_arr.append("overlap")
                
#                 else:
#                     print("Error: Old",ray,j)

#             for j in range(len(comp_dict[new_gen[i]][nion][ray]["col_dens"])+clumps_containted_2): 
#                 if (j) in lonely_new:
#                     new_lone.append(j)
#                 elif (j) in merge:
#                     color_arr.append("merge")
#                 elif (j) in match:
#                     continue
#                 else:
#                     print("Error: New", ray,j)

#             # removing lonely clumps from comparison
#             lonely_new, lonely_old, reduced_uvb_new, reduced_uvb_old = lonely_hunter(comp_dict[new_gen[i]][nion][ray],
#                                                                                      comp_dict[old_gen[i]][nion][ray])
            
#             assert len(lonely_new["col_dens"]) == len(new_lone), "NEW FAIL "+str(len(lonely_new))+":"+str(len(new_lone))
#             assert len(lonely_old["col_dens"]) == len(old_lone), "OLD FAIL "+str(len(lonely_old))+":"+str(len(old_lone))
#             assert len(reduced_uvb_old["col_dens"]) == len(reduced_uvb_new["col_dens"]), "Arrays are different sizes. \n old = "+str(len(reduced_uvb_old["col_dens"]))+ "\n new = "+str(len(reduced_uvb_new["col_dens"]))

#             # keeping track of total number of lonely clumps
#             lonely_new_tot += len(new_lone)
#             lonely_old_tot += len(old_lone)

#             dens_diff =  reduced_uvb_new["col_dens"] - reduced_uvb_old["col_dens"]

#             uvb_dens_diff = np.concatenate((uvb_dens_diff,dens_diff))

#             # storing old/new generation
#             old_gen_dat = np.concatenate((old_gen_dat, reduced_uvb_old["col_dens"]))
#             new_gen_dat = np.concatenate((new_gen_dat, reduced_uvb_new["col_dens"]))

#             assert len(uvb_dens_diff) == len(color_arr), "Arrays diff sizes at ray "+str(ray)+". uvb_dens_diff:"+str(len(uvb_dens_diff))+" color_arr:"+str(len(color_arr))+" "+old_gen[i]+" "+new_gen[i]

#             # adding upper and lower bounds
#             prop_list = ["density","temperature","metallicity"]
#             new_old_col_dens = np.array([reduced_uvb_old["col_dens"],
#                                          reduced_uvb_new["col_dens"]])

#             for k, prop in enumerate(prop_list):
                
#                 ub = np.zeros_like(reduced_uvb_old["col_dens"])
#                 lb = np.zeros_like(reduced_uvb_old["col_dens"])
                
#                 old = reduced_uvb_old[prop]
#                 new = reduced_uvb_new[prop]
#                 key = list(phys_quant.keys())[k]

#                 ub_mask = np.where(old>new)
#                 lb_mask = np.where(old<new)

#                 ub[ub_mask] = old[ub_mask]
#                 ub[lb_mask] = new[lb_mask]

#                 lb[lb_mask] = old[lb_mask]
#                 lb[ub_mask] = new[ub_mask]

#                 old_new_arr = np.array([new,old])
#                 phys_quant[key]["mean"] = np.concatenate((phys_quant[key]["mean"],
#                                                           np.average(old_new_arr, axis = 0,
#                                                                      weights = new_old_col_dens)))
#                 phys_quant[key]["lower"] = np.concatenate((phys_quant[key]["lower"], lb))
#                 phys_quant[key]["upper"] = np.concatenate((phys_quant[key]["upper"], ub))

#         # setting labels for clump categories
#         clump_cat_labels = ["match","diff_size","overlap","merge"]
#         palt = plt.cm.tab10(np.linspace(0,1,4))

#         # casting color_arr as a numpy array
#         color_arr = np.array(color_arr)
        
#         # old col dens vs col dens diff
#         # Line2D([0], [0], color='w', lw=4, label=f"{old_gen[i]}:{lonely_old_tot}\n {new_gen[i]}:{lonely_new_tot}")
#         # for k,cat in enumerate(clump_cat_labels):
#         #     cat_mask = np.where(color_arr == cat)
#         #     plt.plot(old_gen_dat[cat_mask], uvb_dens_diff[cat_mask],
#         #                 linestyle = "None", marker = ".", color = palt[k], label = cat)
#         # plt.xlabel(fr"{old_gen[i]} log10(Column Density) ($cm^{-2}$)", fontsize=18)
#         # plt.ylabel(r"log10($\frac{\sigma_{new}}{\sigma_{old}}$) ($cm^{-2}$)", fontsize=18)
#         # plt.title(f"{ion} Clump By Clump Density Differences of {old_gen[i]} for {nrays} Rays", fontsize=19)
#         # plt.grid()
#         # plt.legend(fontsize="large")
#         # plt.savefig(ion_path+f"/dens_diff_{old_gen[i]}_{ion}.pdf")
#         # plt.clf()

#         # setting up subplots
#         fig, ax = plt.subplots(2,2, figsize = (16,16))
#         ax_lab_size = 15
#         title_size = 15

#         # direct column density comparison
#         match_line = np.linspace(min(np.min(new_gen_dat),np.min(old_gen_dat)),max(np.max(new_gen_dat),np.max(old_gen_dat)))
#         sns.lineplot(x=match_line, y=match_line, color = "black", linestyle = '--', ax = ax[0,0])
#         for k,cat in enumerate(clump_cat_labels):
#             cat_mask = np.where(color_arr == cat)
#             ax[0,0].plot(old_gen_dat[cat_mask], new_gen_dat[cat_mask],
#                         linestyle = "None", marker = ".", color = palt[k], label = cat)
#         ax[0,0].set_xlabel(fr"{old_gen[i]} log10(Column Density) ($cm^{-2}$)", fontsize=ax_lab_size)
#         ax[0,0].set_ylabel(fr"{new_gen[i]} log10(Column Density) ($cm^{-2}$)", fontsize=ax_lab_size)
#         ax[0,0].set_title(f"{ion} Absorber Column Density Compare: {nrays} Rays", fontsize=title_size)
#         ax[0,0].grid()

#         # gas density vs col dens diff
#         legend_labs = [Line2D([0], [0], color='w', lw=4, label=f"{old_gen[i]}:{lonely_old_tot}\n {new_gen[i]}:{lonely_new_tot}")]
#         for k,cat in enumerate(clump_cat_labels):
#             cat_mask = np.where(color_arr == cat)
#             ax[0,1].plot(np.log10(phys_quant["avg_ray_dens"]["mean"][cat_mask]), uvb_dens_diff[cat_mask],
#                         linestyle = "None", marker = ".", color = palt[k], label = cat)
#             ax[0,1].hlines(y=uvb_dens_diff[cat_mask],xmin=np.log10(phys_quant["avg_ray_dens"]["lower"][cat_mask]), 
#                     xmax=np.log10(phys_quant["avg_ray_dens"]["upper"][cat_mask]),
#                     color = palt[k])
            
#             legend_labs.append(Line2D([0], [0], color=palt[k], lw=4, label=cat))
#         ax[0,1].set_xlabel("log10 Mean Gas Density ($cm^{-3}$)", fontsize=ax_lab_size)
#         ax[0,1].set_title(f"{ion} Mean Gas Densities of {old_gen[i]} and {new_gen[i]}: {nrays} Rays", fontsize=title_size)
#         ax[0,1].grid()
#         ax[0,1].legend(handles=legend_labs)

#         # temp vs col dens diff
#         for k,cat in enumerate(clump_cat_labels):
#             cat_mask = np.where(color_arr == cat)
#             ax[1,0].plot(np.log10(phys_quant["avg_ray_temp"]["mean"][cat_mask]), uvb_dens_diff[cat_mask],
#                         linestyle = "None", marker = ".", color = palt[k], label = cat)
#             ax[1,0].hlines(y=uvb_dens_diff[cat_mask],xmin=np.log10(phys_quant["avg_ray_temp"]["lower"][cat_mask]), 
#                     xmax=np.log10(phys_quant["avg_ray_temp"]["upper"][cat_mask]),
#                     color = palt[k])
#         ax[1,0].set_ylabel(r"log10($\frac{\sigma_{old}}{\sigma_{new}}$) ($cm^{-2}$)", fontsize=ax_lab_size)
#         ax[1,0].set_xlabel("Mean Temperature (K)", fontsize=ax_lab_size)
#         ax[1,0].set_title(f"{ion} Mean Temperature of {old_gen[i]} and {new_gen[i]}: {nrays} Rays", fontsize=title_size)
#         ax[1,0].grid()
        
#         # metallicity vs col dens diff
#         for k,cat in enumerate(clump_cat_labels):
#             cat_mask = np.where(color_arr == cat)
#             ax[1,1].plot(np.log10(phys_quant["avg_ray_met"]["mean"][cat_mask]), uvb_dens_diff[cat_mask],
#                         linestyle = "None", marker = ".", color = palt[k], label = cat)
#             ax[1,1].hlines(y=uvb_dens_diff[cat_mask],xmin=np.log10(phys_quant["avg_ray_met"]["lower"][cat_mask]), 
#                     xmax=np.log10(phys_quant["avg_ray_met"]["upper"][cat_mask]),
#                     color = palt[k])
#         ax[1,1].set_ylabel(r"log10($\frac{\sigma_{old}}{\sigma_{new}}$) ($cm^{-2}$)", fontsize=ax_lab_size)
#         ax[1,1].set_xlabel("Mean Metallicity (K)", fontsize=ax_lab_size)
#         ax[1,1].set_title(f"{ion} Mean Metallicity of {old_gen[i]} and {new_gen[i]}: {nrays} Rays", fontsize=title_size)
#         ax[1,1].grid()

#         # hardcoding some lines in because latex and variable strings do not play well together
#         if i==0:
#             ax[0,1].set_ylabel(r"log_{10}($\frac{\sigma_{FG09}}{\sigma_{FG20}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
#             ax[1,0].set_ylabel(r"log_{10}($\frac{\sigma_{FG09}}{\sigma_{FG20}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
#             ax[1,1].set_ylabel(r"log_{10}($\frac{\sigma_{FG09}}{\sigma_{FG20}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
#         elif i==1:
#             ax[0,1].set_ylabel(r"log_{10}($\frac{\sigma_{HM19}}{\sigma_{PCW19}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
#             ax[1,0].set_ylabel(r"log_{10}($\frac{\sigma_{HM19}}{\sigma_{PCW19}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
#             ax[1,1].set_ylabel(r"log_{10}($\frac{\sigma_{HM19}}{\sigma_{PCW19}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
#         elif i==2:
#             ax[0,1].set_ylabel(r"log_{10}($\frac{\sigma_{FG20}}{\sigma_{PCW19}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
#             ax[1,0].set_ylabel(r"log_{10}($\frac{\sigma_{FG20}}{\sigma_{PCW19}})$) ($cm^{-2}$)", fontsize=ax_lab_size)
#             ax[1,1].set_ylabel(r"log_{10}($\frac{\sigma_{FG20}}{\sigma_{PCW19}})$) ($cm^{-2}$)", fontsize=ax_lab_size)

#         # figure settings
#         plt.tight_layout()
#         plt.savefig(ion_path+f"/phys_quantities_{old_gen[i]}_{new_gen[i]}_{ion}.pdf")
#         plt.show()
#         plt.clf()


# # making comparison plots for each ray
# colors = plt.cm.cool([0.3,0.4,0.9,1])

# for ion in ion_list:
#     nion = ion.replace("_"," ")
#     for ray in range(nrays):
#         # loading in ray data
#         n = len(str(nrays))
#         ray_dat = yt.load(rs_path+f"/rays/ray{ray:0{n}d}.h5")
#         ray_pos = ray_dat.r[("gas","l")].to("kpc")

#         for i in range(len(old_gen)):

#             plot_path = rs_path+f"/{old_gen[i]}_{new_gen[i]}_comp/{ion}/ray_dat"
#             if os.path.exists(plot_path) == False:
#                 os.mkdir(plot_path)

#             old = sal_dat[nion][old_gen[i]][sal_dat[nion][old_gen[i]]["lightray_index"] == f"{ray:0{n}d}"].reset_index()
#             new = sal_dat[nion][new_gen[i]][sal_dat[nion][new_gen[i]]["lightray_index"] == f"{ray:0{n}d}"].reset_index()

#             with open(rs_path +f'/{old_gen[i]}/ray_dat/{ion}/ray_{ray:0{n}d}_dens.pickle', "rb") as dens_dat:
#                     old_num_dens = pickle.load(dens_dat)
            
#             with open(rs_path +f'/{new_gen[i]}/ray_dat/{ion}/ray_{ray:0{n}d}_dens.pickle', "rb") as dens_dat:
#                     new_num_dens = pickle.load(dens_dat)
            
#             plt.figure(figsize = [15,8], dpi = 500, facecolor = "white")
#             plt.semilogy(ray_pos, old_num_dens,label = old_gen[i], color = colors[3])
#             plt.semilogy(ray_pos, new_num_dens,label = new_gen[i], color = colors[1])

#             for j in range(len(old["interval_start"])):
                
#                 lb = old["interval_start"][j]
#                 hb = old["interval_end"][j]
#                 rng = [lb,hb]
                
#                 yb = old_num_dens[slice(*rng)]
#                 xb = ray_pos[lb:hb]
                
#                 plt.fill_between(xb,yb, color = colors[2], alpha = 0.3)
            
#             for j in range(len(new["interval_start"])):
#                 lb = new["interval_start"][j]
#                 hb = new["interval_end"][j]
#                 rng = [lb,hb]
                
#                 yb = new_num_dens[slice(*rng)]
#                 xb = ray_pos[lb:hb]
                
#                 plt.fill_between(xb,yb, color = colors[0], alpha = 0.3)
            
#             plt.grid()
#             plt.legend(fontsize = 20)
#             # plt.ylim(10**(-15),10**(-6.5))
#             # plt.xlim(500,900)
#             plt.xlabel("Position Along Ray (kpc)", fontsize = 25)
#             plt.ylabel(r"log10(Density) ($cm^{-3}$)", fontsize = 25)
#             # plt.ylim(1e-12, 2e-6)
#             plt.title(f"UVB Number Density Comparison "+ion, fontsize = 30)
#             plt.savefig(plot_path+f"/UVB_dens_compare_ray_{ray:0{n}d}.pdf")
#             plt.clf()