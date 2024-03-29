"""
Finds the number of lonely clumps her fractional cutoff and
minimum aborber density
"""
# importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import pickle
import sys

sys.path.insert(1, "/mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/uv_sal/pipeline/")

from uvb_abun_pairwise_compare import pairwise_compare

parser = argparse.ArgumentParser(description = "Select cutoff and UVB for analysis")

parser.add_argument('-ion_list', action='store', 
                    required=False, dest='ion_list', default="O_VI", 
                    help='Ions to be analyzed')

parser.add_argument('-nrays', action='store', default=4,
                    required=False, dest='nrays', type = int,
                    help='Number of rays to be used for analysis')

args = parser.parse_args()

clump_dir = "/mnt/scratch/tairaeli/cutoff_bin/clump_data"

hm_name = "HM_2012"
pcw_name = "PCW_2019"

# frac_cutoff params
low_cut = 0.6
hi_cut = 0.9
cut_stepsize = 0.05

# min_dens params
low_dens = 11.0
hi_dens = 15.0
dens_stepsize = 0.5

ion_list = "/mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/uv_sal/error_handling/SALSA_cutoffs/ion_list.txt"
# initializing ion list
init_ion_list = pd.read_csv(ion_list, sep = " ", header=0)

# adjusting formating of ion list
ion_list = []
for i, ion in enumerate(init_ion_list.keys()):
    alt_ion = ion.replace("_"," ")
    ion_list.append(alt_ion)

nrays = args.nrays

for ray in range(nrays):
    print(ray)
    for i,ion in enumerate(ion_list):
        
        # commented out code for debugging
        # print(ion)
        # if ion != "O VI":
        #     continue

        # initializing directory where data will be stored / pulled from
        ion_dir = "/mnt/scratch/tairaeli/cutoff_bin/ray_"+str(ray)+"/"+str(init_ion_list.keys()[i])+"/"

        # where clump data will be extracted from
        clump_dir = ion_dir+"clump_data/"

        cutoff_lonely_counts = pd.DataFrame(index = np.arange(low_cut,
                                                            hi_cut,
                                                            cut_stepsize),
                                            columns = [hm_name, pcw_name, "num_lonely"])
        # performing analysis on fractional cutoffs
        for i in np.arange(low_cut, hi_cut, cut_stepsize):
            # reading in file for given cut
            with (open(clump_dir+hm_name+"/HM_2012_clump_dat_cutoff_frac_"+str(np.round(i,2))+".pickle", "rb")) as df:
                hm_clumps = pickle.load(df)
            
            # code for debugging
            # if i == 14.0:
            #     if ion == "O VI":
            #         print("WA")

            with (open(clump_dir+pcw_name+"/PCW_2019_clump_dat_cutoff_frac_"+str(np.round(i,2))+".pickle", "rb")) as df:
                pcw_clumps = pickle.load(df)
            
            hm_clumps[ion][pcw_name] = pcw_clumps[ion][pcw_name]

            clump_dat =  hm_clumps

            # perfoming clump classification
            comp_dict = pairwise_compare(clump_dat, nrays)

            cutoff_lonely_counts[hm_name][i] = comp_dict[0][ion][ray][-2]

            cutoff_lonely_counts[pcw_name][i] = comp_dict[0][ion][ray][-1]

            cutoff_lonely_counts["num_lonely"][i] = len(cutoff_lonely_counts[hm_name][i]) + \
                len(cutoff_lonely_counts[pcw_name][i])

        min_dens_lonely_counts = pd.DataFrame(index = np.arange(low_dens,
                                                                hi_dens,
                                                                dens_stepsize),
                                            columns = [hm_name, pcw_name, "num_lonely"])
        # performing analysis on min column densities
        for i in np.arange(low_dens, hi_dens, dens_stepsize):
            # reading in file for given cut
            
            # print(ion)
            with (open(clump_dir+hm_name+"/HM_2012_clump_dat_min_dens_"+str(np.round(i,2))+".pickle", "rb")) as df:
                hm_clumps = pickle.load(df)
            
            with (open(clump_dir+pcw_name+"/PCW_2019_clump_dat_min_dens_"+str(np.round(i,2))+".pickle", "rb")) as df:
                pcw_clumps = pickle.load(df)

            hm_clumps[ion][pcw_name] = pcw_clumps[ion][pcw_name]

            clump_dat =  hm_clumps

            # code for debugging
            # if i == 14.0:
            #     if ion == "O VI":
            #         print("WA")

            comp_dict = pairwise_compare(clump_dat, nrays)

            min_dens_lonely_counts[hm_name][i] = comp_dict[0][ion][ray][-2]

            min_dens_lonely_counts[pcw_name][i] = comp_dict[0][ion][ray][-1]

            min_dens_lonely_counts["num_lonely"][i] = len(min_dens_lonely_counts[hm_name][i]) + \
                len(min_dens_lonely_counts[pcw_name][i])

        best_cutoff = np.argmin(cutoff_lonely_counts["num_lonely"])
        best_min_dens = np.argmin(min_dens_lonely_counts["num_lonely"])

        output = open(ion_dir+"best_params.txt", "w")
        
        output.write("Best cutoff: "+str(cutoff_lonely_counts.index[best_cutoff])+" at "+str(cutoff_lonely_counts))
        output.write("Best min density: "+str(min_dens_lonely_counts.index[best_min_dens])+" at "+str(min_dens_lonely_counts))
        
        output.close()

        # print(cutoff_lonely_counts)
        # print(min_dens_lonely_counts)

        print("Best cutoff:",cutoff_lonely_counts.index[best_cutoff])

        print("Best min density:",min_dens_lonely_counts.index[best_min_dens])