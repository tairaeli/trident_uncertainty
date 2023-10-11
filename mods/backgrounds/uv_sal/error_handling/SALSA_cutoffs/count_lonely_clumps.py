"""
Finds the number of lonely clumps her fractional cutoff and
minimum aborber density
"""
# importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import sys

sys.path.insert(1, "/mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/uv_sal/pipeline/")

from uvb_abun_pairwise_compare import pairwise_compare

input_dir = "/mnt/scratch/tairaeli/cutoff_bin/clump_data"

hm_name = "HM_2012"
pcw_name = "PCW_2019"

# min_dens params
low_dens = 11
hi_dens = 17
dens_stepsize = 0.5

# frac_cutoff params
low_cut = 0.6
hi_cut = 0.9
cut_stepsize = 0.05

ion = "O VI"

nrays = 1

cutoff_lonely_counts = pd.DataFrame(index = np.arange(low_cut,
                                                      hi_cut,
                                                      cut_stepsize),
                                    columns = [hm_name, pcw_name])
# performing analysis on fractional cutoffs
for i in np.arange(low_cut, hi_cut, cut_stepsize):
    # reading in file for given cut
    with (open(input_dir+"/"+hm_name+"/HM_2012_clump_dat_cutoff_frac_"+str(np.round(i,2))+".pickle", "rb")) as df:
        hm_clumps = pickle.load(df)
    
    with (open(input_dir+"/"+pcw_name+"/PCW_2019_clump_dat_cutoff_frac_"+str(np.round(i,2))+".pickle", "rb")) as df:
        pcw_clumps = pickle.load(df)
    
    hm_clumps[ion][pcw_name] = pcw_clumps[ion][pcw_name]

    clump_dat =  hm_clumps

    comp_dict = pairwise_compare(clump_dat, nrays)
    # print(comp_dict[0])
    cutoff_lonely_counts[hm_name][i] = comp_dict[0][ion][0][-2]

    cutoff_lonely_counts[pcw_name][i] = comp_dict[0][ion][0][-1]

min_dens_lonely_counts = pd.DataFrame(index = np.arange(low_cut,
                                                        hi_cut,
                                                        cut_stepsize),
                                      columns = [hm_name, pcw_name])
# performing analysis on min column densities
for i in np.arange(low_dens, hi_dens, dens_stepsize):
    # reading in file for given cut
    with (open(input_dir+"/"+hm_name+"/HM_2012_clump_dat_min_dens_"+str(np.round(i,2))+".pickle", "rb")) as df:
        hm_clumps = pickle.load(df)
    
    with (open(input_dir+"/"+pcw_name+"/PCW_2019_clump_dat_min_dens_"+str(np.round(i,2))+".pickle", "rb")) as df:
        pcw_clumps = pickle.load(df)

    hm_clumps[ion][pcw_name] = pcw_clumps[ion][pcw_name]

    clump_dat =  hm_clumps

    comp_dict = pairwise_compare(clump_dat, nrays)

    min_dens_lonely_counts[hm_name][i] = comp_dict[0][ion][0][-2]

    min_dens_lonely_counts[pcw_name][i] = comp_dict[0][ion][0][-1]

print(np.where(cutoff_lonely_counts.sum(axis=1) == np.max(cutoff_lonely_counts.sum(axis=1))))

print(np.where(min_dens_lonely_counts.sum(axis=1) == np.max(min_dens_lonely_counts.sum(axis=1))))