import numpy as np
import matplotlib.pyplot as plt
import pickle
from uvb_abun_pairwise_compare import *

with open("/mnt/scratch/tairaeli/halo2392_pcw_2019/redshift2.0/data/salsa_out_dict.pickle", "rb") as salsa_dat:
    pcw_dat = pickle.load(salsa_dat)

with open("/mnt/scratch/tairaeli/halo2392_hm_2012/redshift2.0/data/salsa_out_dict.pickle", "rb") as salsa_dat:
    hm_dat = pickle.load(salsa_dat)

with open("/mnt/scratch/tairaeli/halo2392_fg_2009/redshift2.0/data/salsa_out_dict.pickle", "rb") as salsa_dat:
    fg_dat = pickle.load(salsa_dat)

pcw_v_hm = pcw_dat
pcw_v_hm["C II"]["HM_2012"] = hm_dat["C II"]["HM_2012"]
pcw_v_hm["C IV"]["HM_2012"] = hm_dat["C IV"]["HM_2012"]
pcw_v_hm["O VI"]["HM_2012"] = hm_dat["O VI"]["HM_2012"]

ion_list = ["C II", "C IV", "O VI"]
nrays = 4

pcw_v_hm_comp = pairwise_compare(pcw_v_hm, ion_list, nrays)

fig, ax = plt.subplots(3,1, figsize = [8,12])

for i, ion in enumerate(ion_list):
    ax[i].plot(pcw_v_hm_comp[1][ion][2],pcw_v_hm_comp[2][ion][2])
    ax[i].set_xlabel("PCW Log(N)")
    ax[i].set_ylabel("HM Log(N)")
    ax[i].set_title(f"{ion} plot of UVB Column Densities")
    

plt.show()