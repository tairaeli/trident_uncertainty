"""
Makes a single data file containing data for four UVBs: Faucher-Gruigere 2009, 
Faucher-Gruigere 2020, Haart and Madau 2012, and Puchwein et al. 2019 for redshifts ___-___.
All intensities the same units of erg/s/cm**2.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd
import unyt as u
import os

def read_cloudy_in(file):
    """
    Reads in cloudy input files as pandas Data Frames

    args:

        file (str) - path to cloudy input file
    
    returns:

        df (Data Frame): the read in / formatted data file
    """
    df = pd.read_csv(file,
                     header = None,
                    skiprows = 5,
                    skipfooter = 3,
                    delim_whitespace = True)
    
    df = df.drop(columns = 0)
    
    for i in range(len(df[1])):
        df[1][i] = float(df[1][i].replace(")","").replace("(",""))
        df[2][i] = float(df[2][i].replace(")","").replace("(",""))
    
    return df

# processing fg09 data
fg09_dir = '/mnt/research/galaxies-REU/tairaeli/fg_2009_uvb_dat'

for filename in os.listdir(fg09_dir):
    rs = 999
    if filename.endswith(".dat"):
        try:
            rs = float(filename[-8:-4])
        except(ValueError):
            rs = float(filename[-7:-4])
        
        fg09_temp = pd.read_csv(fg09_dir+"/"+filename,
                     skiprows=2,delimiter="   ",
                     header=None)
        fg09_ev = fg09_temp[0]*13.605703976
        fg09_int = ((fg09_ev.to_list()*u.eV).to("J")/(u.h)).to_value() # 1/s
        fg09_int = (fg09_temp[1]*10**(-21))*4*np.pi*fg09_int
        fg09 = pd.DataFrame({"E (eV)":fg09_ev,
                             "I (erg/s/cm**2)":fg09_int})

        assert rs != 999, "filename processing failed"

        fg09.to_csv(f"/mnt/scratch/tairaeli/UVB_data_table/fg09/z_{rs}.csv")

# processing the other UVBs
dir_list = ["/mnt/scratch/tairaeli/fg_2020_uvb_dat",
            "/mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/examples/grackle/HM12_UVB",
            "/mnt/scratch/tairaeli/pcw_uvb_dat"]
uvb_names = ["fg20",
             "hm12",
             "pw19"]

for i,dir_loc in enumerate(dir_list):
    for filename in os.listdir(dir_loc):
        rs = 999
        if filename.endswith(".out"):
            rs = float(filename[-14:-4])

            uvb = read_cloudy_in(dir_loc+"/"+filename)
            uvb[1] = uvb[1]*13.605703976
            uvb[3] = ((uvb[1].to_list()*u.eV).to("J")/(u.h)).to_value()
            uvb[4] = (10**uvb[2])*4*np.pi*uvb[3]
            
            out_df = pd.DataFrame({"E (eV)":uvb[1],
                                "I (erg/s/cm**2)":uvb[4]})

            assert rs != 999, "filename processing failed"

            out_df.to_csv(f"/mnt/scratch/tairaeli/UVB_data_table/{uvb_names[i]}/z_{rs}.csv")

# fg20 = read_cloudy_in("/mnt/scratch/tairaeli/fg_2020_uvb_dat/z_2.5000e+00.out")
# fg20[1] = fg20[1]*13.605703976
# fg20[3] = ((fg20[1].to_list()*u.eV).to("J")/(u.h)).to_value()
# fg20[4] = (10**fg20[2])*4*np.pi*fg20[3]

# hm12 = read_cloudy_in("/mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/examples/grackle/HM12_UVB/z_2.5481e+00.out")
# hm12[1] = hm12[1]*13.605703976
# hm12[3] = ((hm12[1].to_list()*u.eV).to("J")/(u.h)).to_value()
# hm12[4] = (10**hm12[2])*4*np.pi*hm12[3]

# pw19 = read_cloudy_in("/mnt/scratch/tairaeli/pcw_uvb_dat/z_2.4790e+00.out")
# pw19[1] = pw19[1]*13.605703976
# pw19[3] = ((pw19[1].to_list()*u.eV).to("J")/(u.h)).to_value()
# pw19[4] = (10**pw19[2])*4*np.pi*pw19[3]

# # reading in fg09 differently so it matches the units of the other UVBs
# fg09_temp = pd.read_csv("/mnt/research/galaxies-REU/tairaeli/fg_2009_uvb_dat/fg_uvb_dec11_z_2.5.dat",
#                      skiprows=2,delimiter="   ",
#                      header=None)
# fg09_ev = fg09_temp[0]*13.605703976
# fg09_int = ((fg09_ev.to_list()*u.eV).to("J")/(u.h)).to_value() # 1/s
# fg09_int = (fg09_temp[1]*10**(-21))*4*np.pi*fg09_int
# fg09 = pd.DataFrame({1:fg09_ev,
#                     "NA2":[0]*261, "NA3":[0]*261,4:fg09_int})


# # masking out range from 10-10^(2.3)
# uvb_list = [fg09,fg20,hm12,pw19]
# for i,uvb in enumerate(uvb_list):
#     uvb_mask = np.where((uvb[1]>10) & (uvb[1]<10**2.3))
#     # print(uvb_mask)
#     uvb = uvb.iloc[uvb_mask]
#     uvb_list[i] = uvb.reset_index()
#     print(uvb.shape)

# out_df = pd.DataFrame({"E_fg09_(eV)":uvb_list[0][1],
#                        "I_fg09_(erg/s/cm**2)":uvb_list[0][4],
#                        "E_fg20_(eV)":uvb_list[1][1],
#                        "I_fg20_(erg/s/cm**2)":uvb_list[1][4],
#                        "E_hm12_(eV)":uvb_list[2][1],
#                        "I_hm12_(erg/s/cm**2)":uvb_list[2][4],
#                        "E_pw19_(eV)":uvb_list[3][1],
#                        "I_pw19_(erg/s/cm**2)":uvb_list[3][4]})

# print(out_df.head())
# out_df.to_csv("/mnt/scratch/tairaeli/UVB_data_table_rs2.5.csv")