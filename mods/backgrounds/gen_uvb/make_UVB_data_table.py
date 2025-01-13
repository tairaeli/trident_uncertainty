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
