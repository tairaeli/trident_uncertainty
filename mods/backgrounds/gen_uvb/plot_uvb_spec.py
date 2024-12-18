"""
Creates a figure showing the intensity of different UVBs at 
different energies
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import seaborn as sn
import scipy as sc
import pandas as pd
import pickle
import unyt as u
import h5py

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

fg2020 = read_cloudy_in("/mnt/scratch/tairaeli/fg_2020_uvb_dat/z_2.5000e+00.out")
fg2020[1] = fg2020[1]*13.605703976
fg2020[3] = ((fg2020[1].to_list()*u.eV).to("J")/(u.h)).to_value()
fg2020[4] = (10**fg2020[2])*4*np.pi*fg2020[3]

hm2012 = read_cloudy_in("/mnt/scratch/tairaeli/hm_2012_uvb_dat/z_2.5481e+00.out")
hm2012[1] = hm2012[1]*13.605703976
hm2012[3] = ((hm2012[1].to_list()*u.eV).to("J")/(u.h)).to_value()
hm2012[4] = (10**hm2012[2])*4*np.pi*hm2012[3]

pcw2019 = read_cloudy_in("/mnt/scratch/tairaeli/pcw_uvb_dat/z_2.4790e+00.out")
pcw2019[1] = pcw2019[1]*13.605703976
pcw2019[3] = ((pcw2019[1].to_list()*u.eV).to("J")/(u.h)).to_value()
pcw2019[4] = (10**pcw2019[2])*4*np.pi*pcw2019[3]

# reading in fg09 differently so it matches the units of the other UVBs
fg2009 = pd.read_csv("/mnt/scratch/tairaeli/fg_2009_uvb_dat/fg_uvb_dec11_z_2.5.dat",
                     skiprows=2,delimiter="   ",
                     header=None)

fg2009_ev = fg2009[0]*13.605703976
fg2009_int = ((fg2009_ev.to_list()*u.eV).to("J")/(u.h)).to_value() # 1/s
fg2009_int = (fg2009[1]*10**(-21))*4*np.pi*fg2009_int

ion_energies = pd.read_csv("./ionization_energies.txt", delimiter = "  ")
imask = ["N_II","N_IV","Ne_VII", "C_II",
         "S_II", "S_III", "S_IV","Mg_II", "Al_III", "O_III", 'Ne_VIII', 'Mg_X']
for ion in imask:
    ion_energies = ion_energies[ion_energies["ion"] != ion]
ion_energies = ion_energies.reset_index()

ion_name_dict = {"H_I":r"H $\mathrm{\i}$",
                 "Si_II":r"Si $\mathrm{\i\i}$",
                 "Si_III":r"Si $\mathrm{\i\i\i}$",
                 "C_III":r"C $\mathrm{\i\i\i}$",
                 "Si_IV":r"Si $\mathrm{\i v}$",
                 "C_IV":r"C $\mathrm{\i v}$",
                 "N_V":r"N $\mathrm{v}$",
                 "O_VI":r"O $\mathrm{v\i}$"}

# creating colors for figure
cmap = colormaps["Dark2"]
colors = cmap(np.linspace(0, 0.7, len(ion_energies["ion"])))

fig, ax = plt.subplots()
trans = ax.get_yaxis_transform()

min_x = 10
max_x = 10**2.3

# plotting vertical lines at ion ionization energies
for i,ion in enumerate(ion_energies["ion"]):
    ax.axvline(ion_energies["ionization energy (eV)"][i],
             linestyle = "-", color = "grey")

# C III
ax.text(ion_energies["ionization energy (eV)"][0]-10**0.63, 10**(-6.7),
                ion_name_dict[ion_energies["ion"][0]], 
                fontsize=12, color = "black", rotation="vertical")

# C IV
ax.text(ion_energies["ionization energy (eV)"][1]-10**0.75, 10**(-7.2),
                ion_name_dict[ion_energies["ion"][1]], 
                fontsize=12, color = "black", rotation="vertical")

# H I
ax.text(ion_energies["ionization energy (eV)"][2]-10**0.09, 10**(-7.2),
                ion_name_dict[ion_energies["ion"][2]], 
                fontsize=12, color = "black", rotation="vertical")

# O VI
ax.text(ion_energies["ionization energy (eV)"][3]-10**1.1, 10**(-7.2),
                ion_name_dict[ion_energies["ion"][3]], 
                fontsize=12, color = "black", rotation="vertical")

# Si II
ax.text(ion_energies["ionization energy (eV)"][4]-10**0.17, 10**(-7.2),
                ion_name_dict[ion_energies["ion"][4]], 
                fontsize=12, color = "black", rotation="vertical")

# Si III
ax.text(ion_energies["ionization energy (eV)"][5]-10**0.47, 10**(-7.2),
                ion_name_dict[ion_energies["ion"][5]], 
                fontsize=12, color = "black", rotation="vertical")

# Si IV
ax.text(ion_energies["ionization energy (eV)"][6]-10**0.57, 10**(-7.2),
                ion_name_dict[ion_energies["ion"][6]], 
                fontsize=12, color = "black", rotation="vertical")

# N V
ax.text(ion_energies["ionization energy (eV)"][7]-10**0.93, 10**(-7.2),
                ion_name_dict[ion_energies["ion"][7]], 
                fontsize=12, color = "black", rotation="vertical")

# plotting spectra
ax.plot(pcw2019[1].astype(np.float),
         pcw2019[4].astype(np.float), label = "PCW 2019")
ax.plot(fg2020[1].astype(np.float),
         fg2020[4].astype(np.float), label = "FG 2020")
ax.plot(fg2009_ev,
         fg2009_int, label = "FG 2009")
ax.plot(hm2012[1].astype(np.float),
         hm2012[4].astype(np.float), label = "HM 2012")

ax.set_xlim(min_x,max_x)
ax.set_ylim(bottom=10**(-7.5))
ax.legend()
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylabel(r"4$\pi$$\nu$$J_{\nu}$ (erg $s^{-1}$ $cm^{-2}$)", fontsize = 16)
ax.set_xlabel("log(E) (eV)", fontsize = 16)
plt.savefig("./uvb_intens_plot_draft.pdf")