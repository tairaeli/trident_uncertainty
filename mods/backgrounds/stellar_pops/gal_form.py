import numpy as np
import unyt as u
from scipy.interpolate import interp1d
import pandas as pd
import os
import argparse
import pickle
import fsps
from astropy.cosmology import FlatLambdaCDM


parser = argparse.ArgumentParser(description = "Pipeline variables and constants for running FSPS")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path where  output data will be stored')
parser.add_argument('--om', nargs='?',action='store', required=True, dest='om_dat', help='Path to Omega+ output data')
parser.add_argument('--pcw', nargs='?',action='store', required=True, dest='pcw_file', help='Path to Putwein et.al. data')
parser.add_argument('--rs', action='store', dest='rs_range', default=[1.2,2.7], type=list, help='Range of redshifts to analyze. Input is a list of 2 values from the lower to upper bound')
parser.add_argument('--d', action='store', dest='d_list', default=[20,50,100,150,200], type=list, help='List of distances from galactic center')

args =parser.parse_args()
dic_args = vars(args)

# loading in array for desired energy bins
pcw_rs = np.genfromtxt(args.pcw_file, max_rows = 1)

# creating a mask to isolate desired redshift data
rs_mask = np.where((pcw_rs>=args.rs_range[0]) & (pcw_rs<=args.rs_range[1]))

# masking putwein redhsifts with desired analysis range
pcw_rs = pcw_rs[rs_mask]

# reading in actual data
pcw_data = np.genfromtxt(args.pcw_file, skip_header=11)

# assigning wave range to variable
pcw_wave = pcw_data[:,0]

pcw_data = pcw_data[:,rs_mask]

# converting wavelengths to energy
Ryd = 2.1798723611035e-18 * u.J
npcw_wave = pcw_wave * u.Angstrom

unq_wv, freq = np.unique(npcw_wave, return_counts=True)
pairs = [(npcw_wave == val).nonzero()[0] for val in unq_wv[freq>1]] # list of arrays w/ indicies

for pair in pairs:
    assert pair.size==2, "More than two duplicates found"
    assert pair[0] < pair[1]

    # adjust wavelenghts by 0.01%
    npcw_wave[pair[0]] -= 0.0001 * npcw_wave[pair[0]]
    npcw_wave[pair[1]] += 0.0001 * npcw_wave[pair[1]]

nu = npcw_wave.to("J", equivalence="spectral") / Ryd
nu = nu.to_value()

    
# rebins spectral data to better match Cloudy's desired input
def rebin(wave,spec,pcw_spec):

    # convert spectral data into luminosities
    lum = spec*(3e8/(wave*1e-10))
    
    # initializing array to store new luminosities
    nlum = np.zeros(len(pcw_wave))
    
    # compare wave data from FSPS to desired wave binning to rebin spec data
    for i in range(pcw_wave.size):
        # accounting for indexing issues in first bin
        if i == 0:
            nlum[i] = sum(lum[np.where((wave<=pcw_wave[i]))])
        else:
            nlum[i] = sum(lum[np.where((wave<=pcw_wave[i])&(wave>pcw_wave[i-1]))])
    
    # converting luminosities back into intensity
    nspec = nlum*((pcw_wave*1e-10)/3e8)
    
    # if type(nspec) != type(pcw_spec):
    #     print(type(nspec))
    #     nspec = nspec.to_value()
    
    # adding the Putwein et.al. intensities to the intensity from FSPS
    nspec += pcw_spec
    
    return nspec
    
# generating stellar population object
sp = fsps.StellarPopulation(zcontinuous=3, imf_type=1, add_agb_dust_model=True,
                        add_dust_emission=True, sfh=3, dust_type=4)

# generating astropy.cosmology object
fl = FlatLambdaCDM(H0=69.5, Om0=0.285, Ob0=0.0461)

# loading in omega+ data and assigning data to variables
om_dat = pd.read_csv(args.om_dat)
ages = np.asarray(om_dat["age"])
sfr_in = np.asarray(om_dat["sfr_in"])
Z = np.asarray(om_dat["metal"])

# setting new floor metalicity based on minimum of stellar population object
nZmin = np.min(sp.zlegend)
Z[Z<=nZmin] = nZmin

# prepping sp object to extract spectrum data
sp.set_tabular_sfh(ages, sfr_in, Z)

# setting lowest alloted intensity
lJ_pad = -50

# TEMPORARY
sp_dat = {}

# iterates through each distance and each redshift to create the 
# input data for cloudy for each distance
for d in args.d_list:
    
    sp_dat[d] = {}
    
    # making a separate directory for each distance
    if not os.path.isdir(args.path+f"/{d}_kpc_dat"):
        os.mkdir(args.path+f"/{d}_kpc_dat")
    
    # do we need a source? either way... I'm just gonna leave this here...
    source = "REDACTED"
    
    # initializing list to store redshift data
    conv_rs = []
    
    # iterating through each redshift
    for irs in range(len(pcw_rs)):
        
        # setting the current redshift
        rs = pcw_rs[irs]
        
        print("Running redshift",str(rs),"at d = "+str(d)+" kpc")
        
        # TEMPORARY
        sp_dat[d][rs] = {}
        
        conv_rs.append(f"{rs:.4e}")
        # conversion of redshift to age (may need more precise value for universe age)
        age = fl.age(rs).to_value()
        
        # generating luminosity at each wavelength
        wave,spec = sp.get_spectrum(tage = age)
        
        # TEMPORARY
        sp_dat[d][rs]["wave"] = wave
        sp_dat[d][rs]["spec"] = spec
        
        # converting distance into centimeters?????
        # d_m = d*3.086e21
        
        # converting spectral luminosities into intensity (solar lum per cm^2)
        d_cm = d*3.086e+21
        spec = spec/(16*np.pi**2*d_cm**2)        
        
        # calling the putwein intensity data for the current redshift
        pcw_spec = pcw_data[:,0,irs]
        
        # rebining data to match desired Cloudy input
        spec = np.log10(rebin(wave,spec,pcw_spec))
        
        # generate interpolation function
        interp = interp1d(nu, spec, fill_value = "extrapolate")
        
        # generating file where data is stored
        fname = args.path+f"/{d}_kpc_dat/z_{rs:.4e}.out"
        with open(fname, "w") as f:
            f.write(f"# {source}\n")
            f.write(f"# z = {rs:.6f}\n")
            f.write("# E [Ryd] log (J_nu)\n")
            
            f.write(f"interpolate ({1e-8:.10f}) ({lJ_pad:.10f})\n")
            f.write(f"continue ({nu[-1]*0.99:.10f}) ({lJ_pad:.10f})\n")
            
            # loop backwards through wavelengths so that lowest energy is first
            for i in range(nu.size-1,-1,-1):
                f.write(f"continue ({nu[i]:.10f}) ({spec[i]:.10f})\n")
                
            f.write(f"continue ({nu[0]*1.01:.10f}) ({lJ_pad:.10f})\n")
            f.write(f"continue ({7.354e6:.10f}) ({lJ_pad:.10f})\n")
            
            x = 10**interp(1)
            f.write(f"f(nu) = {np.log10(x * 4 * np.pi):.10f} at {1:.10f} Ryd\n")
        
    # outputting list of redshifts to another file
    with open('./'+f'd_{d}_kpc_rs.pkl','wb') as f:
        pickle.dump(conv_rs,f)

with open("./spec_dat.pickle","wb") as f:
    pickle.dump(sp_dat,f)
