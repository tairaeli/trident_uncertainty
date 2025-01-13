# Adds in ultraviolet ray intnesity emmited from a galaxy to the 
# ultraviolet background at different distances from the galactic 
# center and at different redshifts

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
parser.add_argument('--uvb', nargs='?',action='store', required=True, dest='uvb_file', help='Path to UV ackground data')
parser.add_argument('--rs', action='store', dest='rs_range', default="1.2,2.7", type=list, help='Range of redshifts to analyze. Input is a list of 2 values from the lower to upper bound')
parser.add_argument('--d', action='store', dest='d_list', default=[20,50,100,150,200], type=list, help='List of distances from galactic center in kpc')

args =parser.parse_args()
dic_args = vars(args)

# loading in array for desired energy bins
uvb_rs = np.genfromtxt(args.uvb_file, max_rows = 1)

# creating a mask to isolate desired redshift data
rs_range = args.rs_range.split(",")
rs_mask = np.where((uvb_rs>=rs_range[0]) & (uvb_rs<=rs_range[1]))

# masking redhsifts with desired analysis range using previously generated mask
uvb_rs = uvb_rs[rs_mask]

# reading in actual data
uvb_data = np.genfromtxt(args.uvb_file, skip_header=11)

# assigning wave range to variable
uvb_wave = uvb_data[:,0]

# masking out uvb data with redshift mask
uvb_data = uvb_data[:,rs_mask]

# labeling wave data (Angstroms)
nuvb_wave = uvb_wave * u.Angstrom

# locating pairs of data within wave data that are duplicates
unq_wv, freq = np.unique(nuvb_wave, return_counts=True)
pairs = [(nuvb_wave == val).nonzero()[0] for val in unq_wv[freq>1]] # list of arrays w/ indicies

# goes through each of the pairs of duplicated data and adjusts the wavelengths very slightly
for pair in pairs:
    assert pair.size==2, "More than two duplicates found"
    assert pair[0] < pair[1]

    # adjust wavelenghts by 0.01%
    nuvb_wave[pair[0]] -= 0.0001 * nuvb_wave[pair[0]]
    nuvb_wave[pair[1]] += 0.0001 * nuvb_wave[pair[1]]

# converting wavelengths to energy
Ryd = 2.1798723611035e-18 * u.J
nu = nuvb_wave.to("J", equivalence="spectral") / Ryd
nu = nu.to_value()

# rebins spectral data to better match Cloudy's desired input
def rebin(Egal,erg,uvb_spec):
    
    # initializing array to store new luminosities
    nerg = np.zeros(len(nu))*u.erg/u.cm**2
    
    erg = np.flip(erg)
    
    # compare wave data from FSPS to desired wave binning to rebin spec data
    iuvb = 0
    for i in range(wave.size):
        
        # accounting for indexing issues in first bin
        if Egal[i]<nu[0]:
            nerg[iuvb] += erg[i]
        else:
            if iuvb <= len(nerg)-2:
                if nu[iuvb+1] <= Egal[i]:
                    iuvb += 1
            
                dist = nu[iuvb+1] - nu[iuvb]
                E_dist = Egal[i] - nu[iuvb]
                norm_E_dist = E_dist/dist
                nerg[iuvb] += erg[i]*(1-norm_E_dist)
                nerg[iuvb+1] += erg[i]*norm_E_dist
            else:
                nerg[iuvb] += erg[i]
    print(f"min nerg at {d} is {np.min(nerg)}")
    # adding the Putwein et.al. intensities to the intensity from FSPS
    nerg += uvb_spec
#     print(uvb_spec[minmask])
    return nerg

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
    # TEMPORARY
    sp_dat[d] = {}
    
    # making a separate directory for each distance
    if not os.path.isdir(args.path+f"/{d}_kpc_dat"):
        os.mkdir(args.path+f"/{d}_kpc_dat")
    
    # do we need a source? either way... I'm just gonna leave this here...
    source = "REDACTED"
    
    # initializing list to store redshift data
    conv_rs = []
    
    # iterating through each redshift
    for irs in range(len(uvb_rs)):
        
        # setting the current redshift
        rs = uvb_rs[irs]
        
        print("Running redshift",str(rs),"at d = "+str(d)+" kpc")
        
        # TEMPORARY
        sp_dat[d][rs] = {}
        
        conv_rs.append(f"{rs:.4e}")
        
        # conversion of redshift to age (may need more precise value for universe age)
        age = fl.age(rs).to_value()
        
        # generating solar luminosities per Hz at each wavelength
        wave,spec = sp.get_spectrum(tage = age)
        
        # TEMPORARY
        sp_dat[d][rs]["wave"] = wave
        sp_dat[d][rs]["spec"] = spec
        
        gal_wave = wave*u.Angstrom
        
        nu_gal = gal_wave.to("J", equivalence="spectral") / Ryd
        nu_gal = np.flip(nu_gal.to_value())
        
        # adding units to spectral data
        spec = spec*u.lsun/u.Hz
        
        # converting distance from galactic center into centimeters
        d_cm = (d*u.kpc).to("cm")
        
        # converting FSPS spectra to units of erg/cm^2
        erg = spec.to("erg")/((4*np.pi)**2*d_cm**2)
        
        # calling the uvb intensity data for the current redshift
        uvb_spec = uvb_data[:,0,irs]*u.erg/u.cm**2
        
        # rebining data to match desired Cloudy input
        spec = rebin(nu_gal,erg,uvb_spec)
        
        assert str(spec.units) == "erg/cm**2", f"""UV Intensities are in incorrect 
                                   units. Got {spec.units} when should have 
                                   gotten erg/cm**2"""
        
        spec = np.log10(spec)
        # TEMPORARY
        sp_dat[d][rs]["rebin_spec"] = spec
        
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

# TEMPORARY: outputing fsps data to error check stuffs
with open("./spec_dat.pickle","wb") as f:
    pickle.dump(sp_dat,f)
