# reformat the information from different UV backgrounds to match the format CLOUDY cooling tools prefers
import numpy as np
import unyt as u
from scipy.interpolate import interp1d
import argparse
import pickle
import os


parser = argparse.ArgumentParser(description = "Pipeline variables and constants for running FSPS")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path where  output data will be stored')
parser.add_argument('--uvb', nargs='?',action='store', required=True, dest='uvb_file', help='Path to UV ackground data')
parser.add_argument('--rs', action='store', dest='rs_range', default="1.2,2.7", type=str, help='Range of redshifts to analyze. Input is a list of 2 values from the lower to upper bound')

args = parser.parse_args()
dic_args = vars(args)

# converting rs argument into a list
rs_range = args.rs_range.split(",")

# loading in array for desired energy bins
uvb_rs = np.genfromtxt(args.uvb_file, max_rows = 1)

# creating a mask to isolate desired redshift data
rs_mask = np.where((uvb_rs>=float(rs_range[0])) & (uvb_rs<=float(rs_range[1])))

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

# converting wavelengths to energy
Ryd = 2.1798723611035e-18 * u.J
nu = nuvb_wave.to("J", equivalence="spectral") / Ryd
nu = nu.to_value()
nu = np.round(nu,6)

# locating pairs of data within wave data that are duplicates
unq_nu, freq = np.unique(nu, return_counts=True)
pairs = [(nu == val).nonzero()[0] for val in unq_nu[freq>1]] # list of arrays w/ indicies

# goes through each of the pairs of duplicated data and adjusts the wavelengths very slightly
for pair in pairs:

    assert pair.size==2, "More than two duplicates found"
    assert pair[0] < pair[1]

    # adjust wavelenghts by 0.0001%
    nu[pair[0]] -= 0.000001 * nu[pair[0]]
    nu[pair[1]] += 0.000001 * nu[pair[1]]

# ensuring the nu array is ordered correctly
nu_msk = None
if len(np.where(np.diff(nu)>=0)) > 0:
    nu_msk = np.argsort(nu)[::-1]
    nu = nu[nu_msk]

assert np.where(np.diff(nu)>=0), "nu array resort failed"

# initializing list to store redshift data
conv_rs = []

# setting a source name
source = "NA"

# setting lowest alloted intensity
lJ_pad = -50

# iterating through each redshift
for irs in range(len(uvb_rs)):

    # setting the current redshift
    rs = uvb_rs[irs]
    
    print("Running redshift",f"{rs:.4e}")
    
    conv_rs.append(f"{rs:.4e}")
    
    # calling the uvb intensity data for the current redshift
    spec = uvb_data[:,0,irs]*u.erg/u.cm**2
    
    assert str(spec.units) == "erg/cm**2", f"""UV Intensities are in incorrect 
                                units. Got {spec.units} when should have 
                                gotten erg/cm**2"""
    
    spec = np.log10(spec)

    assert nu_msk is not None, "nu_msk is not being defined correctly"

    spec = spec[nu_msk]

    # generate interpolation function
    interp = interp1d(nu, spec, fill_value = "extrapolate")
    
    # generating file where data is stored
    fname = args.path+f"/z_{rs:.4e}.out"
    
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
    
    assert os.path.exists(fname), "file did not write"
    