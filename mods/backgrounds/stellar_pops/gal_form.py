import numpy as np
import unyt as u
from scipy.interpolate import interp1d
import pandas as pd
import os
import argparse
import pickle
import fsps
from astropy.cosmology import FlatLambdaCDM # not sure what cosmology to use
                           # so I just picked the one in the example for astropy


parser = argparse.ArgumentParser(description = "Pipeline variables and constants for running FSPS")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path where  output data will be stored')
parser.add_argument('--om', nargs='?',action='store', required=True, dest='om_dat', help='Path to Omega+ output data')
parser.add_argument('--rs', action='store', dest='rs_list', default=[2,2.5,3], type=list, help='List of redshifts')
parser.add_argument('--d', action='store', dest='d_list', default=[20,50,100,150,200], type=list, help='List of distances from galactic center')

args =parser.parse_args()
dic_args = vars(args)

sp = fsps.StellarPopulation(zcontinuous=3, imf_type=1, add_agb_dust_model=True,
                        add_dust_emission=True, sfh=3, dust_type=4)

om_dat = pd.read_csv(args.om_dat)

# not sure if these parameters are correct. May need to adjust later
fl = FlatLambdaCDM(H0=70, Om0=0.3)

# need Omega+ to run
ages = np.asarray(om_dat["age"])
sfr_in = np.asarray(om_dat["sfr_in"])

Z = np.asarray(om_dat["metal"]) # may need to set floor so fsps doesn't complain

# seems like there are multiple negative Z values here.... why?
# removing all z vals less than 0
nZmin = np.min(Z[Z>0])
Z[Z<=0] = nZmin

# prepping sp object to extract spectrum data
sp.set_tabular_sfh(ages, sfr_in, Z)

# setting age of universe in Gyr
univ_age = 13.8 

# preparing value to be included in Cloud input data
lJ_pad = -50

# iterates through each distance and each redshift to create the 
# input data for cloudy for each distance
for d in args.d_list:
    # making a separate directory for each distance
    if not os.path.isdir(args.path+f"/{d}_kpc_dat"):
        os.mkdir(args.path+f"/{d}_kpc_dat")
    
    # do we need a source? either way... I'm just gonna leave this here...
    source = "Taira"
    
    # initializing list to store redshift data
    conv_rs = []
    
    for rs in args.rs_list:
        conv_rs.append(f"{rs:.4e}")
        # conversion of redshift to age (may need more precise value for universe age)
        age = univ_age - fl.lookback_time(rs).to_value()
        
        # generating luminosity at each wavelength
        wave,spec = sp.get_spectrum(tage = age)
        
        # converting wave data into Angstroms
        Ryd = 2.1798723611035e-18 * u.J
        wave = wave * u.Angstrom
        nu = wave.to("J", equivalence="spectral") / Ryd
        
        # converting spectral luminosities into intensity
        spec = np.log(spec/d**2)
        
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
            for i in range(nu.size-1, -1, -1):
                f.write(f"continue ({nu[i]:.10f}) ({spec[i]:.10f})\n")
                
            f.write(f"continue ({nu[0]*1.01:.10f}) ({lJ_pad:.10f})\n")
            f.write(f"continue ({7.354e6:.10f}) ({lJ_pad:.10f})\n")
            
            x = 10**interp(1)
            f.write(f"f(nu) = {np.log10(x * 4 * np.pi):.10f} at {1:.10f} Ryd\n")
    
    with open(args.path+f'd_{d}_kpc_rs.pkl','wb') as f:
        pickle.dump(conv_rs,f)

