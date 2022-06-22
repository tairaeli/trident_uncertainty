import fsps
import time
import numpy as np
import pandas as pd
import pickle
import argparse
from astropy.cosmology import FlatLambdaCDM # Not sure what cosmology to use
                           # So I just picked the one in the example for astropy


parser = argparse.ArgumentParser(description = "Pipeline variables and constants for running FSPS")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path where  output data will be stored')
parser.add_argument('--om', nargs='?',action='store', required=True, dest='om_dat', help='Path to Omega+ output data')
parser.add_argument('--rs', action='store', dest='rs_list', default=[20,18,16], type=list, help='List of redshifts')
parser.add_argument('--d', action='store', dest='d_list', default=[20,50,100,150,200], type=list, help='List of distances from galactic center')

args = parser.parse_args()
dic_args = vars(args)

sp = fsps.StellarPopulation(zcontinuous=3, imf_type=1, add_agb_dust_model=True,
                        add_dust_emission=True, sfh=3, dust_type=4)

om_dat = pd.read_csv(args.om_dat)

# Not sure if these parameters are correct. May need to adjust later
fl = FlatLambdaCDM(H0=70, Om0=0.3)

# Need Omega+ to run
ages = np.asarray(om_dat["age"])
sfr_in = np.asarray(om_dat["sfr_in"])

Z = np.asarray(om_dat["metal"]) # may need to set floor so fsps doesn't complain

# Seems like there are multiple negative Z values here.... why?
# Zmin = np.min(sp.zlegend)

nZmin = np.min(Z[Z>0])

Z[Z<=0] = nZmin

sp.set_tabular_sfh(ages, sfr_in, Z)

univ_age = 13.8 # Gyr

spec_dat = {}
for d in args.d_list:
    spec_dat[d] = {}
    for rs in args.rs_list:
    
        age = univ_age - fl.lookback_time(rs).to_value()
        
        wave,spec = sp.get_spectrum(tage = age)
        
        spec_dat[d][rs] = {'wave':wave,'spec':spec}

with open('broken_spec_dat.pickle','wb') as df:
    pickle.dump(spec_dat,df)