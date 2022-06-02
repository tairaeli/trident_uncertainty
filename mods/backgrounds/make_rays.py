import salsa
from salsa.utils import check_rays
import argparse
import os
import yt 
import numpy as np

parser = argparse.ArgumentParser(description = "Preliminary constants for SALSA pipeline.")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path where rays and output data will be stored. Directory should contain three other directories called "data", "rays", and "visuals" for program to run smoothly.')
parser.add_argument('--nrays', action='store', dest='nrays', default=4, type=int, help='The number of rays to be generated.')
parser.add_argument('--halo', action='store', dest='ds_file', default='/mnt/research/galaxies-REU/sims/FOGGIE/halo_002392/nref11c_nref9f/RD0020/RD0020', help='Path to halo data.')

args = parser.parse_args()

path = os.path.expandvars(os.path.expanduser(args.path))
if not os.path.exists(path+"/data"):
	print(f"Appropriate directories not found on given path. Creating new directories now.")
	# os.makedirs(path)
	os.mkdir(path+"/data")
	os.mkdir(path+"/rays")

#preliminary shenanigans -- load halo data; define handy variables; plant the seed, as it were
ds = yt.load(args.ds_file)
center = ds.arr([23876.757358761424, 23842.452527236022, 22995.717805638298], 'kpc')
other_fields=['density', 'temperature', 'metallicity']
max_impact=15 #kpc
units_dict = dict(density='g/cm**3', metallicity='Zsun')

np.random.seed(11)

# CK: Check that rays already exist, and that the have the additional fields contained
# in the third argument (empty for now; might become a user parameter)
check = check_rays(path, args.nrays, [])
if not check:
    print("WARNING: rays not found. Generating new ones.")
    salsa.generate_lrays(ds, center.to('code_length'), args.nrays, max_impact, ion_list=['H I'], fields=other_fields, out_dir=path+"/rays")


    