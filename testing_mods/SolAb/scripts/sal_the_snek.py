# from sal_the_snake import *
import salsa
from salsa.utils import check_rays
import numpy as np
import pandas as pd
import argparse
import sys
import os
import matplotlib as plt
import yt  
from mpi4py import MPI

comm = MPI.COMM_WORLD


parser = argparse.ArgumentParser(description = "Preliminary constants for SALSA pipeline.")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path where rays and output data will be stored. Directory should contain three other directories called "data", "rays", and "visuals" for program to run smoothly.')
parser.add_argument('--nrays', action='store', dest='nrays', default=4, type=int, help='The number of rays to be generated.')
parser.add_argument('--abun', action='store', dest='file_path', default=argparse.SUPPRESS, help='Path to abundance file, if any. Defaults to solar abundances.')
parser.add_argument('--halo', action='store', dest='ds_file', default='/mnt/research/galaxies-REU/sims/FOGGIE/halo_002392/nref11c_nref9f/RD0020/RD0020', help='Path to halo data.')

args = parser.parse_args()
dic_args = vars(args)

path = os.path.expandvars(os.path.expanduser(args.path))
if not os.path.exists(path):
	os.makedirs(path)
	os.mkdir(path+"/data")
	os.mkdir(path+"/rays")
if 'file_path' in dic_args:
	# abun = args.file_path
	abun = pd.read_csv(args.file_path, delim_whitespace=True)
	nrows = len(abun)
	litty = 'False'
else:
	# abun = 'No file given. Using solar abundances.'
	nrows = 0
	litty = 'True'

def generate_names(length, add=''):

	"""
	Returns a list of generic names for multiplots anda list of generic names for raw data. These lists are typically passed to run_sal()
	
	:length: length of lists to be generated. Do not account for indexing starting at zero.
	
	:add: Additional relevant information for keeping track of multiplots and data later on. Default add=''
	"""
	
	saved_filename_list = []
	
	for i in range(length):
		saved_filename_list.append(f'data_AbundanceRow{i}{add}')
		
	return saved_filename_list

saved = generate_names(nrows)

#preliminary shenanigans -- load halo data; define handy variables; plant the seed, as it were
ds = yt.load(args.ds_file)
center = ds.arr([23876.757358761424, 23842.452527236022, 22995.717805638298], 'kpc')
other_fields=['density', 'temperature', 'metallicity']
max_impact=15 #kpc
units_dict = dict(density='g/cm**3', metallicity='Zsun')
ion_list = ['H I']

# Need this? YES FOR ABSORBER EXTRACTOR
ray_num = f'{0:0{len(str(args.n_rays))}d}'
ray_file=f'{path}/ray{ray_num}.h5'

np.random.seed(69)

#get those rays babyyyy

# CK: Check that rays already exist, and that the have the additional fields contained
# in the third argument (empty for now; might become a user parameter)
check = check_rays(path, args.n_rays, [])
if not check:
    print("WARNING: rays not found. Generating new ones.")
    salsa.generate_lrays(ds, center.to('code_length'), args.n_rays, max_impact, ion_list=ion_list, fields=other_fields, out_dir=path)

# CK: Consider collect_files from salsa.utils
ray_list=[]
for i in range(args.n_rays):
    if len(str(i)) != len(str(args.n_rays)):
        n = len(str(args.n_rays)) - 1
        
        ray_list.append(f'{path}/ray{i: 0{n}d}.h5')
    else:
        ray_list.append(f'{path}/ray{i}.h5')

# CK: Taking a hint from SALSA on how to divvy up the ray list across procs
ray_arr = np.array(ray_list)
ray_files_split = np.array_split(ray_arr, comm.size)
my_rays = ray_files_split[ comm.rank ]


if 'litty' == False:
	for row_num in range(nrows):
		abundances = abun.iloc[row_num].to_dict()
		# return_df = pd.DataFrame()
		abs_ext_civ = salsa.AbsorberExtractor(ds, ray_file, ion_name = i, abundance_table_args = abundances, recalculate=True)
		df_civ = salsa.get_absorbers(abs_ext_civ, my_rays, method='spice', fields=other_fields, units_dict=units_dict)
		df_civ.to_csv(f'{args.path}/data/{saved[row_num]}.txt', sep = ' ')
		# return_df = return_df.append(df_civ)

else:
	# are abundance table args a dictionary or can it be None?
	abs_ext_civ = salsa.AbsorberExtractor(ds, ray_file, ion_name = i, abundance_table_args = abundances, recalculate=True)
	df_civ = salsa.get_absorbers(abs_ext_civ, my_rays, method='spice', fields=other_fields, units_dict=units_dict)
	df_civ.to_csv(f'{args.path}/data/data_SolAb.txt', sep = ' ')