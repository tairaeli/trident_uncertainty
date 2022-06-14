# from sal_the_snake import *
import salsa
from salsa.utils import check_rays
import numpy as np
import pandas as pd
import argparse
import sys
import os
import matplotlib.pyplot as plt
import yt  
from mpi4py import MPI

comm = MPI.COMM_WORLD

print(f"Let's do some math, kids")

parser = argparse.ArgumentParser(description = "Preliminary constants for SALSA pipeline.")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path where rays and output data will be stored. Directory should contain three other directories called "data", "rays", and "visuals" for program to run smoothly.')
parser.add_argument('--nrays', action='store', dest='nrays', default=4, type=int, help='The number of rays to be generated.')
parser.add_argument('--abun', action='store', dest='file_path', default=argparse.SUPPRESS, help='Path to abundance file, if any. Defaults to solar abundances.')
parser.add_argument('--halo_dir', action='store', dest='halo_dir', default='/mnt/research/galaxies-REU/sims/FOGGIE', help='Path to halo data.')
parser.add_argument('--pat', action='store', dest='pat_lis', default=[2392, 2878, 4123, 5016, 5036,8508], type=list, help='List of different halo pattern file IDs')
parser.add_argument('--rshift', action='store', dest='r_lis', default=[20,18,16], type=list, help='List of different redshift file IDs')

args = parser.parse_args()
dic_args = vars(args)

def generate_names(length, add=''):
        
        	"""
        	Returns a list of generic names for multiplots anda list of generic names for raw data. These lists are typically passed to run_sal()
        	
        	:length: length of lists to be generated. Do not account for indexing starting at zero.
        	
        	:add: Additional relevant information for keeping track of multiplots and data later on. Default add=''
        	"""
        	ndigits = len(str(length))
        	saved_filename_list = []
        	
        	for i in range(length): ##made this so that it would sort correctly for making plots
        		m= i+1
        		n_len = len(str(m))
        		n_zeros = ndigits - n_len
        		k = "0" * n_zeros + str(m)
        		saved_filename_list.append(f'data_AbundanceRow{k}{add}')
        		
        	return saved_filename_list

# iterates through each halo pattern at each redshift
for halo in args.pat_lis:
    for rshift in args.r_lis:
        
        # creating variable names for data bin locations
        path = os.path.expandvars(os.path.expanduser(args.path))
        ray_path = path+'/halo'+f'{halo}'+'/redshift'+f'{rshift}'+'/rays'
        dat_path = path+'/halo'+f'{halo}'+'/redshift'+f'{rshift}'+'/data'
        vis_path = path+'/halo'+f'{halo}'+'/redshift'+f'{rshift}'+'/visuals'
        
        # creating directories for storing data
        os.mkdir(path+'/halo'+f'{halo}')
        os.mkdir(path+'/halo'+f'{halo}'+'/redshift'+f'{rshift}')
        os.mkdir(ray_path) 
        os.mkdir(dat_path)
        os.mkdir(vis_path)
        
        #preliminary shenanigans -- load halo data; define handy variables; plant the seed, as it were
        
        ds = yt.load(f'args.halo_dir/halo_00{halo}/nref11c_nref9f/RD00{rshift}/RD00{rshift}')
        center = ds.arr([23876.757358761424, 23842.452527236022, 22995.717805638298], 'kpc')
        gal_vel = ds.arr([-0.02410298432958413, -136.9259851493111, -147.88486721075668], 'km/s')
        other_fields=['density', 'temperature', 'metallicity', 'radius']
        max_impact=15 #kpc
        units_dict = dict(density='g/cm**3', metallicity='Zsun')
        
        ray_num=f'{0:0{len(str(args.nrays))}d}'
        ray_file=f'{ray_path}/ray{ray_num}.h5'
        
        np.random.seed(11)
        
        #get those rays babyyyy
        
        # CK: Check that rays already exist, and that the have the additional fields contained
        # in the third argument (empty for now; might become a user parameter)
        check = check_rays(ray_path+"/rays", args.nrays, [])
        if not check:
            print("WARNING: rays not found. Generating new ones.")
            salsa.generate_lrays(ds, center.to('code_length'), args.nrays, max_impact, length=600, field_parameters={'bulk_velocity', gal_vel}, ion_list=['H I'], fields=other_fields, out_dir=(ray_path))
        
        ray_list=[]
        for i in range(args.nrays):
            if len(str(i)) != len(str(args.nrays)):
                n = len(str(args.nrays)) - 1
                
                ray_list.append(f'{ray_path}/ray{i: 0{n}d}.h5')
            else:
                ray_list.append(f'{ray_path}/ray{i}.h5')
        
        # CK: Taking a hint from SALSA on how to divvy up the ray list across procs
        ray_arr = np.array(ray_list)
        ray_files_split = np.array_split(ray_arr, comm.size)
        my_rays = ray_files_split[ comm.rank ]
        
        ion_list = ['C II', 'C IV', 'O VI']
        
        # if 'use_SolAb' == False:
        if 'file_path' in dic_args:
        	abun = pd.read_csv(args.file_path, delim_whitespace=True)
        	nrows = len(abun)
        	saved = generate_names(nrows)
        	for row_num in range(nrows):
        		for i in ion_list:
        			abundances = abun.iloc[row_num].to_dict()
        			abs_ext = salsa.AbsorberExtractor(ds, ray_file, ion_name = i, velocity_res =20, abundance_table = abundances, calc_missing=True)
        			df = salsa.get_absorbers(abs_ext, my_rays, method='spice', fields=other_fields, units_dict=units_dict)
        			df.to_csv(f'{dat_path}/{saved[row_num]}_{i.replace(" ", "_")}.txt', sep = ' ')
        			print("Go look at your data!")
        
        else:
        	nrows = 0
        	saved = generate_names(nrows)
        	for i in ion_list:
        		abs_ext = salsa.AbsorberExtractor(ds, ray_file, ion_name = i, abundance_table = None, calc_missing=True)
        		df = salsa.get_absorbers(abs_ext, my_rays, method='spice', fields=other_fields, units_dict=units_dict)
        		df.to_csv(f'{dat_path}/data_SolAb_{i.replace(" ", "_")}.txt', sep = ' ')
        		print("Go look at your data!")