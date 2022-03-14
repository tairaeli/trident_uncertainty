#MAY STILL NEED TO ADD READING_FUNC_ARGS TO SALSA.GENERATE_CATALOG

import yt
import salsa
from salsa.utils import check_rays
import numpy as np
import pandas as pd
import matplotlib as plt
# import pdb
from mpi4py import MPI

comm = MPI.COMM_WORLD
	
#print("let's do some math, kids", flush=True)

def sal(ds_file='/mnt/research/galaxies-REU/sims/FOGGIE/halo_002392/nref11c_nref9f/RD0020/RD0020', ray_dir='rays', 
		n_rays=4, center_list=[23876.757358761424, 23842.452527236022, 22995.717805638298], 
		ion_list = ['H I', 'C IV', 'O VI'], df_type = 'cat', **kwargs):
	"""
	Does all the dirty work. 
	Uses yt to load nifty halo dataset. Uses a list of cool ions and SALSA to generate trident LightRay objects and extract absorbers from them. Returns a pandas dataset that contains info on absorbers from singular ray file, many ray files, or a catalog of all ray files with all absorbers. Catalog is the default. 

	:ds_file: Dataset with halo information. Uses HireIsolatedGalaxy by default

	:ray_dir: Directory path where ray.h5 files with be stored

	:n_rays: Number of trident LightRay objects to create. Default is 4

	:ray_num: Used for indexing and saving ray.h5 files so they are not continuously overridden. Default is 0, so saved and/or indexed file will read "ray0.h5"

	:center_list: Defined center of galaxy (I think...). Default is x,y,z = 0.53

	:ion_list: Ions to add to generated spectrum

	:df_type: Preferred data to be returned -- info on absorbers from a single ray file, many ray files, or a catalog of all ray files with all absorbers. Catalog is the default.

	:kwargs: Contains necessary information for returning information oabsorbers frm many ray files (specify as a dictionary named "mult" in kwargs, great place to add an abundance table to be 		 passed to salsa.AbsorberExtractor(), etc.
	"""
	
	if 'mult' in kwargs:
		mult = kwargs['mult']
	else:
		print('mult_args not given')
		mult = {}

	if 'reading_func_args' in kwargs:
		funky_args = kwargs['reading_func_args']
		mult['reading_func_args'] = funky_args
	else:
		print('reading_func_args not given')
		funky_args = {}
	
	
	def mult_salsa(ds, ray_directory, ray_file, units_dict, field, n_rays, ion_list, **mult):
	
		if 'reading_func_args' in mult:
			funky_args = mult['reading_func_args']
		else:
			funky_args = {}
		
		# CK: Consider collect_files from salsa.utils
		ray_list=[]
		for i in range(n_rays):
			if len(str(i)) != len(str(n_rays)):
				n = len(str(n_rays)) - 1
				
				ray_list.append(f'{ray_directory}/ray{i: 0{n}d}.h5')
			# elif len(str(i)) == 2: 
			# 	ray_list.append(f'{ray_directory}/ray0{i}.h5')
			else:
				ray_list.append(f'{ray_directory}/ray{i}.h5')
		
		print(f"RAY LIST: {ray_list}")
		
		# CK: Taking a hint from SALSA on how to divvy up the ray list across procs
		ray_arr = np.array(ray_list)
		ray_files_split = np.array_split(ray_arr, comm.size)
		my_rays = ray_files_split[ comm.rank ]

		return_df = pd.DataFrame()
		
		for i in ion_list:
			abs_ext_civ = salsa.AbsorberExtractor(ds, ray_file, ion_name = i, abundance_table_args = funky_args, recalculate=True)
			df_civ = salsa.get_absorbers(abs_ext_civ, my_rays, method='spice', fields=other_fields, units_dict=units_dict)
			return_df = return_df.append(df_civ)
		
		return return_df
		
	#look for handy args in kwargs
	
	#preliminary shenanigans -- load halo data; define handy variables; plant the seed, as it were
	ds = yt.load(ds_file)
	
	center = ds.arr(center_list, 'kpc')

	other_fields=['density', 'temperature', 'metallicity']
	max_impact=15 #kpc
	units_dict = dict(density='g/cm**3', metallicity='Zsun')

	ray_num = f'{0:0{len(str(n_rays))}d}'
	ray_file=f'{ray_dir}/ray{ray_num}.h5'

	np.random.seed(69)

	#get those rays babyyyy

	# CK: Check that rays already exist, and that the have the additional fields contained
	# in the third argument (empty for now; might become a user parameter)
	check = check_rays(ray_dir, n_rays, [])
	if not check:
		print("WARNING: rays not found. Generating new ones.")
		salsa.generate_lrays(ds, center.to('code_length'), n_rays, max_impact, ion_list=ion_list, fields=other_fields, out_dir=ray_dir)

	#get absorbers -- either many, singular, or catalog returned
	if df_type == 'multiple':
	
		# if len(str(ray_num)) != len(str(n_rays)):
		# 	n = len(str(n_rays)) - 1
		# 	new_ray_num = f'{ray_num: 0{n}d}'
		# else:
		# 	new_ray_num = ray_num
			
		# ray_file=f'{ray_dir}/ray{ray_num}.h5'
		spicy = mult_salsa(ds=ds, ray_directory=ray_dir, ray_file=ray_file, units_dict=units_dict, field=other_fields, n_rays=n_rays, ion_list=ion_list, **mult)

		return spicy
	if df_type == 'single': 
		# ray_file=f'{ray_dir}/ray{ray_num}.h5'
		abs_ext=salsa.AbsorberExtractor(ds, ray_file, ion_name='H I', recalculate=True, **reading_func_args)
		spicy = abs_ext.get_spice_absorbers(other_fields, units_dict=units_dict)

		return spicy
	if df_type == 'cat':
		catalog = salsa.generate_catalog(ds, n_rays, ray_dir, ion_list, fields=other_fields, center = center_list, impact_param_lims=(0, max_impact), method='spice', units_dict=units_dict, reading_func_args = funky_args)

		return catalog
