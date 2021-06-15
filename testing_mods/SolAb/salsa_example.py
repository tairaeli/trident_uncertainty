print("let's do some math, kids")

def single_sal(ds_file='IsolatedGalaxy/galaxy0030/galaxy0030', ray_dir='my_rays', n_rays=4, ray_num=0, shii=18, center_list=[0.53, 0.53, 0.53], single=True, **mult):
	import yt
	import salsa
	import numpy as np
	import pandas as pd
	
	def mult_salsa(ray_directory, ray_file, ray_list_num = 4, **mult):
		
		ray_list=[]
		for i in range(ray_list_num):
			ray_list.append(f'{ray_dir}/ray{i}.h5')
		
		abs_ext_civ = salsa.AbsorberExtractor(ds, ray_file, ion_name = 'C IV')
		df_civ = salsa.get_absorbers(abs_ext_civ, ray_list, method='spice', fields=other_fields, units_dict=units_dict)
		return df_civ
	
	ds = yt.load(ds_file)
	
	
	ion_list=['H I', 'C IV', 'O VI']
	other_fields=['density', 'temperature', 'metallicity']
	max_impact=15 #kpc
	units_dict = dict(density='g/cm**3', metallicity='Zsun')
	
	np.random.seed(shii)
	salsa.generate_lrays(ds, center_list, n_rays, max_impact, ion_list=ion_list, fields=other_fields, out_dir=ray_dir)
	
	ray_file=f'{ray_dir}/ray{ray_num}.h5'
	
	if single is not True:
		return mult_salsa(ray_directory=ray_dir, ray_file=ray_file, **mult)
	
	else:
		
		abs_ext=salsa.AbsorberExtractor(ds, ray_file, ion_name='H I')
		return abs_ext.get_spice_absorbers(other_fields, units_dict=units_dict)


