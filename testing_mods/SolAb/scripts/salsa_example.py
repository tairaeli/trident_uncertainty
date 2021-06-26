print("let's do some math, kids")

def sal(ds_file='HiresIsolatedGalaxy/DD0044/DD0044', ray_dir='more_rays', n_rays=4, ray_num=0, center_list=[0.53, 0.53, 0.53], ion_list=['H I', 'C IV', 'O VI'], df_type = 'multiple', vis=True, trid_rays = False, **kwargs):
	import yt
	import salsa
	import numpy as np
	import pandas as pd
	import matplotlib as plt
	"""
	def make_rays(ds_file, Data_filename, line_list, **rad_rays):
		ds = yt.load(ds_file)
		ray_start = ds.domain_left_edge
		ray_end = ds.domain_right_edge
	
		if 'rad_nums' in rad_rays:
			reading_func_args = rad_rays['rad_nums']
    			field_params = {'reading_func_args':reading_func_args}
    			ray = trident.make_simple_ray(ds, start_position=ray_start, end_position=ray_end, data_filename = Data_Filename, lines = line_list, ftype = 'gas', field_parameters = field_params)
    		
		else:
    			ray = trident.make_simple_ray(ds, start_position=ray_start, end_position=ray_end, data_filename = Data_Filename, lines = line_list, ftype = 'gas')

		return ray
	"""
	
	def visualize(ds_file, center_list, ray_num, name='example_multiplot.png', num_dense_min = 1e-11, num_dense_max=1e-5, **vis_args):
		ray = yt.load(f'{ray_dir}/ray{ray_num}.h5')
		plotter = salsa.AbsorberPlotter(ds_file, ray, "H I", center_gal=center_list, use_spectacle=False, plot_spectacle=False, plot_spice=True, num_dense_max=num_dense_max, num_dense_min=num_dense_min)
		
		fig, axes = plotter.create_multi_plot(outfname=name)
	
	
	def mult_salsa(ds, ray_directory, ray_file, units_dict, field, n_rays, **mult):
		
		ray_list=[]
		for i in range(n_rays):
			ray_list.append(f'{ray_directory}/ray{i}.h5')
		
		abs_ext_civ = salsa.AbsorberExtractor(ds, ray_file, ion_name = 'C II', **reading_func_args)
		df_civ = salsa.get_absorbers(abs_ext_civ, ray_list, method='spice', fields=other_fields, units_dict=units_dict)
		return df_civ
	
	if 'mult' in kwargs:
		mult = kwargs['mult']
	if 'vis_args' in kwargs:
		visual = kwargs['vis_args']
		
	if 'reading_func_args' in kwargs:
		reading_func_args = kwargs['reading_func_args']
	
	ds = yt.load(ds_file)
	
	other_fields=['density', 'temperature', 'metallicity']
	max_impact=15 #kpc
	units_dict = dict(density='g/cm**3', metallicity='Zsun')
	
	np.random.seed(69)
	ray_file=f'{ray_dir}/ray{ray_num}.h5'
	
	if trid_rays is not True:
		salsa.generate_lrays(ds, center_list, n_rays, max_impact, ion_list=ion_list, fields=other_fields, out_dir=ray_dir)
		
	else:
		ray = make_rays(ds_file, Data_filename=ray_file, line_list=ion_list, **rad_rays)
	
	
	if df_type == 'multiple':
		spicy = mult_salsa(ds=ds, ray_directory=ray_dir, ray_file=ray_file, units_dict=units_dict, field=other_fields, **mult)
		if vis == True:
			visualize(ds_file, center_list, ray_num, **vis_args)
		#print(f'SPICY: \n {spicy}')
		return spicy
	if df_type == 'single': 
		abs_ext=salsa.AbsorberExtractor(ds, ray_file, ion_name='H I', **reading_func_args)
		spicy = abs_ext.get_spice_absorbers(other_fields, units_dict=units_dict)
		if vis == True:
			visualize(ds_file, center_list, ray_num, **vis_args)
			
		#print(f'SPICY: \n {spicy}')	
		return spicy
	if df_type == 'cat':
		catalog = salsa.generate_catalog(ds, n_rays, ray_dir, ion_list, fields=other_fields, center = center_list, impact_param_lims=(0, max_impact), method='spice', units_dict=units_dict)
		if vis == True:
			visualize(ds_file, center_list, ray_num, **vis_args)	
		#print(f'CATALOG: \n {catalog}')
		return catalog
