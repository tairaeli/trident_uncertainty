from salsa_example import *
import numpy as np

#print("let's do some math, kids")



def visualize(ds_file, center_list, ray_dir, ray_num, name='example_multiplot.png', num_dense_min = 1e-11, num_dense_max=1e-5, **vis_args):
	if len(str(ray_num)) == 1:
		new_ray_num = f'0{ray_num}'
	else:
		new_ray_num=ray_num
		
	ray = yt.load(f'{ray_dir}/ray{new_ray_num}.h5')
	plotter = salsa.AbsorberPlotter(ds_file, ray, "H I", center_gal=center_list, use_spectacle=False, plot_spectacle=False, plot_spice=True, num_dense_max=num_dense_max, num_dense_min=num_dense_min)

	fig, axes = plotter.create_multi_plot(outfname=name)



def change(List, change_index, Change):
	new_list = List
	new_list[change_index] = Change
	
	return new_list



def generate_names(length, add='_'):
	vis_name_list = []
	saved_filename_list = []
	
	for i in range(length):
		vis_name_list.append(f'multiplot{add}{i}')
		saved_filename_list.append(f'data{add}{i}')
		
	return vis_name_list, saved_filename_list



def run_pipe(vis_name_list, saved_filename_list, ray_dir, visual=True, **kwargs):


	def pipe(vis_name, saved_filename, vis_tf, ray_num, pipe_add='_', saved_args = '_', **kwargs):
		if kwargs:
			print('KWARGS PASSED TO PIPE')
			
		if 'pipe_add' in kwargs:
			pipe_add = kwargs['pipe_add']

		if 'saved_add' in kwargs:
			saved_add = kwargs['saved_add']

		new_vis_name = vis_name+pipe_add
		vis_args = dict(name = f'{new_vis_name}.png')

		catalog = sal(ray_dir=ray_dir, ray_num=ray_num, n_rays = 50, **kwargs)
		
		new_saved_filename = saved_filename+saved_add
		catalog.to_csv(f'{new_saved_filename}.txt', sep = ' ')
		catalog.to_csv(f'{new_saved_filename}.csv', sep = ' ')
		
		if vis_tf == True:
			visualize(ds_file='HiresIsolatedGalaxy/DD0044/DD0044', center_list=[0.53, 0.53, 0.53], ray_dir=ray_dir, ray_num=ray_num, **vis_args)
	
	
	for i in range(len(saved_filename_list)):
		pipe(vis_name_list[i], saved_filename_list[i], ray_num=i, vis_tf=visual, **kwargs)
	
	print('go look at your data!')


	
list1 = ['Ne VIII', 'Mg X', 'O VI', 'S IV', 'Si III', 'C II', 'N I']	
list2 = change(list1, 5, 'Fe II')
list3 = change(list2, 3, 'Si IV')
list4 = change(list1, 3, 'Si IV')

ions = np.array([[list1], [list2], [list3], [list4]])

kwargs = dict(reading_func_args=dict(filename='~/git_env/research/oshea/trident_modifications/testing_mods/SolAb/abundances/cgm_abundances_2eb.txt', ratios=False), ray_dir='pipetest_rays')

for i in range(len(ions)):
	kwargs['ion_list'] = list(ions[i][0])
	kwargs['pipe_add'] = f'_ionlist{i}'
	kwargs['saved_add'] = f'_ionlist{i}'
	vis, saved = generate_names(4, add='_ionlist{i}_')
	run_pipe(vis, saved, **kwargs)
