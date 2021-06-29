from salsa_example import *
import numpy as np

print("let's do some math, kids")

list1 = ['O VIII', 'Mg X', 'O VII', 'O VI', 'Ne IV', 'Fe III', 'C II', 'N I']
def change(List, change_index, Change):
	new_list = List
	new_list[change_index] = Change
	return new_list
	
list2 = change(list1, 6, 'Fe II')
list3 = change(list2, 4, 'Si V')
list4 = change(list1, 4, 'Si V')
list5 = change(list1, 4, 'Si IV')
list6 = change(list2, 4, 'Si IV')

ions = np.array([[list1], [list2], [list3], [list4], [list5], [list6]])

def generate_names(length, add='_'):
	vis_name_list = []
	saved_filename_list = []
	for i in range(length):
		vis_name_list.append(f'multiplot{add}{i}')
		saved_filename_list.append(f'data{add}{i}')
	return vis_name_list, saved_filename_list

def run_pipe(vis_name_list, saved_filename_list, ray_dir, **kwargs):

	def pipe(vis_name, saved_filename, ray_num, pipe_add='_', saved_args = '_', **kwargs):
		if 'pipe_add' in kwargs:
			pipe_add = kwargs['pipe_add']
		if 'saved_add' in kwargs:
			saved_add = kwargs['saved_add']
			
		new_vis_name = vis_name+pipe_add
		kwargs['vis_args'] = dict(name = f'{new_vis_name}.png')

		catalog = sal(ray_dir=ray_dir, ray_num=ray_num, n_rays = 50, **kwargs)
		
		
		new_saved_filename = saved_filename+saved_add
		catalog.to_csv(f'{new_saved_filename}.txt', sep = ' ')
		catalog.to_csv(f'{new_saved_filename}.csv', sep = ' ')
	

	for i in range(len(saved_filename_list)):
		pipe(vis_name_list[i], saved_filename_list, ray_num=i, **kwargs)
	
	print('go look at your data!')
	
kwargs = dict(reading_func_args=dict(filename='~/git_env/research/oshea/trident_modifications/testing_mods/SolAb/abundances/cgm_abundances_2eb.txt', ratios=False), ray_dir='pipetest_rays')

for i in range(len(ions)):
	kwargs['ion_list'] = list(ions[i][0])
	kwargs['pipe_add'] = f'_ionlist{i}'
	kwargs['saved_add'] = f'_ionlist{i}'
	vis, saved = generate_names(4, add='_ionlist{i}_')
	run_pipe(vis, saved, **kwargs)

	
"""
Things to chuck in kwargs:
- ray_dir
- reading_func_args
- ions
"""
