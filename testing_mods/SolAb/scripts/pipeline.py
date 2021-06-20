run scripts/salsa_example.py
#kwargs = dict()
#vis_name_list = [] #bunch of names
#saved_filename_list = [] #bunch of associated filenames for saved data

def generate_names(length):
	vis_name_list = []
	saved_filename_list = []
	for i in range(length):
		vis_name_list.append(f'multiplot{i}')
		saved_filename_list.append(f'data{i}')
	return vis_name_list, saved_filename_list

def run_pipe(vis_name_list, saved_filename_list, ray_dir='rays', **kwargs):

	def pipe(vis_name, saved_filename, ray_num, **kwargs):
		kwargs['vis_args'] = dict(name = f'{vis_name}.png')

		catalog = sal(ray_dir=ray_dir, ray_num=ray_num, n_rays = 50, **kwargs)
	
		catalog.to_csv(f'{saved_filename}.txt', sep = ' ')
		catalog.to_csv(f'{saved_filename}.csv', sep = ' ')
	

	for i in range(len(saved_filename_list)):
		pipe(vis_name_list[i], saved_filename_list, ray_num=i, **kwargs)
	
	print('go look at your data!')
	
	
"""
Things to chuck in kwargs:
- ray_dir
- reading_func_args
"""
