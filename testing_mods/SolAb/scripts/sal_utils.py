from sal_the_snake import *
import numpy as np
import pandas as pd
import argparse
import sys

def strikie(text):
    result = ''
    for c in text:
        result = result + c + '\u0336'
    return result


print(f"Let's do some {strikie('meth ')}math, kids")

parser = argparse.ArgumentParser(description = "Preliminary constants for SALSA pipeline.")
parser.add_argument('-ds', nargs='?', action='store', required=True, dest='path', help='Path where rays and output data will be stored. Directory should contain three other directories called "data", "rays", and "visuals" for program to run smoothly.')
parser.add_argument('--nrays', action='store', dest='nrays', default=4, type=int, help='The number of rays to be generated.')
parser.add_argument('--abun', action='store', dest='file_path', default=argparse.SUPPRESS, help='Path to abundance file, if any. Defaults to solar abundances.')
# ask claire about how to use nsc flag without interaction (want to specify just the halo)
# parser.add_argument('-nsc', action='store_true', dest='nsc', help='nsc stands for "not super computer". Use this flag if you would like to manually type out the entire path to your halo data. Otherwise, the interpreter assumes you are using the HPCC FOGGIE simulation data.')
parser.add_argument('--halo', action='store', dest='ds_file', default='/mnt/research/galaxies-REU/sims/FOGGIE/halo_002392/nref11c_nref9f/RD0020/RD0020', help='Path to halo data.')

#args = parser.parse_args()
args = parser.parse_args()
dic_args = vars(args)

# define some variables 
path_list = ['/mnt/research/fuhrmane/test_sal/', args.path]
path = ''.join(path_list)
if 'file_path' in dic_args:
	abundances = args.file_path
	nrows = args.nrays
	litty = 'False'
else:
	abundances = 'No file given. Using solar abundances.'
	df = pd.read_csv(args.file_path, delim_whitespace=True)
	nrows = len(df)
	litty = 'True'

cont = input(f"RAY AND OUTPUT DATA WILL BE STORED HERE: {path}\nNUMBER OF RAYS TOO BE GENERATED: {args.nrays}\nUSING HALO DATA CONTAINED HERE: {args.ds_file}\nABUNDANCE DATA TO BE USED: {abundances}\n\ncontinue?[y/n] ")
if cont == 'y':
    print("yeet")

else:
    sys.exit("Pull your shit together, Becky.")

# path = '/mnt/home/fuhrmane/test_sal/test3/'

def visualize(ds_file, center_list, ray_dir, n_rays, ray_num, ion='O VI', name='example_multiplot', num_dense_min = 1e-11, num_dense_max=1e-5, **vis_args):
	
	"""
	Uses salsa.AbsorberPlotter() to generate multiplots of data produced by salsa.AbsorberExtractor()
	
	:ds_file: Halo dataset used in generating data
	
	:center_list: x, y, and z coordinates of center of galaxy on plot
	
	:ray_dir: Directory where ray.h5 files are stored (same one that's passed to sal())
	
	:ray_num: Same as ray_num from run_sal(). Here, it is used to index ray.h5 file for information
	
	:name: String used to name multiplot file
	
	:num_dense_min: To be passed to salsa.AbsorberPlotter(). Default 1e-11
	
	:num_dense_max: To be passed to salsa.AbsorberPlotter(). Default 1e-5
	
	:vis_args: Mainly used by run_sal() to pass a multiplot name other than "example_multiplot.png". Not meant to be manually passed by user. 
	"""
	
	# if len(str(ray_num)) == 1:
	# 	new_ray_num = f'0{ray_num}'
	# else:
	# 	new_ray_num=ray_num
	
	# newname = f'{name}_ray{ray_num}.png'
	
	if len(str(ray_num)) != len(str(n_rays)):
		n = len(str(n_rays)) - 1
		new_ray_num = f'{ray_num: 0{n}d}'
		
	else: 
		new_ray_num = ray_num

	newname = f'{name}_ray{new_ray_num}.png'	
	
	ray = yt.load(f'{ray_dir}/ray{new_ray_num}.h5')
	plotter = salsa.AbsorberPlotter(ds_file, ray, ion, center_gal=center_list, use_spectacle=False, plot_spectacle=False, plot_spice=True, num_dense_max=num_dense_max, num_dense_min=num_dense_min)

	fig, axes = plotter.create_multi_plot(outfname=newname)



def update_ionlist(list_name, change_index, ion):
	
	"""
	Used to more efficiently generate lists of ions to be passed to run_sal() (later passed to sal())
	
	:list_name: Original list to be updated
	
	:change_index: Integer -- index of ion to be replaced
	
	:ion: Ion that will replace the old ion
	"""
	
	new_list = list_name.copy()
	new_list[change_index] = ion
	
	return new_list



def generate_names(length, add='_'):

	"""
	Returns a list of generic names for multiplots anda list of generic names for raw data. These lists are typically passed to run_sal()
	
	:length: length of lists to be generated. Do not account for indexing starting at zero.
	
	:add: Additional relevant information for keeping track of multiplots and data later on. Default add=''
	"""
	
	vis_name_list = []
	saved_filename_list = []
	
	for i in range(length):
		vis_name_list.append(f'multiplot_row{i}{add}')
		saved_filename_list.append(f'data_row{i}{add}')
		
	return vis_name_list, saved_filename_list
	
def run_sal(vis_name, saved_filename, vis_tf, ray_dir, path, n_rays, ds_file, vis_add='_', saved_add = '_', **kwargs):
	
	"""
	Calls sal() and visualize() functions and saves data. 
	
	:vis_name: Indexed from vis_name_list and used to name generated multiplot of data according to iteration number (and other naming characteristics which are passed through vis_add)
	
	:saved_filename: Indexed from saved_filename_list and used to name saved data according to iteration number (and other naming characteristics which are pass through saved_add)
	
	:vis_tf: If True, uses visualize() to generate multiplots to accompany data
	
	:ray_num: Used to keep track of iteration value from for loop in run_sal() to aid in naming multiplots and data
	
	:path: Path to directories where data is stored. Default path=''
	
	:vis_add: Used to add relevant information to the names of generated multiplots. Should be passed through **kwargs.
	
	:saved_add: Used to add relevant information to the names of generated data. Should be passed through **kwargs.
	
	:kwargs: See description from run_sal()
	"""
	
	#if 'ion_list' in kwargs:
	#	print("ION LIST GIVEN")
	#	ion_list = kwargs['ion_list']
	#else:
	#	ion_list = None
	#print(f"ION LIST: {kwargs['ion_list']}")
	print(f'KWARGS: {kwargs}')



	catalog = sal(ds_file=ds_file, ray_dir=ray_dir, n_rays = n_rays, df_type='multiple', **kwargs)

	new_saved_filename = saved_add+saved_filename
	catalog.to_csv(f'{path}{new_saved_filename}.txt', sep = ' ')
	catalog.to_csv(f'{path}{new_saved_filename}.csv', sep = ' ')

	if vis_tf == True:
		new_vis_name = vis_add+vis_name
		vis_args = dict(name = f'{path}{new_vis_name}')
		for r in range(n_rays):
			visualize(ds_file=ds_file, center_list=[0.53, 0.53, 0.53], ray_dir=ray_dir, ray_num=r, **vis_args)

	print('go look at your data!')


	
list1 = ['Ne VIII', 'Mg X', 'O VI', 'S IV', 'Si III', 'C II', 'N I']	
vis, saved = generate_names(nrows)

#list2 = update_ionlist(list1, 5, 'Fe II')
#list3 = update_ionlist(list2, 3, 'Si IV')
#list4 = update_ionlist(list1, 3, 'Si IV')

#ions = np.array([[list1], [list2], [list3], [list4]])

if litty == 'True':
	kwargs = dict(reading_func_args = dict(filename = args.file_path, ratios=False), ray_dir=f'{path}rays', ion_list = list1)
	selr = 0
	for i in range(nrows):
		#kwargs['ion_list'] = list(ions[i][0])
		kwargs['reading_func_args']['select_row'] = selr
		run_sal(vis[i], saved[i], vis_tf=False, path=path, n_rays=args.nrays, ds_file=args.ds_file, saved_add=f'data/', **kwargs)
		selr += 1

else:
	kwargs = dict(ray_dir=f'{path}rays', ion_list = list1)
	for i in range(nrows):
		#kwargs['ion_list'] = list(ions[i][0])
		#kwargs['reading_func_args']['select_row'] = selr
		run_sal(vis[i], saved[i], vis_tf=False, path=path, n_rays=args.nrays, ds_file=args.ds_file, saved_add=f'data/', **kwargs)
		#selr += 1



# selr = 0
# # nrays = 10
# for i in range(nrows):
# 	#kwargs['ion_list'] = list(ions[i][0])
# 	#kwargs['reading_func_args']['select_row'] = selr
# 	run_sal(vis[i], saved[i], vis_tf=False, ray_num=i, path=path, n_rays=nrays, saved_add=f'data/', **kwargs)
# 	#selr += 1
