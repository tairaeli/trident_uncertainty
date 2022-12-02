import salsa
from salsa.utils import check_rays
import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import yt  
from mpi4py import MPI
import pickle
from scipy import stats
import trident
import configparser

# set which analysis to perform
from abundance_compare import abundance_compare

# set statistics method
from saltistics import abundance_stats

comm = MPI.COMM_WORLD

print(f"Let's do some math, kids")

sal_args = configparser.ConfigParser()
sal_args.read("./sal_params.par")

def get_true_rs(val): ##define how to get actual rshift numbers
    if val == 20:
        true_rs = '2.0'
    elif val == 18:
        true_rs = '2.5'
    return true_rs

halo_names_dict = sal_args["halo_names"]

def get_halo_names(num):
    if str(num) in halo_names_dict.keys():
        halo_name = halo_names_dict[str(num)]
    return halo_name

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

# defining analysis parameters
# Note: these dictionaries are temporary and should most likely be included in the arguments at some point

def weighted_av(values, weights): ##define functions necessary in making statistics
    weighted_sum = []
    for value, weight in zip(values, weights):
        weighted_sum.append(value * weight)

    return sum(weighted_sum) / sum(weights)

def make_full_list(list_in, list_out):
    for element in list_in:
        list_out.append(element)
    return list_out

# indicate location of Foggie directory
foggie_dir = sal_args["base_settings"]["foggie_dir"]

# set desired halo pattern
halo = sal_args["galaxy_settings"]["gal_pattern"]
rs = sal_args["galaxy_settings"]["redshift"]

# gets the true rs needed
true_rs = get_true_rs(int(rs)) 

# takes in the foggie halo info directory
# outputs a dictionary of galactic center locations/velocities for all redshifts in each halo pattern
# NOTE: this function is temporary and has some hard-coded variables that will need to be changed
def foggie_defunker(foggie_dir):
    # initializing dictionary to store all of the galactic center data
    center_dat = {}

    # creating branch for each halo
    center_dat[halo] = {}
    # some hardcoded pipelies that will need to be changed
    cen_dat = pd.read_csv(f"{foggie_dir}/00{halo}/nref11c_nref9f/halo_c_v", sep = '|', names = ['null','redshift','name','xc','yc','zc','xv','yv','zv','null2'])
    # making some fixes specific to these files
    cen_dat = cen_dat.drop(0)
    cen_dat = cen_dat.drop(columns = ['null','null2'])
    # creating branch for each redshift in each halo 
    center_dat[halo][rs] = {}
    # isolating data to a specific redshift 
    rs_dat = cen_dat[cen_dat['name'] == ' RD00'+str(rs)+' ']
    # making 2 more branches to store the position and velocity data of the galactic center
    center_dat[halo][rs]['pos'] = [float(rs_dat["xc"]),float(rs_dat["yc"]),float(rs_dat["zc"])]
    center_dat[halo][rs]['vel'] = [float(rs_dat["xv"]),float(rs_dat["yv"]),float(rs_dat["zv"])]
     
    return center_dat

# fetching the galactic center data for all halo patterns and redshifts
center_dat = foggie_defunker(foggie_dir)

# identifying the path argument as a variable
out_file = sal_args["base_settings"]["output_file"]
path = os.path.expandvars(os.path.expanduser(out_file))

#  creating variable names for data bin locations
halo_path = path+'/halo'+f'{halo}'
rs_path = halo_path + '/redshift'+f'{true_rs}'
ray_path = rs_path +'/rays'
dat_path = rs_path +'/data'
stat_path = rs_path +'/stats'

# creating dictionaries to store all of our data (if they don't already exist)
if os.path.exists(halo_path) == False:
    os.mkdir(path+'/halo'+f'{halo}')
    
if os.path.exists(rs_path) == False:
    os.mkdir(rs_path)
    os.mkdir(ray_path) 
    os.mkdir(dat_path)
    os.mkdir(stat_path)
    
# load halo data
ds = yt.load(f'{sal_args["base_settings"]["halo_directory"]}/halo_00{halo}/nref11c_nref9f/RD00{rs}/RD00{rs}')

# define desired ions analyzed
ion_list = sal_args["galaxy_settings"]["ions"].split(" ")

# loading in list of different ionization tables for different UVBs
uvb_list = sal_args["uvb_analysis"]["uvb_files"].split(" ")

# define desired uvbs analyzed
uvb_names = sal_args["uvb_analysis"]["uvb_names"].split(" ")

# defining analysis parameters
# Note: these dictionaries are temporary and should most likely be included in the arguments at some point
center = ds.arr(center_dat[halo][rs]['pos'], 'kpc')
gal_vel = ds.arr(center_dat[halo][rs]['vel'], 'km/s')

# other_fields=['density', 'temperature', 'metallicity', ('index', 'radius')]
# max_impact=15 #kpc
# units = dict(density='g/cm**3', metallicity='Zsun')

# loading in data on other fields 
field_items = sal_args["galaxy_settings"]["field_items"].split(" ")
field_types = sal_args["galaxy_settings"]["field_types"].split(" ")
units = sal_args["galaxy_settings"]["field_units"].split(" ")

assert len(field_items) == len(field_types), f"Number of field items ({field_items}) do not match their labels ({field_types})"
assert len(field_items) == len(units), f"Number of field items ({field_items}) do not match their units ({units})"
other_fields = []
field_units = {}
# constructing field list to put into SALSA
for i, item in enumerate(field_items):
    
    if field_types[i] == "NA":
        field_tup = item
    else:
        field_tup = (field_types[i],item)
        
    other_fields.append(field_tup)
    
    if units[i] == "NA":
        continue
    else:
        field_units[f"{item}"] = units[i]

print(other_fields)

# defining some other parameters to include into SALSA
max_impact = int(sal_args["galaxy_settings"]["max_impact"])
nrays = int(sal_args["galaxy_settings"]["nrays"])
ray_num=f'{0:0{len(str(nrays))}d}'
ray_file=f'{ray_path}/ray{ray_num}.h5'

np.random.seed(13)

#get those rays babyyyy
# CK: Check that rays already exist, and that the have the additional fields contained
# in the third argument (empty for now; might become a user parameter)
check = check_rays(ray_path, nrays, [])
if not check:
    print("WARNING: rays not found. Generating new ones.")
    salsa.generate_lrays(ds, center.to('code_length'), nrays, max_impact, length=600, field_parameters={'bulk_velocity':gal_vel}, ion_list=['H I'], fields=other_fields, out_dir=ray_path)

ray_list=[]
for i in range(nrays):
    if len(str(i)) != len(str(nrays)):
        n = len(str(nrays)) - 1
        
        ray_list.append(f'{ray_path}/ray{i: 0{n}d}.h5')
    else:
        ray_list.append(f'{ray_path}/ray{i}.h5')

# CK: Taking a hint from SALSA on how to divvy up the ray list across procs
ray_arr = np.array(ray_list)
ray_files_split = np.array_split(ray_arr, comm.size)
my_rays = ray_files_split[ comm.rank ]

# add "_" in ion list for later for file writing
alt_ion_list = []
for i, ion in enumerate(ion_list):
    
    alt_ion = ion.replace("_"," ")
    alt_ion_list.append(alt_ion)
    
# set path to abundance table
abun_path = sal_args["base_settings"]["abundance_file"]

# if 'file_path' in dic_args:
abun = pd.read_csv(abun_path, delim_whitespace=True)
nrows = len(abun)

salsa_out_dict = {}
    
# for ion in alt_ion_list:
    
#     salsa_out_dict[ion] = {}
    
#     for iuvb, uvb in enumerate(uvb_list):

# do we need to check for this?
# if 'file_path' in dic_args:

abun_file = sal_args["base_settings"]["abundance_file"]

abun = pd.read_csv(abun_file, delim_whitespace=True)
nrows = len(abun)
saved = generate_names(nrows)
    
for ion in alt_ion_list:
    
    salsa_out_dict[ion] = {}
    
    for row_num in range(nrows):
    
        saved = generate_names(nrows)

        try:
            abundances = abun.iloc[row_num].to_dict()
            
            abs_ext = salsa.AbsorberExtractor(ds, ray_file, ion_name = ion, velocity_res =20, abundance_table = abundances, calc_missing=True)
            
            df = salsa.get_absorbers(abs_ext, my_rays, method='spice', fields=other_fields, units_dict=field_units).drop(columns='index')
            
            filename = f'{dat_path}/{saved[row_num]}_{ion.replace(" ", "_")}.txt'
            
            salsa_out_dict[ion][row_num] = df
            
            df.to_csv(filename, sep = ' ', index = False)
            
            print("Go look at your data!")
            
            
        except AttributeError: ##handles if there are no clumps in a halo
        
            df = pd.DataFrame(columns =['name', 'wave', 'redshift', 'col_dens', 'delta_v', 'vel_dispersion', 'interval_start', 'interval_end', 'density', 'temperature', 'metallicity', 'radius', 'lightray_index'], index = ['0'] )
            
            filename = f'{dat_path}/{saved[row_num]}_{ion.replace(" ", "_")}.txt'
            
            salsa_out_dict[ion][row_num] = df
            
            df.to_csv(f'{dat_path}/{saved[row_num]}_{i.replace(" ", "_")}_null.txt')
            
# dumping SALSA output to pickle file for debugging           
salsa_out_file = open(f'{dat_path}/salsa_out_dict.pickle',"wb") ##saves the dictonaries so that they can be accesssed later
pickle.dump(salsa_out_dict, salsa_out_file, protocol=3)	
salsa_out_file.close()

raynum = int(sal_args["galaxy_settings"]["nrays"])

compare_dict = abundance_compare(salsa_out_dict, alt_ion_list, raynum)

stat_dict = abundance_stats(compare_dict, salsa_out_dict, alt_ion_list, raynum, stat_path, halo, true_rs)


