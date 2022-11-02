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

from uvb_abun_pairwize_compare import pairwise_compare

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

# set to 0 for now, but could easily buid for loop to include multiple rows
# may need to expand the salsa_out_dict
row_num = 0

salsa_out_dict = {}
    
for ion in alt_ion_list:
    
    salsa_out_dict[ion] = {}
    
    for iuvb, uvb in enumerate(uvb_list):
        # impliments the ionization table for each different UVB model
        trident.ion_balance.add_ion_fields(ds, ions = alt_ion_list, ionization_table = uvb)
        saved = generate_names(nrows,uvb_names[iuvb])
        
    
        try:
            abundances = abun.iloc[row_num].to_dict()
            
            abs_ext = salsa.AbsorberExtractor(ds, ray_file, ion_name = ion, velocity_res =20, abundance_table = abundances, calc_missing=True)
            
            df = salsa.get_absorbers(abs_ext, my_rays, method='spice', fields=other_fields, units_dict=field_units).drop(columns='index')
            
            filename = f'{dat_path}/{saved[row_num]}_{ion.replace(" ", "_")}.txt'
            
            salsa_out_dict[ion][uvb_names[iuvb]] = df
            
            df.to_csv(filename, sep = ' ', index = False)
            
            print("Go look at your data!")
            
            
        except AttributeError: ##handles if there are no clumps in a halo
        
            df = pd.DataFrame(columns =['name', 'wave', 'redshift', 'col_dens', 'delta_v', 'vel_dispersion', 'interval_start', 'interval_end', 'density', 'temperature', 'metallicity', 'radius', 'lightray_index'], index = ['0'] )
            
            filename = f'{dat_path}/{saved[row_num]}_{ion.replace(" ", "_")}.txt'
            
            salsa_out_dict[ion][uvb_names[iuvb]] = df
            
            df.to_csv(f'{dat_path}/{saved[row_num]}_{i.replace(" ", "_")}_null.txt')
            
            
pickling_match = open(f'{dat_path}/salsa_out_dict.pickle',"wb") ##saves the dictonaries so that they can be accesssed later
pickle.dump(salsa_out_dict, pickling_match, protocol=3)	
pickling_match.close()

compare_dict = pairwise_compare(salsa_out_dict, ion_list, raynum)


    for r in range(raynum):
            
            for i in range(datanum):
                m = i+1
                n_len = len(str(m))
                n_zeros = ndigits - n_len
                p = "0" * n_zeros + str(m)
                row_data = pd.read_csv(dat_path+f"/data_AbundanceRow{p}_{ion}.txt", delim_whitespace=True) ##read in data files
                row_work = row_data[row_data["lightray_index"]==r] ##filter to only one ray
                df = row_work.reset_index().drop(columns="index") ##make indexing work
                var_rows.append(df)
    
            sup_st = [] ##get indexes of superclumps
            sup_en = []
            for i in range(1, len(super_clumps)):
                n= i-1
                if super_clumps[n]<super_clumps[i]: ##start of a super clump
                    if super_clumps[n]== 2:
                        sup_st.append(n-1)
                    else: 
                        sup_st.append(n)
                elif super_clumps[n]>super_clumps[i]: ##end of a super clump
                    if super_clumps[n] ==2:
                        sup_en.append(n-1)
                    else:
                        sup_en.append(n)
            num_spl_sho = 0
    
            for k in range(len(sup_st)): ##depending on which category each clump belongs to in super_clumps
                col_density_match = [] ##lists for col_density
                col_density_split = []
                col_density_short = []
                match_done = [] ##keep track of which super clumps already have a representative
                short_done = []
                split_done = []
    
                for row, index in match.items(): ##first, see what matches the super clump
    
                    for j in range(len(index)):
    
                        if (index[j][0]>=sup_st[k]) and (index[j][1]<=sup_en[k]):
    
                            ds = var_rows[row-1]
                            indexq = np.where((index[j][0]) == (ds["interval_start"]))
                            if row != 1: ##exclude solar abundance
                                col_density_match.append(ds["col_dens"][int(indexq[0])]) ##get column density
                            elif row == 1: ##keep frack of what the solar abundance is 
                                sol_ab_col_dens = ds["col_dens"][int(indexq[0])]
                            if sup_st[k] not in match_done: ##if this is the first one done for a super clump, get all the other data and make this clump a "representative" of the super clump
                                print('Match:', row, index[j][0], sup_st[k])
                                distances.append(ds["radius"][int(indexq[0])])
                                central_v.append(ds["delta_v"][int(indexq[0])])
                                vel_dispersions.append(ds["vel_dispersion"][int(indexq[0])])
                                densities.append(ds["density"][int(indexq[0])])
                                temperatures.append(ds["temperature"][int(indexq[0])])
                                rows_of_rep_clumps.append(row)
                                cat_rep_clump.append('match')
                                match_done.append(sup_st[k])
    
    
                for rows, indexs in short.items(): ##exactly the same method as match, but do it second bc if there is a match, we want to use that as our representataive
    
                    for j in range(len(indexs)):
    
                        if (indexs[j][0]>=sup_st[k]) and (indexs[j][1]<=sup_en[k]):
                            num_spl_sho += 1
                            ds = var_rows[rows-1]
                            indexq = np.where((indexs[j][0]) == (var_rows[rows-1]["interval_start"]))
                            if rows != 1:
                                col_density_short.append(ds["col_dens"][int(indexq[0])]) ##get column density
                            elif rows == 1:
                                sol_ab_col_dens = ds["col_dens"][int(indexq[0])]
    
                            if (len(match_done) == 0) and (sup_st[k] not in short_done):
                                print('Short:', row, indexs[j][0], sup_st[k])
                                distances.append(ds["radius"][int(indexq[0])])
                                central_v.append(ds["delta_v"][int(indexq[0])])
                                vel_dispersions.append(ds["vel_dispersion"][int(indexq[0])])
                                densities.append(ds["density"][int(indexq[0])])
                                temperatures.append(ds["temperature"][int(indexq[0])])
                                rows_of_rep_clumps.append(row)
                                cat_rep_clump.append('short')
                                short_done.append(sup_st[k])
    
                for rowm, indexm in split.items(): ##split is a bit weird bc there are multiple in a single super clump
                    temp_col_dens =[] ##create temporary lists 
                    temp_delta_v = []
                    temp_vel_dis = []
                    temp_dens = []
                    temp_temp = []
                    temp_rad = []
                    col_dens_for_weights = []
    
                    for j in range(len(indexm)):
    
                        if (indexm[j][0]>=sup_st[k]) and (indexm[j][1]<=sup_en[k]):
                            ds = var_rows[rowm-1]
                            indexq = np.where((indexm[j][0]) == (var_rows[rowm-1]["interval_start"]))
                                
                            if len(list(indexq[0])) > 0:
                                temp_col_dens.append(10 ** ds["col_dens"][int(indexq[0])]) ##col_dens is on a log scale, so to add it all togeher must make it the exponent of 10
                                col_dens_for_weights.append(ds["col_dens"][int(indexq[0])]) 
                                temp_delta_v.append(ds["delta_v"][int(indexq[0])])## append the other dictionaries with the temporary values
                                temp_vel_dis.append(ds["vel_dispersion"][int(indexq[0])])
                                temp_dens.append(ds["density"][int(indexq[0])])
                                temp_temp.append(ds["temperature"][int(indexq[0])])
                                temp_rad.append(ds["radius"][int(indexq[0])])
    
                                if (len(match_done) == 0) and (len(short_done) == 0) and (sup_st[k] not in split_done):
                                    print('Split:', row, indexm[j][0], sup_st[k])
                                    distances.append(weighted_av(temp_rad, col_dens_for_weights)) ##use the col_dens as a weight for each value so we get a weight average if a split is chosen as a "represenatative"
                                    central_v.append(weighted_av(temp_delta_v, col_dens_for_weights))
                                    vel_dispersions.append(weighted_av(temp_vel_dis, col_dens_for_weights))
                                    densities.append(weighted_av(temp_dens, col_dens_for_weights))
                                    temperatures.append(weighted_av(temp_temp, col_dens_for_weights))
                                    rows_of_rep_clumps.append(row)
                                    split_done.append(sup_st[k])
                                    cat_rep_clump.append('split')
    
    
                    if len(temp_col_dens) != 0 and rowm != 1: ##finally, get one value for the col_dens of the whole thing
                        num_spl_sho +=1
                        log_sum_dens = np.log10(sum(temp_col_dens))
                        col_density_split.append(log_sum_dens)
                    elif len(temp_col_dens) != 0 and rowm == 1:
                        num_spl_sho +=1
                        sol_ab_col_dens = np.log10(sum(temp_col_dens))
    
              ##make lists to put into dictionaries
                ray_nums.append(r)
                super_cl_nums.append(k)
    
                full_col_density = []
                make_full_list(col_density_match, full_col_density)
                make_full_list(col_density_split, full_col_density)
                make_full_list(col_density_short, full_col_density)
                med_col_dens.append(np.median(full_col_density))
                mad_for_med.append(stats.median_abs_deviation(full_col_density))
                
                num_clumps.append(len(full_col_density))
                freq_split_short.append(num_spl_sho)
                num_spl_sho = 0
                
                if sol_ab_col_dens != 0: ##handle the case where there is no col_density for the solar abundance
                    diff_from_sol.append(np.log10((10 ** sol_ab_col_dens) -(10 ** np.median(full_col_density))))
                else:
                    diff_from_sol.append(np.NaN)

        ##make and fill the dictionary that will be convered into a csv file
        clump_stats = {}
        clump_stats["ray_num"] = ray_nums
        clump_stats["super_clump_number"] = super_cl_nums 
        clump_stats["median_col_desnity"] = med_col_dens 
        clump_stats["mad_for_col_desnity"] = mad_for_med
        # clump_stats["diff_from_solar_abun"] = diff_from_sol
        # clump_stats["distance_from_galaxy"] = distances 
        # clump_stats["central_velocity"] = central_v 
        # clump_stats["vel_dispersion"] = vel_dispersions 
        # clump_stats["density"] = densities 
        # clump_stats["temperature"] = temperatures
        # clump_stats["num_split_or_short"] = freq_split_short
        # clump_stats["num_of_clumps"] = num_clumps
        # clump_stats["rep_clump_row"] = rows_of_rep_clumps
        # clump_stats["category_rep_clump"] = cat_rep_clump
    
        print(len(ray_nums), len(super_cl_nums), len(med_col_dens), len(mad_for_med), len(central_v), len(vel_dispersions))
        
        df = pd.DataFrame.from_dict(clump_stats)
        df.to_csv(f"{stat_path}/{halo}_z{true_rs}_{ion}_abun_all-model-families_all-clumps.csv" ,sep = ' ', na_rep = 'NaN') ##save the files to scratch

        ##as we're done with each file, delete it so we don't get residual data we don't need
        os.remove(f"{stat_path}/Match_{ion}_Ray{r}.pickle")
        os.remove(f"{stat_path}/Short_{ion}_Ray{r}.pickle")
        os.remove(f"{stat_path}/Split_{ion}_Ray{r}.pickle")
        os.remove(f'{stat_path}/super_clumps_array_{ion}_ray{r}.npy')
