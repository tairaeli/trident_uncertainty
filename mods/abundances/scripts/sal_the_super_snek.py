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
import pickle
from scipy import stats

comm = MPI.COMM_WORLD

print(f"Let's do some math, kids")

parser = argparse.ArgumentParser(description = "Preliminary constants for SALSA pipeline.")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path where rays and output data will be stored. Directory should contain three other directories called "data", "rays", and "visuals" for program to run smoothly.')
parser.add_argument('--nrays', action='store', dest='nrays', default=4, type=int, help='The number of rays to be generated.')
parser.add_argument('--abun', action='store', dest='file_path', default=argparse.SUPPRESS, help='Path to abundance file, if any. Defaults to solar abundances.')
parser.add_argument('--halo_dir', action='store', dest='halo_dir', default='/mnt/research/galaxies-REU/sims/FOGGIE', help='Path to halo data.')
parser.add_argument('--pat', action='store', dest='pattern', default=2392, type=int, help='Desired halo pattern file ID')
parser.add_argument('--rshift', action='store', dest='rs', default=20, type=int, help='Redshift file IDs')

args = parser.parse_args()
dic_args = vars(args)

def get_true_rs(val): ##define how to get actual rshift numbers
    if val == 20:
        true_rs = '2.0'
    elif val == 18:
        true_rs = '2.5'
    return true_rs

halo_names_dict = {'2392'  :  'Hurricane' ,'2878'  :  'Cyclone' , '4123'  :  'Blizzard' , '5016'  :  'Squall' ,'5036'  :  'Maelstrom' , '8508'  :  'Tempest'}

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

# EDIT THIS LINE TO LOCAL FOGGIE LOCATION

foggie_dir = "/mnt/home/f0104093/foggie/foggie/halo_infos"

# set desired halo pattern
halo = args.pattern
rs = args.rs
true_rs = get_true_rs(rs) ##gets the true rs needed

# takes in the foggie halo info directory
# outputs a dictionary of galactic center locations/velocities for all redshifts in each halo pattern
# NOTE: this function is temporary and has some hard-coded variables that will need to be changed
def foggie_defunker(foggie_dir):
    # initializing dictionary to store all of the galactic center data
    center_dat = {}

    # creating branch for each halo
    center_dat[halo] = {}
    # some hardcoded pipelies that will need to be changed
    cen_dat = pd.read_csv(f"/mnt/home/f0104093/foggie/foggie/halo_infos/00{halo}/nref11c_nref9f/halo_c_v", sep = '|', names = ['null','redshift','name','xc','yc','zc','xv','yv','zv','null2'])
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
path = os.path.expandvars(os.path.expanduser(args.path))

#  creating variable names for data bin locations
halo_path = path+'/halo'+f'{halo}'
rs_path = halo_path + '/redshift'+f'{true_rs}'
ray_path = rs_path +'/rays'
dat_path = rs_path +'/data'
stat_path = rs_path +'/stats'

# creating dictionaries to store all of our data (if they don't already exist)
if os.path.exists(halo_path) == "False":
    os.mkdir(path+'/halo'+f'{halo}')
    
if os.path.exists(rs_path) == "False":
    os.mkdir(ray_path) 
    os.mkdir(dat_path)
    os.mkdir(stat_path)
    

# load halo data
ds = yt.load(f'{args.halo_dir}/halo_00{halo}/nref11c_nref9f/RD00{rs}/RD00{rs}')
print(type(ds))
# defining analysis parameters
# Note: these dictionaries are temporary and should most likely be included in the arguments at some point
center = ds.arr(center_dat[halo][rs]['pos'], 'kpc')
gal_vel = ds.arr(center_dat[halo][rs]['vel'], 'km/s')
other_fields=['density', 'temperature', 'metallicity', ('index', 'radius')]
max_impact=15 #kpc
units_dict = dict(density='g/cm**3', metallicity='Zsun')

ray_num=f'{0:0{len(str(args.nrays))}d}'
ray_file=f'{ray_path}/ray{ray_num}.h5'

np.random.seed(13)

#get those rays babyyyy
# CK: Check that rays already exist, and that the have the additional fields contained
# in the third argument (empty for now; might become a user parameter)
check = check_rays(ray_path, args.nrays, [])
if not check:
    print("WARNING: rays not found. Generating new ones.")
    salsa.generate_lrays(ds, center.to('code_length'), args.nrays, max_impact, length=600, field_parameters={'bulk_velocity':gal_vel}, ion_list=['H I'], fields=other_fields, out_dir=ray_path)

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
new_ion_list = ['C_II', 'C_IV', 'O_VI']

# if 'file_path' in dic_args:
# 	abun = pd.read_csv(args.file_path, delim_whitespace=True)
# 	nrows = len(abun)
# 	saved = generate_names(nrows)
# 	for row_num in range(nrows):
# 		for i in ion_list:
# 			abundances = abun.iloc[row_num].to_dict()
# 			abs_ext = salsa.AbsorberExtractor(ds, ray_file, ion_name = i, velocity_res =20, abundance_table = abundances, calc_missing=True)
# 			df = salsa.get_absorbers(abs_ext, my_rays, method='spice', fields=other_fields, units_dict=units_dict).drop(columns='index')
# 			df.to_csv(f'{dat_path}/{saved[row_num]}_{i.replace(" ", "_")}.txt', sep = ' ')
# 			print("Go look at your data!")

# else:
# 	nrows = 0
# 	saved = generate_names(nrows)
# 	for i in ion_list:
# 		abs_ext = salsa.AbsorberExtractor(ds, ray_file, ion_name = i, abundance_table = None, calc_missing=True)
# 		df = salsa.get_absorbers(abs_ext, my_rays, method='spice', fields=other_fields, units_dict=units_dict)
# 		df.to_csv(f'{dat_path}/data_SolAb_{i.replace(" ", "_")}.txt', sep = ' ').drop(columns='index')
# 		print("Go look at your data!")
            
            
for ion in new_ion_list:
    ##categorize clumps##
    datanum = len(pd.read_csv(args.file_path, delim_whitespace=True)) ##number of rows on the abundance table
    ndigits= len(str(datanum))
    raynum = args.nrays ##number of rays used, modify as needed

    for r in range(raynum):
        rowlist = []
        for i in range(datanum):
            m = i+1
            n_len = len(str(m))
            n_zeros = ndigits - n_len
            k = "0" * n_zeros + str(m)
            row_data = pd.read_csv(dat_path+f"/data_AbundanceRow{k}_{ion}.txt", delim_whitespace=True) ##read in data files
            row_work = row_data[row_data["lightray_index"]==r] ##filter to only ray1
            df = row_work.reset_index().drop(columns="index") ##make indexing work
            rowlist.append(df)

        mx= 0  ##find how long each array should be
        for ds in rowlist: #find the cell index of the furthest clump
            if len(ds['interval_end']) == 0: ##handles if there are no clumps in a row
                break
            else:
                row_mx = max(ds["interval_end"])
                if row_mx>mx:
                    mx=row_mx
        
    
        super_clumps = np.zeros(int(mx))
        clmaps = []

        problems = []
        hassles = {}

        row_tracker = 0
        for ds in rowlist: ##make masks for each row and form super_clumps
            ds_clump_loc = np.zeros(int(mx))
            row_tracker +=1
            hassles_ind = []
            for j in range(ds.shape[0]):
                ds_clump_loc[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1
                    #print(f"Row Start: {int(ds['interval_start'][j])}, Row End: {int(ds['interval_end'][j])}")

                if j-1 == -1:
                    super_clumps[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1

                elif int(ds["interval_end"][j-1]) == int(ds["interval_start"][j]):

                        #check for wether or not two clumps are within thier standard deviations
                    if not (float(ds["delta_v"][j-1]) - 1.5 * (float(ds["vel_dispersion"][j-1]))) <= float(ds["delta_v"][j]) <= 1.5 * (float(ds["delta_v"][j-1]) + (float(ds["vel_dispersion"][j-1]))):
                        ds_clump_loc[int(ds["interval_start"][j])] = 2
                    #edge case handling
                        if int(ds["interval_end"][j-1]) not in problems:
                            problems.append(int(ds["interval_end"][j-1])) 
                             #to make sure the bigger clumps stays in superclumps:
                        if (super_clumps[int(ds["interval_end"][j-1])] == 0) or (super_clumps[int(ds["interval_end"][j-1])] == 2):
                            super_clumps[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1
                            super_clumps[int(ds["interval_start"][j])] = 2
                    else:
                        hassles_ind.append([int(ds["interval_start"][j-1]), int(ds["interval_end"][j-1])])
                        hassles_ind.append([int(ds["interval_start"][j]), int(ds["interval_end"][j])])                
                        super_clumps[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1
                        print(f"Hassles: {row_tracker, hassles_ind}")

                else:
                    super_clumps[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1  

            hassles[row_tracker]=hassles_ind
            clmaps.append(ds_clump_loc)

            for index in problems:
                for row in clmaps:
                    if index in row:
                        if ((row[index-1] == 0 and row[index+1] != 0) or (row[index+1] == 0 and row[index-1] != 0)) and (row[index] != 0): ##edge case handling
                            super_clumps[index] = 2

        super_clumps=np.append(0, super_clumps)  ##make indexing work  
        super_clumps=np.append(super_clumps, 0)
        np.save(f'super_clumps_array_{ion}_ray{r}', super_clumps) ##save super_clumps for future reference
        match = {} ##create dictionaries to store indexes of clumps in the row that correspond to one another, keys will be row numbers and values will be indecies except for the lonlies
        short = {}
        split = {}
        lonely = {}
        maybe_lonely = {}
        rownum = 0

        for row in clmaps:  ##start by iterating over a whole row
            row = np.append(0,row) # adding an extra element to prevent booleans from failing
            row = np.append(row,0)
            rownum += 1 ##define which row we're working on
            ##make all the variables and lists necessary
            row_st_cnt = 0 ##count how many starts of row clumps there have been within a super clump
            row_st_ind = []  ##keep track of start location(s) of a row clump

            row_en_cnt = 0 ##how many ends of row clumps there have been within a super clump
            row_en_ind = [] ##keep track of end location(s) of a row clump
            row_match = []
            row_short = []
            row_split = []
            sup_st_true =[]
            weird_index = 0

            for i in range(len(row)):

                if super_clumps[i-1]<super_clumps[i]: ##start of a super clump
                    sup_st = i-1 ##keep track of start location of a super clump
                    sup_st_true.append(sup_st)

                if row[i-1]<row[i] and (len(row_st_ind) == 0 or len(row_st_ind) ==1):##start of a clump in the row
                    if super_clumps[i-1] == 2: ##if oops is printed, there is another edge case, debugging must begin again
                        print("oops")
                    elif super_clumps[i] != 2 or len(row_st_ind)!=0: ##normal clumps look like this
                        row_st_ind.append(i-1)
                        row_st_cnt += 1
                        print(row_st_ind, rownum)
                    else: ##edge case handling ##FIXME
                        weird_index = i-1

                elif row[i-1]>row[i]: ##end of a clump in row
                    row_en_cnt += 1
                    if row[i-1]==2: ## edge case handling
                        row_en_ind.append(i-2)
                    else:
                        row_en_ind.append(i-1)

                if super_clumps[i-1]>super_clumps[i]: ##end of a super clump
                    if super_clumps[i-1]==2:
                        sup_en = i-2
                        sup_st = sup_st_true[0] ##edge case handling

                    if super_clumps[i-1] == 0 or super_clumps[i-1] == 1:
                        sup_en = i-1 ##keep track of the location of the end of a super clump



                    if (row_st_cnt == 1) or (row_en_cnt == 1): ##check for if there is only one row clump in the super clump
                        if weird_index != 0: ##edge case handling
                            row_st_ind.append(weird_index)

                        if (len(row_en_ind) == 0 and len(row_st_ind) != 0) or (len(row_st_ind) == 0 and len(row_en_ind) != 0):

                            break

                        hassle_st = []
                        hassle_en = []

                        if rownum in hassles.keys(): ##even more edge case handling
                            for value in hassles.values():
                                for t in value:
                                    hassle_st.append(t[0])
                                    hassle_en.append(t[1])

                            if row_st_ind[0] in hassle_st:
                                for j in range(len(hassle_st)):
                                    if [hassle_st[j], hassle_en[j]] not in row_split:
                                        row_split.append([hassle_st[j], hassle_en[j]])
                            elif row_en_ind[0] in hassle_en:

                                continue

                        if (row_st_ind[0] == sup_st) and (row_en_ind[0] == sup_en) and (row_st_ind[0] not in hassle_st) and (row_en_ind[0] not in hassle_en): ##if the starts and ends match, the clumps are identical

                            row_match.append([row_st_ind[0],row_en_ind[0]]) ##thus, start and end indecies appended to a list of them

                        elif (row_st_ind[0] not in hassle_st) and (row_en_ind[0] not in hassle_en):
                            row_short.append([row_st_ind[0],row_en_ind[0]]) ##if not, then the row clump must be shorter and the start and end indecies are appended to the appropraite list

                    elif (row_st_cnt == 0) and (row_en_cnt == 0): ##check if there is nothing in the row that matches the super clump

                        if str([sup_st,sup_en]) in maybe_lonely.keys(): ##check if we have already seen this superclump, if not make the entry in the dictionary
                            maybe_lonely[str([sup_st,sup_en])] += 1
                        else:
                            maybe_lonely[str([sup_st,sup_en])] = 1

                    else: ##only other senario is there there was a split
                        for j in range(min([len(row_en_ind), len(row_st_ind)])): ##organize the indecies to make the list in order
                            row_split.append([row_st_ind[j],row_en_ind[j]]) 


                    

                    ##reset variables      
                    row_st_cnt = 0
                    row_st_ind = []
                    row_en_cnt = 0
                    row_en_ind = []
                    sup_st_true = []
                    match[rownum] = row_match
                    short[rownum] = row_short
                    split[rownum] = row_split

                    if super_clumps[i-1] == 2:
                        sup_st = i-2
                        if row[i-1]>row[i]:
                            row_st_ind.append(i-2)
                            row_st_cnt += 1

    #for clump in maybe_lonely: ##later must set a limits on the minimumm number of times something has to appear in maybe lonely for it to actually be considered lonely

        pickling_match = open(f"Match_{ion}_Ray{r}.pickle","wb") ##saves the dictonaries so that they can be accesssed later
        pickle.dump(match, pickling_match, protocol=3)	
        pickling_match.close()

        pickling_split = open(f"Split_{ion}_Ray{r}.pickle","wb")
        pickle.dump(split, pickling_split, protocol=3)
        pickling_split.close() 

        pickling_short = open(f"Short_{ion}_Ray{r}.pickle","wb")
        pickle.dump(short, pickling_short, protocol=3)
        pickling_short.close()

        pickling_maybe_lonely = open(f"MaybeLonely_{ion}_Ray{r}.pickle","wb")
        pickle.dump(maybe_lonely, pickling_maybe_lonely, protocol=3)
        pickling_maybe_lonely.close()
    
    ####make stats files for each super_clump#####
    
    ray_nums = [] ##define lists that will eventually go into pandas sets
    super_cl_nums = []
    med_col_dens = []
    mad_for_med = []
    distances = []
    central_v = []
    vel_dispersions =[]
    densities = []
    temperatures = []
    num_clumps = []
    freq_split_short = []
    rows_of_rep_clumps = []
    cat_rep_clump = []
    diff_from_sol = []
    sol_ab_col_dens = 0

    
    for r in range(raynum):
        if os.path.exists(f'/mnt/scratch/f0104093/cgm_abundance_variance/Match_{ion}_Ray{r}.pickle'):
            pickle_match_off = open(f"Match_{ion}_Ray{r}.pickle", 'rb') ##get previously made data
            match = pickle.load(pickle_match_off)
    
            pickle_split_off = open(f"Split_{ion}_Ray{r}.pickle", 'rb')
            split = pickle.load(pickle_split_off)
    
            pickle_short_off = open(f"Short_{ion}_Ray{r}.pickle", 'rb')
            short = pickle.load(pickle_short_off)
    
            super_clumps = np.load(f'super_clumps_array_{ion}_ray{r}.npy')
    
            var_rows = []
            
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
        clump_stats["diff_from_solar_abun"] = diff_from_sol
        clump_stats["distance_from_galaxy"] = distances 
        clump_stats["central_velocity"] = central_v 
        clump_stats["vel_dispersion"] = vel_dispersions 
        clump_stats["density"] = densities 
        clump_stats["temperature"] = temperatures
        clump_stats["num_split_or_short"] = freq_split_short
        clump_stats["num_of_clumps"] = num_clumps
        clump_stats["rep_clump_row"] = rows_of_rep_clumps
        clump_stats["category_rep_clump"] = cat_rep_clump
    
        print(len(ray_nums), len(super_cl_nums), len(med_col_dens), len(mad_for_med), len(central_v), len(vel_dispersions))
        
        df = pd.DataFrame.from_dict(clump_stats)
        df.to_csv(f"{stat_path}/{halo}_z{true_rs}_{ion}_abun_all-model-families_all-clumps.csv" ,sep = ' ', na_rep = 'NaN') ##save the files to scratch
        ##as we're done with each file, delete it so we don't get residual data we don't need
   
        os.remove(f"Match_{ion}_Ray{r}.pickle")
        os.remove(f"Short_{ion}_Ray{r}.pickle")
        os.remove(f"Split_{ion}_Ray{r}.pickle")
        os.remove(f'super_clumps_array_{ion}_ray{r}.npy')
