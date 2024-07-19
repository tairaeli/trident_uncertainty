"""
Makes comparison between the SALSA output from 2 different Ultraviolet Background (UVB) models 
"""

import numpy as np
import pandas as pd
import pickle
import configparser
import argparse
import os

from condense_clumps import condense_pairwise_data

def find_max_length(uvb_list):
    """
    Finds which uvb data has the largest number of indices (gas clumps)

    args:

        uvb_list (List[Dataset]) - list containing datasets of uvb data
    
    returns:

        mx (int) - maximum number of gas clumps out of all of the elements in the uvb_list
    """
    mx= 0  ##find how long each array should be
    for ds in uvb_list: #find the cell index of the furthest clump
        if len(ds['interval_end']) == 0: ##handles if there are no clumps in a row
            break
        else:
            uvb_mx = len(ds["interval_end"])
            if uvb_mx>mx:
                mx=uvb_mx
    
    return mx

def is_shorter(start_1, start_2, end_1, end_2, clump_error):
    '''
    Function for determining whether or not a given clump in the second abundance
    data set is shorter than the first abundance data set.

    args:

        start_1 (int) - start position of the current clump 1

        start_2 (int) - start position of the current clump 2

        end_1 (int) - end position of current clump 1

        end_2 (int) - end position of current clump 2

        clump_error (int) - margin of error by which we consider corresponding clumps to be at equivalent indices
    
    returns:

        shorter (bool) - whether clump 2 is shorter than clump 1
    '''
    shorter = False
    
    if (start_1 <= start_2 + clump_error) and (end_1 >= end_2 - clump_error):
        shorter = True    
            
    return shorter

def is_longer(start_1, start_2, end_1, end_2, clump_error):
    '''
    Function for determining whether or not a given clump in the second abundance
    data set is longer than the first abundance data set.

    args:

        start_1 (int) - start position of the current clump 1

        start_2 (int) - start position of the current clump 2

        end_1 (int) - end position of current clump 1

        end_2 (int) - end position of current clump 2

        clump_error (int) - margin of error by which we consider corresponding clumps to be at equivalent indices
    
    returns:

        longer (bool) - whether clump 2 is longer than clump 1
    '''
    longer = False
    
    if (start_1 >= start_2 - clump_error) and (end_1 <= end_2 + clump_error):
        longer = True    
            
    return longer

def is_overlap(start_1, start_2, end_1, end_2, clump_error):
    '''
    Function for determining whether or not a given clump in the second abundance
    data set overlaps, but doesn't quite line up with the first abundace set.

    args:

        start_1 (int) - start position of the current clump 1

        start_2 (int) - start position of the current clump 2

        end_1 (int) - end position of current clump 1

        end_2 (int) - end position of current clump 2

        clump_error (int) - margin of error by which we consider corresponding clumps to be at equivalent indices
    
    returns:

        overlap (bool) - whether clump 2 overlaps with clump 1
    '''
    overlap = False
    
    if (end_2 >= start_1 - clump_error) and (start_2 <= start_1 + clump_error):
        overlap = True   

    elif (end_1 >= start_2 - clump_error) and (start_1 <= start_2 + clump_error):
        overlap = True    
            
    return overlap


def check_split(end_big, small_ends, id_small, clump_error):
    '''
    Checks if the next few clumps in the second abundance data set all fall within
    the range of the first abundance data set

    args:

        end_big (int) - end position of the clump that is assumed to be the larger of the two clumps

        small_ends (Series[int]) - Series of end positions of a series of clumps that may fit within the
                                   bounds of the supposed "larger" clump

        id_small - index along small_ends to begin

        clump_error (int) - margin of error by which we consider corresponding clumps to be at equivalent indices

    returns:

        clumps_within (int) - number of smaller clumps within the larger clump 
    '''
    clumps_within = 0
    # iterates through initial clump + number of clumps after this clump that 
    # are still within the bounds of clump 1

    while(((id_small + clumps_within) < len(small_ends)) and (small_ends[id_small + clumps_within] <= end_big - clump_error)):
        clumps_within +=1
    
    return clumps_within
    
def actually_lonely(id_comp, id_ques, uvb_comp, uvb_ques, clump_error):
    """
    Checks if the gas clump in question actually belongs to the lonely category

    args:

        id_comp (int) - index of the clump being compared to the clump in question

        id_ques (int) - index of the clump in question

        uvb_comp (Dataset) - dataset containing data to be compared to the clump in question

        uvb_ques (Dataset) - dataset containing the information on the clump in question

        clump_error (int) - margin of error by which we consider corresponding clumps to be at equivalent indices
    
    returns:

        True - if clump in question matches the comparison clump

        False - if the two clumps do not match
    """

    if (uvb_comp["interval_start"][id_comp] - clump_error) <= uvb_ques["interval_start"][id_ques] <= (uvb_comp["interval_start"][id_comp] + clump_error) or \
       (uvb_comp["interval_end"][id_comp] - clump_error) <= uvb_ques["interval_end"][id_ques] <= (uvb_comp["interval_end"][id_comp] + clump_error):
        
        return False
    else:
        return True

def pairwise_compare(salsa_out, nrays):
    '''
    Compares two different SALSA abudnance tables to one another created from two
    different Ultraviolet Backgrounds
    
    args:

        salsa_out (Dict) - output dictionary from SALSA run in 'sal_the_super_uvb.py'. 
                           Contains info on the gas clump identification process from SALSA for each UVB
        
        nrays (int) - represents the number of rays that are within the salsa_dict
        
    returns:
        
        compare_dict (Dict) - dictionary containing info on compared data for each ion for each ray
        
        col_dens1 (Dict) - contains 
    '''
    
    compare_dict = {}
    col_dens_1 = {}
    col_dens_2 = {}
    
    # should turn into argument
    clump_error = 10
    
    for ion in salsa_out.keys():
        
        ion_dat = salsa_out[ion]
        
        compare_dict[ion] = {}
        
        for ray in range(nrays):
            str_ray = None
            
            if len(str(ray)) != len(str(nrays)):
                n = len(str(nrays))
                
                str_ray = f'{ray:0{n}d}'
            else:
                str_ray = f'{ray}'

            uvb_list = []
            
            for uvb_name in ion_dat.keys():
                
                uvb_dat = ion_dat[uvb_name]

                ray_dat = uvb_dat[uvb_dat["lightray_index"]==str_ray].reset_index()
                
                uvb_list.append(ray_dat)
                
                # when including row data, add another for loop here
            
            uvb1 = uvb_list[0]
            uvb2 = uvb_list[1]
            
            match = {} ##create dictionaries to store indexes of clumps in the row that correspond to one another, keys will be row numbers and values will be indices except for the lonlies
            shorter = {}
            longer = {}
            overlap = {}
            split = {}
            merge = {}
            lonely_1 = []
            lonely_2 = []
            
            id1 = 0
            id2 = 0
            
            mx = find_max_length(uvb_list)
            # print("mx",mx)
            while ((id1 < mx) and (id2 < mx)):

                if id1 > len(uvb1["interval_start"])-1:

                    id1 = len(uvb1["interval_start"])-1
                
                if id2 > len(uvb2["interval_start"])-1:

                    id2 = len(uvb2["interval_start"])-1 

                # checks for instances where length of data is 0
                if id1 < 0:
                    print(f"Something weird with ray: {ray}, ion: {ion} in uvb 1")
                    break
                
                if id2 < 0:
                    print(f"Something weird with ray: {ray}, ion: {ion} in uvb 2")
                    break

                # print("id1",id1)
                # print("id2",id2)

                if len(lonely_1) != 0 and len(lonely_2) != 0:
                    # if a lonely clump was just added, check to see if it's actually lonely
                    if lonely_1[-1] == id1-1:

                        if not actually_lonely(id2, id1-1, uvb2, uvb1, clump_error):
                            lonely_1.remove(id1-1)
                            id1-=1
                    
                    if lonely_2[-1] == id2-1:
                        
                        if not actually_lonely(id1, id2-1, uvb1, uvb2, clump_error):
                            lonely_2.remove(id2-1)
                            id2-=1
     
                start_1 = uvb1["interval_start"][id1]
                start_2 = uvb2["interval_start"][id2]
                
                end_1 = uvb1["interval_end"][id1]
                end_2 = uvb2["interval_end"][id2]

                # checks to see if clumps are the same size
                if (((start_2-clump_error) <= start_1 <= (start_2+clump_error)) and 
                     ((end_2 - clump_error) <= end_1 <= (end_2 + clump_error))):
                    match[id1] = id2
                    id1+=1
                    id2+=1
                
                # checks if clump 2 shorter than clump 1
                elif is_shorter(start_1, start_2, end_1, end_2, clump_error):
                        
                    clumps_within = check_split(end_1, uvb2["interval_end"], id2, clump_error)
                    if clumps_within:
                        
                        split[id1] = []
                        
                        # puts current clump id and all of the clumps
                        # contained within clump 1 bounds into split dict
                        for clump in range(clumps_within):
                            
                            split[id1].append(id2+clump)

                        id1+=1
                        
                        # accounting for current clump and all clumps that are
                        # within clump 1
                        id2+= (1 + clumps_within)
                            
                    else:
                        shorter[id1] = id2
                        id1+=1
                        id2+=1    
                            
                # checks if clump 2 is longer than clump 1
                elif is_longer(start_1, start_2, end_1, end_2, clump_error):
                    
                    clumps_within = check_split(end_2, uvb1["interval_end"], id1, clump_error)
                    if clumps_within:
                        
                        merge[id2] = []
                        
                        # puts current clump id and all of the clumps
                        # contained within clump 1 bounds into split dict
                        for clump in range(clumps_within):
                            
                            merge[id2].append(id1+clump)
                        
                        # accounting for current clump and all clumps that are
                        # within clump 1
                        id1+= (1 + clumps_within)
                        
                        id2+=1
                    
                    else:
                        # clump 2 is longer than clump 1
                        longer[id1] = id2
                        id1+=1
                        id2+=1

                # checks for a rare overlapping case where a clump falls into neiter of
                # the above categories
                elif is_overlap(start_1, start_2, end_1, end_2, clump_error):
                    
                    overlap[id1] = id2
                    id1+=1
                    id2+=1 

                else:
                    # if the clumps do not fall into any other category, temorarily
                    # store them in lonely and evaluate in the next iteration
                    lonely_1.append(id1)
                    
                    lonely_2.append(id2)
                    
                    id1+=1
                    id2+=1
            
            # creating list of sorted catagories         
            sorted_list = [match, shorter, longer, split, merge, lonely_1, lonely_2, overlap]
            
            # storing list into output dict
            compare_dict[ion][ray] = sorted_list
    
    return compare_dict, col_dens_1, col_dens_2

def get_true_rs(val): ##define how to get actual rshift numbers
    if val == 20:
        true_rs = '2.0'
    elif val == 18:
        true_rs = '2.5'
    return true_rs

# reading in arguments
parser = argparse.ArgumentParser(description = "Select cutoff and UVB for analysis")

parser.add_argument('-uvb_path1', action='store', 
                    required=False, dest='uvb1', 
                    help='Path to UVB file')

parser.add_argument('-uvb_name1', action='store', 
                    required=False, dest='uvb_name1', 
                    help='Label to assign to uvb.')

parser.add_argument('-uvb_path2', action='store', 
                    required=False, dest='uvb2', 
                    help='Path to UVB file')

parser.add_argument('-uvb_name2', action='store', 
                    required=False, dest='uvb_name2', 
                    help='Label to assign to uvb.')

args = parser.parse_args()
dic_args = vars(args)

uvb_names = [args.uvb_name1, args.uvb_name2]
uvb_paths = [args.uvb1,args.uvb2]

sal_args = configparser.ConfigParser()
sal_args.read("./sal_params.par")

# set desired halo pattern
halo = sal_args["galaxy_settings"]["gal_pattern"]
rs = sal_args["galaxy_settings"]["redshift"]
nrays = int(sal_args["galaxy_settings"]["nrays"])

# gets the true rs needed
true_rs = get_true_rs(int(rs))

# identifying the path argument as a variable
out_file = sal_args["base_settings"]["output_file"]
path = os.path.expandvars(os.path.expanduser(out_file))
halo_path = path+'/halo'+f'{halo}'
rs_path = halo_path + '/redshift'+f'{true_rs}'

sal_dat = None

nuvb = 2
for i in range(nuvb):
    dat_path = rs_path +f'/{uvb_names[i]}/data'
    stat_path = rs_path +f'/{uvb_names[i]}/stats'

    if i == 0:
        with open(f'{dat_path}/salsa_out_dict.pickle',"rb") as dat:
            sal_dat = pickle.load(dat)
    else:
        with open(f'{dat_path}/salsa_out_dict.pickle',"rb") as dat:
            uvb_dict = pickle.load(dat)

            for ion in uvb_dict.keys():
        
                sal_dat[ion][uvb_names[i]] = uvb_dict[ion][uvb_names[i]]


print("Starting comparison", flush=True)
compare_dict, col_dens_1, col_dens_2 = pairwise_compare(sal_dat, nrays)

print("Starting condensing", flush=True)

ion_list = sal_args["galaxy_settings"]["ions"].split(" ")

uvb_dens_1 = {}
uvb_dens_2 = {}
dens_comp_dict = {}
for ion in ion_list:
    
    nion = ion.replace("_"," ")

    uvb_dens_1[nion] = {}
    uvb_dens_2[nion] = {}
    dens_comp_dict[nion] = {}
    for ray in range(nrays):
        
        uvb_dens_1[nion][ray], uvb_dens_2[nion][ray], dens_comp_dict[nion][ray] = condense_pairwise_data(compare_dict[nion][ray], 
                                                                                                         sal_dat[nion], 
                                                                                                         ray, 
                                                                                                         nrays)

out_dict = {uvb_names[0]:uvb_dens_1, uvb_names[1]:uvb_dens_2}

out_file = open(f'{rs_path}/uvb_clump_labels_{uvb_names[0]}_{uvb_names[1]}.pickle',"wb") ##saves the dictonaries so that they can be accesssed later
pickle.dump(compare_dict, out_file, protocol=3)	
out_file.close()

out_file = open(f'{rs_path}/uvb_compare_{uvb_names[0]}_{uvb_names[1]}.pickle',"wb") ##saves the dictonaries so that they can be accesssed later
pickle.dump(out_dict, out_file, protocol=3)	
out_file.close()

