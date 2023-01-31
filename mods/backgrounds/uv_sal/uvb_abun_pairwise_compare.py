import numpy as np
import pandas as pd

from sal_uvb_stats import gen_pairwise_data

def find_max_length(uvb_list):
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
    data set was shorter or longer than the first abundance data set.
    '''
    shorter = False
    
    if (start_1 < start_2 - clump_error) or (end_1 > end_2 + clump_error):
        shorter = True    
            
    return shorter
    

def check_split(end_big, small_ends, id_small, clump_error):
    '''
    Checks if the next few clumps in the second abundance data set all fall within
    the range of the first abundance data set
    '''
    clumps_within = 0

    # iterates through initial clump + number of clumps after this clump that 
    # are still within the bounds of clump 1

    while((id_small + clumps_within < len(small_ends)) and (small_ends[id_small + clumps_within] <= end_big - clump_error)):
        clumps_within +=1
    
    return clumps_within
    


def pairwise_compare(salsa_out, ion_list, nrays):
    '''
    Compares two different SALSA abudnance tables to one another created from two
    different Ultraviolet Backgrounds
    '''
    
    compare_dict = {}
    col_dens_1 = {}
    col_dens_2 = {}
    
    # should turn into argument
    clump_error = 10
    
    for ion in ion_list:
        
        ion_dat = salsa_out[ion]
        
        compare_dict[ion] = {}
        col_dens_1[ion] = {}
        col_dens_2[ion] = {}
        
        for ray in range(nrays):
            
            uvb_list = []
            
            for uvb_name in ion_dat.keys():
                
                uvb_dat = ion_dat[uvb_name]
                
                ray_dat = uvb_dat[uvb_dat["lightray_index"]==str(ray)].reset_index()
                
                uvb_list.append(ray_dat)
                
                # when including row data, add another for loop here
            
            uvb1 = uvb_list[0]
            uvb2 = uvb_list[1]
            
            match = {} ##create dictionaries to store indexes of clumps in the row that correspond to one another, keys will be row numbers and values will be indices except for the lonlies
            shorter = {}
            longer = {}
            split = {}
            merge = {}
            lonely_1 = []
            lonely_2 = []
            
            id1 = 0
            id2 = 0
            
            mx = find_max_length(uvb_list)
            print("mx",mx)
            while ((id1 < mx) and (id2 < mx)):

                if id1 >= len(uvb1["interval_start"])-1:

                    id1 = len(uvb1["interval_start"])-1
                
                if id2 >= len(uvb2["interval_start"])-1:

                    id2 = len(uvb2["interval_start"])-1 

                print("id1",id1)
                print("id2",id2)
                
                start_1 = uvb1["interval_start"][id1]
                start_2 = uvb2["interval_start"][id2]
                
                end_1 = uvb1["interval_end"][id1]
                end_2 = uvb2["interval_end"][id2]
                
                # checks to see if clumps are the same
                if (((start_2-clump_error) <= start_1 <= (start_2+clump_error)) and 
                     ((end_2 - clump_error) <= end_1 <= (end_2 + clump_error))):
                    match[id1] = id2
                    id1+=1
                    id2+=1
                    
                    continue
                
                # checks if clump2 is either shorter or longer than the other 
                elif (((start_2-clump_error) <= start_1 <= (start_2+clump_error)) or 
                     ((end_2 - clump_error) <= end_1 <= (end_2 + clump_error))):
                    
                    short_true = is_shorter(start_1, start_2, end_1, end_2, clump_error)
                    
                    # if clump 2 is shorter than clump 1, check if it is just a
                    # split up version of clump 1
                    if short_true:
                        
                        clumps_within = check_split(end_1, uvb2["interval_end"], id2, clump_error)
                        if clumps_within:
                            
                            split[id1] = []
                            
                            # puts current clump id and all of the clumps
                            # contained within clump 1 bounds into split dict
                            for clump in range(clumps_within+1):
                                
                                split[id1].append(id2+clump)

                            id1+=1
                            
                            # accounting for current clump and all clumps that are
                            # within clump 1
                            id2+= (1 + clumps_within)
                            
                            
                        else:
                            shorter[id1] = id2
                            id1+=1
                            id2+=1    
                            
                    # if clump 2 is longer than clump 1, check if clump 2 is 
                    # just a split up version of clump 1    
                    else:
                        
                        clumps_within = check_split(end_1, uvb2["interval_end"], id2, clump_error)
                        if clumps_within:
                            
                            merge[id2] = []
                            
                            # puts current clump id and all of the clumps
                            # contained within clump 1 bounds into split dict
                            for clump in range(clumps_within+1):
                                
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
                
                else:
                    
                    if (start_1 - end_2) > 0 :
                        
                        lonely_1.append(id1)
                        
                    else:
                        
                        lonely_2.append(id1)
                    
                    id1+=1
                    id2+=1
                       
            sorted_list = [match, shorter, longer, split, merge, lonely_1, lonely_2]
            
            compare_dict[ion][ray] = sorted_list
            
            col_dens_1[ion][ray], col_dens_2[ion][ray] = gen_pairwise_data(sorted_list, uvb_list, ion, ray)
    
    return compare_dict, col_dens_1, col_dens_2
    
