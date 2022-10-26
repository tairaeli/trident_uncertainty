import numpy as np
import pandas as pd


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

def is_shorter(start_1, start_2, end_1, end_2):
    '''
    Function for determining whether or not a given clump in the second abundance
    data set was shorter or longer than the first abundance data set.
    '''
    shorter = False
    
    if start_1 < start_2:
        shorter = True
        
    elif start_1 > start_2:
        shorter = False
        
    elif end_1 < end_2:
        shorter = False
    
    # the case where end_1 > end_2
    else: 
        shorter = True
            
    return shorter
    

def check_split(end_big, uvb_small, id_small):
    '''
    Checks if the next few clumps in the second abundance data set all fall within
    the range of the first abunfance data set
    '''
    clumps_within = 0
    
    # iterates through initial clump + number of clumps after this clump that 
    # are still within the bounds of clump 1
    while(uvb_small["interval_end"][id_small + clumps_within] <= end_big):
        clumps_within +=1
    
    return clumps_within
    
    
def pairwize_compare(salsa_out, ion_list, nrays):
    
    for ion in ion_list:
        
        ion_dat = salsa_out[ion]
        
        for ray in range(nrays):
            
            uvb_list = []
            
            for uvb_name in ion_dat.keys():
                
                uvb_dat = ion_dat[uvb_name]
                
                ray_dat = uvb_dat[uvb_dat["lightray_index"]==ray].reset_index()
                
                uvb_list.append(ray_dat)
                
                # when including row data, add another for loop here
            
            uvb1 = uvb_list[0]
            uvb2 = uvb_list[1]
            
            match = {} ##create dictionaries to store indexes of clumps in the row that correspond to one another, keys will be row numbers and values will be indecies except for the lonlies
            shorter = {}
            longer = {}
            split = {}
            merge = {}
            lonely_1 = []
            lonely_2 = []
            
            id1 = 0
            id2 = 0
            
            mx = find_max_length(uvb_list)
            
            while ((id1 < mx) and (id2 < mx)):
                
                start_1 = uvb1["interval_start"][id1]
                start_2 = uvb2["interval_start"][id2]
                
                end_1 = uvb1["interval_end"][id1]
                end_2 = uvb2["interval_end"][id2]
                
                # checks to see if clumps are the same
                if (start_1 == start_2) and (end_1 == end_2):
                    match[id1] = id2
                    id1+=1
                    id2+=1
                    
                    continue
                
                # checks if clump2 is either shorter or longer than the other 
                elif (start_1 == start2) or (end_1 == end_2):
                    
                    short_true = is_shorter(start_1, start_2, end_1, end_2)
                    
                    # if clump 2 is shorter than clump 1, check if it is just a
                    # split up version of clump 1
                    if short_true:
                        
                        clumps_within = check_split(end_1, uvb2, id2)
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
                        
                        clumps_within = check_split(end_1, uvb2, id2)
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
                            longer[id1] = id2
                            id1+=1
                            id2+=1
                
                else:
                    
                    if (start_1 - end_2) > 0 :
                        
                        lonely_1.append(id1)
                        
                    else:
                        
                        lonely_2.append(id1)

    
    return match, shorter, longer, split, merge, lonely_1, lonely_2
    
    
    
    
    
    
