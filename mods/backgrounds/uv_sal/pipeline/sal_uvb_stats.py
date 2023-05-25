import numpy as np
import pandas as pd

def gen_pairwise_data(sorted_list, uvb_list, ion, ray):
    
    # unpacking data
    match = sorted_list[0]
    shorter = sorted_list[1]
    longer = sorted_list[2]
    split = sorted_list[3]
    merge = sorted_list[4]
    lonely_1 = sorted_list[5]
    lonely_2 = sorted_list[6]
    
    uvb1 = uvb_list[0]
    uvb2 = uvb_list[1]
    
    id1 = 0
    id2 = 0    
    
    # initializing output data storage bins
    uvb1_col_dens = []
    uvb2_col_dens = []

    while (id1 < len(uvb1["col_dens"])) and (id2 < len(uvb2["col_dens"])):
        
        # checking for instances where uvb2 has split clumps relative to uvb1
        if id1 in split:

            uvb1_col_dens.append(uvb1["col_dens"][id1])
            
            clump_piece_ids = split[id1]
            
            uvb2_col_dens_avg = np.mean(uvb2["col_dens"].iloc[clump_piece_ids])
            
            uvb2_col_dens.append(uvb2_col_dens_avg)
            
            id2 += (len(clump_piece_ids) - 1)
        
        # checking for instances where uvb2 has merged clumps relative to uvb1
        elif id2 in merge:
            
            uvb2_col_dens.append(uvb2["col_dens"][id2])
            
            clump_piece_ids = merge[id2]
            
            uvb1_col_dens_avg = np.mean(uvb1["col_dens"].iloc[clump_piece_ids])
            
            uvb1_col_dens.append(uvb1_col_dens_avg)
            
            id1 += (len(clump_piece_ids) - 1)
        
        elif id1 in lonely_1:
            
            uvb1_col_dens.append(uvb1["col_dens"][id1])
            
            uvb2_col_dens.append(0)
    
        elif id2 in lonely_2:
            
            uvb2_col_dens.append(uvb2["col_dens"][id2])
            
            uvb1_col_dens.append(0)
        
        # if none of the above cases pass, then the clump indices line up
        else:
            
            uvb1_col_dens.append(uvb1["col_dens"][id1])
            
            uvb2_col_dens.append(uvb2["col_dens"][id2])
            
            
        id1+=1
        id2+=1
    
    # If everything worked correctly, the two output arrays should be the same length
    assert len(uvb1_col_dens) == len(uvb2_col_dens), f"column density arrays are different sizes. \n uvb1 = {len(uvb1_col_dens)} \n uvb2 = {len(uvb2_col_dens)}"

    return uvb1_col_dens, uvb2_col_dens
    