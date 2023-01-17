import numpy as np
import pandas as pd

def gen_pairwise_data(sorted_list, uvb_list, ion, ray):

    
    match = sorted_list[0]
    shorter = sorted_list[1]
    longer = sorted_list[2]
    split = sorted_list[3]
    merge = sorted_list[4]
    lonely_1 = sorted_list[5]
    lonely_2 = sorted_list[6]
    
    uvb1 = uvb_list[0]
    uvb2 = uvb_list[0]
    
    id1 = 0
    id2 = 0    
    
    uvb1_col_dens = []
    uvb2_col_dens = []
    
    while (id1 < len(uvb1["col_dens"]) & id2 < len(uvb2["col_dens"]))
        
        # checking for instances where uvb2 has split clumps relative to uvb1
        if id1 in split:
            
            uvb1_col_dens.append(uvb1["col_dens"][id1])
            
            clump_piece_ids = split[id1]
            
            num_pieces = len(clump_piece_ids)
            
            uvb2_col_dens_avg = sum(uvb2["col_dens"][clump_piece_ids])/num_pieces
            
            uvb2_col_dens.append(uvb2_col_dens_avg)
            
            id2 += (num_pieces - 1)
        
        # checking for instances where uvb2 has merged clumps relative to uvb1
        elif id2 in merge:
            
            uvb2_col_dens.append(uvb2["col_dens"][id2])
            
            clump_piece_ids = split[id2]
            
            num_pieces = len(clump_piece_ids)
            
            uvb1_col_dens_avg = sum(uvb1["col_dens"][clump_piece_ids])/len(clump_piece_ids)
            
            uvb1_col_dens.append(uvb1_col_dens_avg)
            
            id1 += (num_pieces - 1)
        
        elif id1 in lonely_1:
            
            uvb1_col_dens.append(uvb1["col_dens"][id1])
            
            uvb2_col_dens.append(np.nan)
    
        elif id2 in lonely_2:
            
            uvb2_col_dens.append(uvb2["col_dens"][id2])
            
            uvb2_col_dens.append(np.nan)
        
        else:
            
            uvb1_col_dens.append(uvb1["col_dens"][id1])
            
            uvb2_col_dens.append(uvb2["col_dens"][id2])
            
            
        id1+=1
        id2+=1