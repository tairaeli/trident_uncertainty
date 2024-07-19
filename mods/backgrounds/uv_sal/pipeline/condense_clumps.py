import numpy as np
import pandas as pd

def condense_pairwise_data(sorted_list, ion_dat, ray, nrays):
    """
    Takes in data from a pairwise comparison between the clump
    data from the two different UVB models

    args:

        sorted_list (List) - list containing data from the sorted
                             SALSA gas clumps for a given ion and 
                             ray
        
        ion_dat (Dict) - contains information on a given ion. Should 
                         be from SALSA output
        
    returns:

        uvb1_col_dens (Dict) - column densities of the newly condensed
                               clump data for whichever data is associated
                               with UVB model 1

        uvb2_col_dens (Dict) - column densities of the newly condensed
                               clump data for whichever data is associated
                               with UVB model 2
    """
    
    # unpacking data
    match = sorted_list[0]
    shorter = sorted_list[1]
    longer = sorted_list[2]
    split = sorted_list[3]
    merge = sorted_list[4]
    lonely_1 = sorted_list[5]
    lonely_2 = sorted_list[6]
    overlap = sorted_list[7]

    # 
    nsplit = {}
    nmerge = {}

    uvb_list = []
    
    for uvb_name in ion_dat.keys():
        
        str_ray = None
        if len(str(ray)) != len(str(nrays)):
            n = len(str(nrays))
            
            str_ray = f'{ray:0{n}d}'
        else:
            str_ray = f'{ray}'
        
        uvb_dat = ion_dat[uvb_name]

        ray_dat = uvb_dat[uvb_dat["lightray_index"]==str_ray].reset_index()
        
        uvb_list.append(ray_dat)

    uvb1 = uvb_list[0]
    uvb2 = uvb_list[1]
    id1 = 0
    id2 = 0    
    
    # initializing output data storage bins

    uvb1_out = {}
    uvb1_out["col_dens"] = []
    uvb1_out["density"] = []
    uvb1_out["temperature"] = []
    uvb1_out["metallicity"] = []

    uvb2_out = {}
    uvb2_out["col_dens"] = []
    uvb2_out["density"] = []
    uvb2_out["temperature"] = []
    uvb2_out["metallicity"] = []

    # uvb1_col_dens = []
    # uvb2_col_dens = []

    while (id1 < len(uvb1["col_dens"])) and (id2 < len(uvb2["col_dens"])):
        
        # checking for instances where uvb2 has split clumps relative to uvb1
        if id1 in split:
            
            for key in uvb1_out.keys():
                uvb1_out[key].append(uvb1[key][id1])

            # uvb1_col_dens.append(uvb1["col_dens"][id1])
            
            clump_piece_ids = split[id1]
            
            nsplit[id1] = len(clump_piece_ids)

            for key in uvb2_out.keys():
                uvb2_col_dens_avg = np.mean(uvb2[key].iloc[clump_piece_ids])
                uvb2_out[key].append(uvb2_col_dens_avg)
            
            # uvb2_col_dens_avg = np.mean(uvb2["col_dens"].iloc[clump_piece_ids])
            # uvb2_col_dens.append(uvb2_col_dens_avg)
            
            id2 += (len(clump_piece_ids) - 1)
        
        # checking for instances where uvb2 has merged clumps relative to uvb1
        elif id2 in merge:
            
            for key in uvb2_out.keys():
                uvb2_out[key].append(uvb2[key][id2])

            # uvb2_col_dens.append(uvb2["col_dens"][id2])
            
            clump_piece_ids = merge[id2]
            
            nmerge[id2] = len(clump_piece_ids)

            for key in uvb1_out.keys():
                uvb1_col_dens_avg = np.mean(uvb1[key].iloc[clump_piece_ids])
                uvb1_out[key].append(uvb1_col_dens_avg)

            # uvb1_col_dens_avg = np.mean(uvb1["col_dens"].iloc[clump_piece_ids])
            # uvb1_col_dens.append(uvb1_col_dens_avg)
            
            id1 += (len(clump_piece_ids) - 1)
        
        # handling lonely cases
        elif id1 in lonely_1:
            
            for key in uvb1_out.keys():
                uvb1_out[key].append(uvb1[key][id1])
                uvb2_out[key].append(0)

            # uvb1_col_dens.append(uvb1["col_dens"][id1])
            # uvb2_col_dens.append(0)
    
        elif id2 in lonely_2:
            
            for key in uvb1_out.keys():
                uvb2_out[key].append(uvb2[key][id2])
                uvb1_out[key].append(0)
            
            # uvb2_col_dens.append(uvb2["col_dens"][id2])
            # uvb1_col_dens.append(0)
        
        # if none of the above cases pass, then the clump indices line up
        else:

            for key in uvb1_out.keys():
                uvb1_out[key].append(uvb1[key][id1])
                uvb2_out[key].append(uvb2[key][id2])
            
            # uvb1_col_dens.append(uvb1["col_dens"][id1])
            # uvb2_col_dens.append(uvb2["col_dens"][id2])
             
        id1+=1
        id2+=1
    
    # If everything worked correctly, the two output arrays should be the same length
    for key in uvb1_out.keys():
        assert len(uvb1_out[key]) == len(uvb2_out[key]), f"{key} arrays are different sizes. \n uvb1 = {len(uvb1_out[key])} \n uvb2 = {len(uvb2_out[key])}"
    
    out_sort_list = sorted_list
    out_sort_list[3] = nsplit
    out_sort_list[4] = nmerge

    return uvb1_out, uvb2_out, out_sort_list
    