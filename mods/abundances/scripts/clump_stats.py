import pickle
import numpy as np
import pandas as pd
from scipy import stats

raynum = 4 ##number of rays we have
ion = "C_IV" ##which ion to analyse

##necessary lists
ray_nums = []
super_cl_nums = []
med_col_dens = []
mad_for_med = []
#distances = []
central_v = []
vel_dispersions =[]
densities = []
temperatures = []
num_clumps = []
rows_of_rep_clumps = []

##necessary functions
def weighted_av(values, weights):
    weighted_sum = []
    for value, weight in zip(values, weights):
        weighted_sum.append(value * weight)
    
    return sum(weighted_sum) / sum(weights)

def make_full_list(list_in, list_out):
        for element in list_in:
            list_out.append(element)
        return list_out

for r in range(raynum):
    pickle_match_off = open(f"Match_{ion}_Ray{r}.pickle", 'rb') ##get data from id_clumps
    match = pickle.load(pickle_match_off)

    pickle_split_off = open(f"Split_{ion}_Ray{r}.pickle", 'rb')
    split = pickle.load(pickle_split_off)

    pickle_short_off = open(f"Short_{ion}_Ray{r}.pickle", 'rb')
    short = pickle.load(pickle_short_off)

    super_clumps = np.load(f'super_clumps_array_{ion}_ray{r}.npy')

    var_rows = []
    path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/"
    datanum = 26 ##number of rows on the abundance table, modify as needed
    ndigits= len(str(datanum))
    for i in range(datanum):
        m = i+1 ##correct filenames
        n_len = len(str(m))
        n_zeros = ndigits - n_len
        p = "0" * n_zeros + str(m)
        row_data = pd.read_csv(path+f"data_AbundanceRow{p}_{ion}.txt", delim_whitespace=True) ##read in data files
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
                

    for k in range(len(sup_st)): ##depending on which category each clump belongs to in super_clumps, append its column density to a list
        col_density_match = []
        col_density_split = []
        col_density_short = []
        match_done = [] ##keep track of where things have and haven't been accounted for
        short_done = []
        split_done = []

        for row, index in match.items():
        
            for j in range(len(index)): ##iterate over the list of indecies in the dictionary

                if (index[j][0]>=sup_st[k]) and (index[j][1]<=sup_en[k]): ##how to differentaite between super_clumps
                    ds = var_rows[row-1]
                    indexq = np.where((index[j][0]) == (ds["interval_start"])) ##find the correct row in the dataset
                    col_density_match.append(ds["col_dens"][int(indexq[0])]) 
                
                    if sup_st[k] not in match_done: ##for each super_clump, get each of these quantites once
                        print(row, index[j][0], sup_st[k])
                        #distances.append(ds["radius"][int(indexq[0])])
                        central_v.append(ds["delta_v"][int(indexq[0])])
                        vel_dispersions.append(ds["vel_dispersion"][int(indexq[0])])
                        densities.append(ds["density"][int(indexq[0])])
                        temperatures.append(ds["temperature"][int(indexq[0])])
                        rows_of_rep_clumps.append(row)
                        match_done.append(sup_st[k])
                        

        for rows, indexs in short.items(): ##same procedure as match
            
            for j in range(len(indexs)):
                
                if (indexs[j][0]>=sup_st[k]) and (indexs[j][1]<=sup_en[k]):
                    ds = var_rows[rows-1]
                    indexq = np.where((indexs[j][0]) == (var_rows[rows-1]["interval_start"]))
                    col_density_short.append(ds["col_dens"][int(indexq[0])])
                    
                    if len(match_done) == 0 and (sup_st[k] not in short_done): ##only do this if there isn't a match
                        print(rows, indexs[j][0], sup_st[k])
                        #distances.append(ds["radius"][int(indexq[0])])
                        central_v.append(ds["delta_v"][int(indexq[0])])
                        vel_dispersions.append(ds["vel_dispersion"][int(indexq[0])])
                        densities.append(ds["density"][int(indexq[0])])
                        temperatures.append(ds["temperature"][int(indexq[0])])
                        rows_of_rep_clumps.append(row)
                        short_done.append(sup_st[k])

        for rowm, indexm in split.items(): ##split is a bit weird
        ##if necessary to use a split as a "representative, these lists are how we make the weighted average
            temp_col_dens =[]
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
                        temp_col_dens.append(10 ** ds["col_dens"][int(indexq[0])]) ##column density must be summed to give a coreect value
                        col_dens_for_weights.append(ds["col_dens"][int(indexq[0])])
                        temp_delta_v.append(ds["delta_v"][int(indexq[0])])
                        temp_vel_dis.append(ds["vel_dispersion"][int(indexq[0])])
                        temp_dens.append(ds["density"][int(indexq[0])])
                        temp_temp.append(ds["temperature"][int(indexq[0])])
                        #temp_rad.append(ds["radius"][int(indexq[0])])
						
                        
                        if len(match_done) == 0 and len(short_done) == 0 and sup_st[k] not in split_done: ##only do this if no match or shorter
                            #distances.append(weighted_av(temp_rad, col_dens_for_weights))
                            central_v.append(weighted_av(temp_delta_v, col_dens_for_weights))
                            vel_dispersions.append(weighted_av(temp_vel_dis, col_dens_for_weights))
                            densities.append(weighted_av(temp_dens, col_dens_for_weights))
                            temperatures.append(weighted_av(temp_temp, col_dens_for_weights))
                            rows_of_rep_clumps.append(row)
                            split_done.append(sup_st[k])
                            
                    if len(temp_col_dens) != 0:
                        log_sum_dens = np.log10(sum(temp_col_dens))
                        col_density_split.append(log_sum_dens)

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
        
##make and fill the dictionary that becomes our dataframe
clump_stats = {}
clump_stats["ray_num"] = ray_nums
clump_stats["super_clump_number"] = super_cl_nums 
clump_stats["median_col_desnity"] = med_col_dens 
clump_stats["mad_for_col_desnity"] = mad_for_med
#clump_stats["distance_from_galaxy"] = distances 
clump_stats["central_velocity"] = central_v 
clump_stats["vel_dispersion"] = vel_dispersions 
clump_stats["density"] = densities 
clump_stats["temperature"] = temperatures
clump_stats["num_of_clumps"] = num_clumps
clump_stats["rep_clumps_row"] = rows_of_rep_clumps

df = pd.DataFrame.from_dict(clump_stats)
df.to_csv(f"halo_rshift_{ion}_abun_model-family_all-clumps.csv" ,sep = ' ')
