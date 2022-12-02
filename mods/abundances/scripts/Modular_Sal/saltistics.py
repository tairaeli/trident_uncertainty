import numpy as np
import pandas as pd
from scipy import stats

def make_full_list(list_in, list_out):
    for element in list_in:
        list_out.append(element)
    return list_out


def collect_ray_dat(salsa_out_dict, ion, r):
    
    """
    Looks through SALSA output and isolates data on a specific ray for each row
    Puts that data into a list
    Each entry in list is a certain row containing only the info from
    the specified ray
    """
    
    var_rows = []
    
    num_rows = len(salsa_out_dict[ion].keys())
    
    for row in range(num_rows):
    
        row_dat = salsa_out_dict[ion][row]
        
        ray_dat = row_dat[row_dat["lightray_index"] == str(r)].reset_index().drop(columns="index")
        
        var_rows.append(ray_dat)
    
    return var_rows
    
def weighted_av(values, weights): ##define functions necessary in making statistics
    weighted_sum = []
    for value, weight in zip(values, weights):
        weighted_sum.append(value * weight)

    return sum(weighted_sum) / sum(weights)


def abundance_stats(compare_dict, salsa_out_dict, ion_list, raynum, stat_path, halo, true_rs):
    
    stat_dict = {}
    
    for ion in ion_list:
        
        stat_dict[ion] = {}
            
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
            
            match = compare_dict[ion][r]["match"]
            short = compare_dict[ion][r]["short"]
            split = compare_dict[ion][r]["split"]
            
            super_clumps = compare_dict[ion][r]["super_clump"]
            
            var_rows = collect_ray_dat(salsa_out_dict, ion, r)
            
            
            
            # for i in range(datanum):
            #     m = i+1
            #     n_len = len(str(m))
            #     n_zeros = ndigits - n_len
            #     p = "0" * n_zeros + str(m)
            #     row_data = pd.read_csv(dat_path+f"/data_AbundanceRow{p}_{ion}.txt", delim_whitespace=True) ##read in data files
            #     row_work = row_data[row_data["lightray_index"]==r] ##filter to only one ray
            #     df = row_work.reset_index().drop(columns="index") ##make indexing work
            #     var_rows.append(df)
    
            sup_st = [] ##get indexes of superclumps
            sup_en = []
            for i in range(1, len(super_clumps)):
                n = i-1
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
                            
                            # could this pose issues at row = 0 ?
                            ds = var_rows[row-1]
                            indexq = np.where((index[j][0]) == (ds["interval_start"]))
                            
                            if row != 1: ##exclude solar abundance
                            
                                col_density_match.append(ds["col_dens"][int(indexq[0])]) ##get column density
                                
                            elif row == 1: ##keep frack of what the solar abundance is 
                                # print(ds["col_dens"])
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
                                print(ds["col_dens"])
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
        
        stat_dict[ion] = pd.DataFrame(clump_stats)
        stat_dict[ion].to_csv(f"{stat_path}/{halo}_z{true_rs}_{ion}_abun_all-model-families_all-clumps.csv" ,sep = ' ', na_rep = 'NaN') ##save the files to scratch
        
    return stat_dict
        
        
