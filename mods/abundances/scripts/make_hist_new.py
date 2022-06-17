import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

r = 3 ##which ray to work with

pickle_match_off = open(f"MatchCIIRay{r}.pickle", 'rb')
match = pickle.load(pickle_match_off)

pickle_split_off = open(f"SplitCIIRay{r}.pickle", 'rb')
split = pickle.load(pickle_split_off)

pickle_short_off = open(f"ShortCIIRay{r}.pickle", 'rb')
short = pickle.load(pickle_short_off)

super_clumps = np.load(f'super_clumps_array_C_II_ray{r}.npy')

var_rows = [] ##list of all the data sets
path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/"
datanum = 26 ##number of rows on the abundance table, modify as needed
ndigits= len(str(datanum))
for i in range(datanum):
    m = i+1
    n_len = len(str(m))
    n_zeros = ndigits - n_len
    p = "0" * n_zeros + str(m)
    row_data = pd.read_csv(path+f"data_AbundanceRow{p}_C_II.txt", delim_whitespace=True) ##read in data files
    row_work = row_data[row_data["lightray_index"]==r] ##filter to only one ray
    df = row_work.reset_index().drop(columns="index") ##make indexing work
    var_rows.append(df)

sup_st = [] ##get indexes of superclumps
sup_en = []
for i in range(1, len(super_clumps)):
    n= i-1
    if super_clumps[n]<super_clumps[i]: ##start of a super clump
        if super_clumps[n]== 2: ##edge case handling
            sup_st.append(n-1)
        else: 
            sup_st.append(n)
    elif super_clumps[n]>super_clumps[i]: ##end of a super clump
        if super_clumps[n] ==2:
            sup_en.append(n-1)
        else:
            sup_en.append(n)

	
for k in range(len(sup_st)): ##depending on which category each clump belongs to in super_clumps, append its column density to a list, make one plot per super clump
    col_density_match = [] ##lists to hold column densities 
    col_density_split = []
    col_density_short = []

    for row, index in match.items(): 
        for j in range(len(index)): ##over all the indices in a row in the dictionary
            
            if (index[j][0]>=sup_st[k]) and (index[j][1]<=sup_en[k]): ##filter to within the clump
                
                ds = var_rows[row-1] ##find the dataset
                indexq = np.where((index[j][0]) == (ds["interval_start"])) ##get the row of the clump within the dataset    
                col_density_match.append(ds["col_dens"][int(indexq[0])]) ##append to list

    for rowm, indexm in split.items(): ##split is a bit weird so we have to sum the densities, otherwise just like match
        temp_col_dens =[] ##define the list to sum
        for j in range(len(indexm)):
            if (indexm[j][0]>=sup_st[k]) and (indexm[j][1]<=sup_en[k]):
                ds = var_rows[rowm-1]
                indexq = np.where((indexm[j][0]) == (var_rows[rowm-1]["interval_start"]))
                if len(list(indexq[0])) > 0:
                    temp_col_dens.append(10 ** ds["col_dens"][int(indexq[0])]) ##have to add the actual numbers, not the log of the numbers (which is what salsa automaticlly gives you
        if len(temp_col_dens) != 0:
            log_sum_dens = np.log10(sum(temp_col_dens)) ##convert back to log to be able to compare
            col_density_split.append(log_sum_dens)

    for rows, indexs in short.items(): ##just like match
        for j in range(len(indexs)):
            if (indexs[j][0]>=sup_st[k]) and (indexs[j][1]<=sup_en[k]):
                ds = var_rows[rows-1]
                indexq = np.where((indexs[j][0]) == (var_rows[rows-1]["interval_start"]))
                col_density_short.append(ds["col_dens"][int(indexq[0])])

    
   ##plot the results##
    plt.hist((col_density_match, col_density_split, col_density_short), histtype='barstacked', label=['col_density_match', 'col_density_split', 'col_density_short'])
    plt.legend()
    plt.title(f"C_II -- LightRay Index {r} -- Super Clump {k}")
    plt.savefig(f"/mnt/scratch/f0104093/condensed_pipeline_tests/visuals/super_clump_hist/Hist_CII_RayIndex{r}_SuperClump{k}_.png")
    plt.close()
    print("Plotted!")


print("Go Look!")
