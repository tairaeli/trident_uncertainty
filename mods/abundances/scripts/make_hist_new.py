import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

pickle_match_off = open("MatchCIIRay0.pickle", 'rb')
match = pickle.load(pickle_match_off)

pickle_split_off = open("SplitCIIRay0.pickle", 'rb')
split = pickle.load(pickle_split_off)

pickle_short_off = open("ShortCIIRay0.pickle", 'rb')
short = pickle.load(pickle_short_off)

super_clumps = np.load('super_clumps_array_C_II_ray0.npy')

var_rows = []
path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/"
datanum = 26 ##number of rows on the abundance table, modify as needed
ndigits= len(str(datanum))
for i in range(datanum):
  m = i+1
  n_len = len(str(m))
  n_zeros = ndigits - n_len
  p = "0" * n_zeros + str(m)
  row_data = pd.read_csv(path+f"data_AbundanceRow{p}_C_II.txt", delim_whitespace=True) ##read in data files
  row_work = row_data[row_data["lightray_index"]==0] ##filter to only one ray
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

    for row, index in match.items():
        for j in range(len(index)):
            if (index[j][0]==sup_st[k]) and (index[j][1]==sup_en[k]):
                ds = var_rows[row-1]
                indexq = np.where((index[j][0]) == (var_rows[row-1]["interval_start"]))
                col_density_match.append(ds["col_dens"][int(indexq[0])])

    for rowm, indexm in split.items(): ##split is a bit weird so we have to average the densities maybe should sum though?
        temp_col_dens =[]
        for j in range(len(indexm)):
            if (indexm[j][0]>=sup_st[k]) and (indexm[j][1]<=sup_en[k]):
                ds = var_rows[rowm-1]
                indexq = np.where((indexm[j][0]) == (var_rows[rowm-1]["interval_start"]))
                temp_col_dens.append(10 ** ds["col_dens"][int(indexq[0])])
        if len(temp_col_dens) != 0:
            log_sum_dens = np.log10(sum(temp_col_dens))
            col_density_split.append(log_sum_dens)

    for rows, indexs in short.items():
        for j in range(len(indexs)):
            if (indexs[j][0]>=sup_st[k]) and (indexs[j][1]<=sup_en[k]):
                ds = var_rows[rows-1]
                indexq = np.where((indexs[j][0]) == (var_rows[rows-1]["interval_start"]))
                col_density_short.append(ds["col_dens"][int(indexq[0])])

    
   ##plot the results##
    plt.hist((col_density_match, col_density_split, col_density_short), histtype='barstacked', label=['col_density_match', 'col_density_split', 'col_density_short', 'col_density_false_split'])
    plt.legend()
    plt.title(f"C_II -- LightRay Index 0 -- Super Clump {k}")
    plt.savefig(f"/mnt/scratch/f0104093/condensed_pipeline_tests/visuals/super_clump_hist/Hist_CII_RayIndex0_SuperClump{k}_.png")
    plt.close()
    print("Plotted!")


print("Go Look!")
