import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

pickle_match_off = open("match.pickle", 'rb')
match = pickle.load(pickle_match_off)

pickle_merge_off = open("merge.pickle", 'rb')
merge = pickle.load(pickle_merge_off)

pickle_short_off = open("short.pickle", 'rb')
short = pickle.load(pickle_short_off)

super_clumps = np.load('super_clumps_array.npy')

var_rows = []
path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/"
datanum = 26 ##number of rows on the abundance table, modify as needed
ndigits= len(str(datanum))
for i in range(datanum):
  m = i+1
  n_len = len(str(m))
  n_zeros = ndigits - n_len
  p = "0" * n_zeros + str(m)
  row_data = pd.read_csv(path+f"data_AbundanceRow{p}_C_IV.txt", delim_whitespace=True) ##read in data files
  row_work = row_data[row_data["lightray_index"]==1] ##filter to only ray1
  df = row_work.reset_index().drop(columns="index")
  ver_rows.append(df)

sup_st = []
sup_en = []
for i in range(1, len(super_clumps)):
    n= i-1
    if super_clumps[n]<super_clumps[i]: ##start of a super clump
        sup_st.append(n) 
    elif super_clumps[n]>super_clumps[i]: ##end of a super clump
        sup_en.append(n)
	
for k in range(len(sup_st)):
    col_density_match = []
    col_density_merge = []
    col_density_short = []
    for row, index in match.items():
        for j in range(len(index)):
            if (index[j][0]==sup_st[k]) and (index[j][1]==sup_en[k]):
                ds = var_rows[row-1]
                col_density_match.append(ds["col_dens"][j])
				
    for rowm, indexm in merge.items():
        for j in range(len(indexm)):
            if (indexm[j][0]>=sup_st[k]) and (indexm[j][1]<=sup_en[k]):
                ds = var_rows[rowm-1]
                col_density_merge.append(ds["col_dens"][j])
				
    for rows, indexs in short.items():
        for j in range(len(indexs)):
            if (indexs[j][0]>=sup_st[k]) and (indexs[j][1]<=sup_en[k]):
                ds = var_rows[rows-1]
                col_density_short.append(ds["col_dens"][j])
    sum_merge = sum(col_density_merge)
    plt.hist((col_density_match, sum_merge, col_density_short), histtype='barstacked', label=['col_density_match', 'col_density_merge', 'col_density_short'])
    plt.legend()
    plt.title(f"C_IV -- LightRay Index 1 -- Super Clump {k}")
    plt.savefig(f"/mnt/scratch/f0104093/condensed_pipeline_tests/visuals/super_clump_hist/Hist_CIV_RayIndex1_SuperClump{k}_.png")
    plt.close()
    print("Plotted!")

print("Go Look!")
