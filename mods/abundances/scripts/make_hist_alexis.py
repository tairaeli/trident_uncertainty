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

path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/" # modify as needed
data1 = pd.read_csv(path + "data_AbundanceRow09_C_IV.txt", sep = " ") ##read in data files
data2 = pd.read_csv(path + "data_AbundanceRow10_C_IV.txt", sep =" ")
df1_work=data1[data1["lightray_index"]==1].reset_index().drop(columns="index")  ##filter to only ray1
df2_work=data2[data2["lightray_index"]==1].reset_index().drop(columns="index") 

sup_st = []
sup_en = []
for i in range(1, len(super_clumps)):
    n= i-1
    if super_clumps[n]<super_clumps[i]: ##start of a super clump
        sup_st.append(n) 
    elif super_clumps[n]>super_clumps[i]: ##end of a super clump
        sup_en.append(n)
	
var_rows = [df1_work, df2_work]
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
    plt.hist((col_density_match, col_density_merge, col_density_short), histtype='barstacked', label=['col_density_match', 'col_density_merge', 'col_density_short'])
    plt.legend()
    plt.title(f"C_IV -- LightRay Index 1 -- Super Clump {k}")
    plt.savefig(f"/mnt/scratch/f0104093/condensed_pipeline_tests/visuals/super_clump_hist/Hist_CIV_RayIndex1_SuperClump{k}_.png")
    plt.close()
    print("Plotted!")

print("Go Look!")
