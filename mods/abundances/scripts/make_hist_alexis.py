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
for i in range(1, len(super_clumps)+1):
	n= i-1
        if super_clumps[n]<super_clumps[i]: ##start of a super clump
            sup_st.append(n) 
	elif super_clumps[n]>super_clumps[i]: ##end of a super clump
            sup_en.append(n)
	
var_rows = [df1_work, df2_work]
col_density_match = []
for row, index in match.items():
 	for j in range(len(index)):
 		if index[j][0]>=sup_st[0] and index[j][1]<=sup_st[0]:
			ds = var_rows[row-1]
			col_density_match.append(ds["col_dens"][j])

##to create stacked histogram, really do a stacked bar plot, https://www.tutorialexample.com/matplotlib-create-stacked-histogram-a-beginner-guide-matplotlib-tutorial/
# col_density_match = []
# 
# 			name = ##figure out how to make the naming things work from string to variable
# 			col_density_match.append(df{row}_work["col_dens"][row_clump])

plt.hist(col_density_match)
plt.title(f"C_IV -- LightRay Index 1 -- Matches")
plt.savefig(f"/mnt/scratch/f0104093/condensed_pipeline_tests/data/Hist_CIV_RayIndex1_Matches.png")

            
