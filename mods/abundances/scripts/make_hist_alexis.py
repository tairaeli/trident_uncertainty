import pickle
import matplotlib as plt
import numpy as np
import pandas as pd

pickle_match_off = open("match.pickle", 'rb')
match = pickle.load(pickle_match_off)

pickle_merge_off = open("merge.pickle", 'rb')
merge = pickle.load(pickle_merge_off)

pickle_short_off = open("short.pickle", 'rb')
short = pickle.load(pickle_short_off)

path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/" # modify as needed
data1 = pd.read_csv(path + "data_AbundanceRow09_C_IV.txt", sep = " ") ##read in data files
data2 = pd.read_csv(path + "data_AbundanceRow10_C_IV.txt", sep =" ")
df1_work=data1[data1["lightray_index"]==1] ##filter to only ray1
df2_work=data2[data2["lightray_index"]==1]

col_density_match = []
for row, index in match:
	row_clump = np.where(index == list((df1_work["interval_start"]), (df1_work["interval_end"])))
	coldensity.append(df1_work["col_dens"][row_clump])
plt.hist(col_density_match)
plt.title(f"C_IV -- LightRay Index 1 -- Matches")
plt.savefig(f"/mnt/scratch/f0104093/condensed_pipeline_tests/data/Hist_CIV_RayIndex1_Matches.png")

            
