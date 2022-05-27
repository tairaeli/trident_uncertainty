import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import argparse
# %matplotlib inline
import sys
import os

df1=pd.read_csv("/mnt/scratch/f0104093/condensed_pipeline_data/data_AbundanceRow09_C_IV.txt", sep = " ") ##read in data files
df2=pd.read_csv("/mnt/scratch/f0104093/condensed_pipeline_data/data_AbundanceRow10_C_IV.txt", sep =" ")

df1_work=df1[df1["lightray_index"]==1] ##filter to only ray1
df2_work=df2[df2["lightray_index"]==1]
df1_clumps = df1_work[["interval_start","interval_end"]] ##filter to only indexes
df2_clumps = df2_work[["interval_start","interval_end"]]
df1_st = np.asarray(df1_clumps["interval_start"]) ##make each interval start and end into lists for easy indexing
df2_st = np.asarray(df2_clumps["interval_start"])
df1_en = np.asarray(df1_clumps["interval_end"])
df2_en = np.asarray(df2_clumps["interval_end"])

mx= -np.inf 
rowlist = [df1_clumps, df2_clumps]
for ds in rowlist: #find the cell index of the furthest clump
  row_mx = max(ds["interval end"])
  if row_mx>mx:
    mx=row_mx
  
super_clumps = np.zeros(int(mx))

clmaps = []
for ds in rowlist:
  ds_clump_loc = np.zeros(int(mx))
  for i in range(ds.shape[0]):
    ds_clump_loc[int(ds["interval_start"][i]):int(ds["interval_end"][i])] = 1
    super_clumps[int(ds["interval_start"][i]):int(ds["interval_end"][i])] = 1
  clmaps.append(ds_clump_loc)
super_clumps.append(0, super_clumps)    

one_locs = np.where(super_clumps == 1)[0]
start_list = []
end_list = []
for i in range(len(one_locs)):
  if one_locs[i+1] - one_locs[i] > 1:
    start_list.append(one_locs[i+1])
    end_list.append(one_locs[i])
                    
match = {} ##create dictionaries to store indexes of clumps in the row that correspond to one another, keys will be row numbers and values will be indecies
short = {}
merge = {}
lonely = {}
maybe_lonely = {}

rownum = 0
for row in clmaps:
  rownum += 1
  row=np.append(0, row)
  
  row_st_count = 0
  row_end_cnt = 0
  row_st_ind = []
  row_en_ind = []
  
  row_match=[]
  row_short = []
  row_merge = []
  
  if super_clumps[i-1] < super_clumps[i]:
    sp_st = i-1
    
  if row[i-1] < row[i]:
    row_st_ind.append(i-1)
    row_st_cnt += 1
    
  elif row[i-1] > row[i]:
    row_en_ind.append(i-1)
    row_en_cnt += 1
  
  if super_clumps[i-1] > super_clumps[i]:
    sp_en = i-1
    
    if (row_st_count == 1) and (row_en_cnt == 1):
      if (row_st_ind[0] == sp_st) and (row_en_ind[0] == sp_en):
        row_match.append([row_st_ind[0], row_en_ind[0]])
      else:
        row_short.append([row_st_ind[0], row_end_ind[0]])
        
    elif (row_st_cnt == 0) and (row_en_cnt ==0):
      if str([sp_st, sp_en]) in maybe_lonely.keys():
        maybe_lonely[str([sp_st, sp_en])] += 1
      else:
        maybe_lonely[str([sp_st, sp_en])] = 1
        
    else:
      for j in range(len(row_st_ind)):
        row_merge.append([row_st_ind[j], row_en_ind[j]
    
    row_st_count = 0
    row_en_count = 0
    row_st_ind = []
    row_en_ind = []

  match[rownum]=row_match
  short[rownum]=row_short
  merge[rownum]=row_merge
                          
for clump in maybe_lonely:
   
  
                         
  
##err = 10 ##set allowable amount of error
# n = 0
# for j in range(len(df2_clumps)):
#   split_inst = []
#   merge_inst =[]
#   if df1_en[n] == df2_en[j] and df1_st[n] == df2_st[j]: ##conditions for identical matches
#     match[str(n)] = j
#     n += 1
#   elif ((df2_en[j] - err) <= df1_en[n] <= (df2_en[j] + err)) or ((df2_st[j] - err) <= df1_st[n] <= (df2_st[j] + err)): ##conditions for slight errors
#     slight_off[str(n)] = j
#     n += 1
#   ##elif ((df1_st[n] == df2_st[j]) and (df1_en[n] != df2_en[j]) or (df1_en[n] == df2_en[j] and df1_st[n] != df2_st[j]) or ## commented out so still able to be seen
#   elif df2_st[j] >= df1_st[n] and df2_en[j] <= df1_en[n]: ##conditions for a split 
#     split_inst.append(j)
#     if df1_en[n] == df2_en[j]:
#        split[str(n)]= split_inst
#        split_inst =[]
#        n+=1
#     else:
#        continue
#   elif df1_st[n] >= df2_st[i] and df1_en[n] <= df2_en[i]: ##conditions for a merge 
#     merge_inst.append(j)
#     if df1_en[n] == df2_en[j]:
#        merge[str(n)]= merge_inst
#        merge_inst =[]
#        n+=1
#     else:
#        continue
#   elif :  
