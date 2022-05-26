import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import argparse
# %matplotlib inline
import sys
import os

df1=pd.read_csv("/mnt/scratch/f0104093/condensed_pipeline_data/data_AbundanceRow9_C_IV.txt", sep = " ") ##read in data files
df2=pd.read_csv("/mnt/scratch/f0104093/condensed_pipeline_data/data_AbundanceRow10_C_IV.txt", sep =" ")

df1_work=df1[df1["lightray_index"]==1] ##filter to only ray1
df2_work=df2[df2["lightray_index"]==1]

df1_clumps = df1_work[["interval_start","interval_end"]] ##filter to only indexes
df2_clumps = df2_work[["interval_start","interval_end"]]
df1_st = np.asarray(df1_clumps["interval_start"]) ##make each interval start and end into lists for easy indexing
df2_st = np.asarray(df2_clumps["interval_start"])
df1_en = np.asarray(df1_clumps["interval_end"])
df2_en = np.asarray(df2_clumps["interval_end"])

match = {} ##create dictionaries to store indexes of clumps that correspond to one another
slight_off = {}
merge = {}
split = {}

err = 10 ##set allowable amount of error
n = 0
for j in range(len(df2_clumps)):
  split_inst = []
  merge_inst =[]
  if df1_en[n] == df2_en[j] and df1_st[n] == df2_st[j]: ##conditions for identical matches
    match[str(n)] = j
    n += 1
  elif ((df2_en[j] - err) <= df1_en[n] <= (df2_en[j] + err)) or ((df2_st[j] - err) <= df1_st[n] <= (df2_st[j] + err)): ##conditions for slight errors
    slight_off[str(n)] = j
    n += 1
  ##elif ((df1_st[n] == df2_st[j]) and (df1_en[n] != df2_en[j]) or (df1_en[n] == df2_en[j] and df1_st[n] != df2_st[j]) or ## commented out so still able to be seen
  elif df2_st[j] >= df1_st[n] and df2_en[j] <= df1_en[n]: ##conditions for a split 
    split_inst.append(j)
    if df1_en[n] == df2_en[j]:
       split[str(n)]= split_inst
       split_inst =[]
       n+=1
    else:
       continue
  elif df1_st[n] >= df2_st[i] and df1_en[n] <= df2_en[i]: ##conditions for a merge 
    merge_inst.append(j)
    if df1_en[n] == df2_en[j]:
       merge[str(n)]= merge_inst
       merge_inst =[]
       n+=1
    else:
       continue
  elif :  
