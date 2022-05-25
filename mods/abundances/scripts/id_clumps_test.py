import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import argparse
# %matplotlib inline
import sys
import os

df1=pd.read_csv("/mnt/scratch/f0104093/condensed_pipeline_data/data_AbundanceRow14_C_IV.txt", sep = " ") ##read in data files
df2=pd.read_csv("/mnt/scratch/f0104093/condensed_pipeline_data/data_AbundanceRow16_C_IV.txt", sep =" ")

df1_work=df1[df1["lightray_index"]==1] ##filter to only ray1
df2_work=df2[df2["lightray_index"]==1]

df1_clumps = df1_work[["interval_start","interval_end"]] ##filter to only indexes
df2_clumps = df2_work[["interval_start","interval_end"]]
df1_st = np.asarray(df1_clumps["interval_start"])
df2_st = np.asarray(df2_clumps["interval_start"])
df1_en = np.asarray(df1_clumps["interval_end"])
df1_en = np.asarray(df1_clumps["interval_end"])

match = {}
slight_off = {}
no_pair = {}

err = 10
n=0
for j in range(len(df2_clumps)):
  if df1_en[n] == df2_en[j] and df1_st[n] == df2_st[j]:
    match[str(n)] = j
    n += 1
  elif ((df2_en[j] - err) <= df1_en[n] <= (df2_en[j] + err)) or ((df2_st[j] - err) <= df1_st[n] <= (df2_st[j] + err)):
    slight_off[str(n)]=j
    n += 1
  elif 
