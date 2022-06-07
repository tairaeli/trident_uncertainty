import numpy as np
import pandas as pd

df1 = pd.read_csv("./data/data_AbundanceRow08_C_IV.txt", delim_whitespace=True)
df2 = pd.read_csv("./data/data_AbundanceRow09_C_IV.txt", delim_whitespace=True)

df1_work=df1[df1["lightray_index"]==1] ##filter to only ray1
df2_work=df2[df2["lightray_index"]==1]

df1_clumps = df1_work[["interval_start","interval_end"]].reset_index().drop(columns = "index") ##filter to only clump locations
df2_clumps = df2_work[["interval_start","interval_end"]].reset_index().drop(columns = "index")

rowlist = [df1_clumps,df2_clumps]

mx = -np.inf
for ds in rowlist:
    row_mx = max(ds["interval_end"])
    if row_mx > mx:
        mx = row_mx
        
super_clumps = np.zeros(int(mx))

clmaps = []
for ds in rowlist:
    ds_clump_loc = np.zeros(int(mx))
    for i in range(ds.shape[0]):
        ds_clump_loc[int(ds['interval_start'][i]):int(ds['interval_end'][i])] = 1
        super_clumps[int(ds['interval_start'][i]):int(ds['interval_end'][i])] = 1
    clmaps.append(ds_clump_loc)

super_clumps = np.append(0,super_clumps)
super_clumps = np.append(super_clumps,0)
super_clumps

match = {} ##create dictionaries to store indexes of clumps that correspond to one another
shorter = {}
merge = {}
lonely = {}

maybe_lonely = {}

rownum = 0

for row in clmaps: 

    row = np.append(0,row) # adding an extra element to prevent booleans from failing
    row = np.append(row,0)
    rownum += 1

    row_st_cnt = 0
    row_st_ind = []
    
    row_en_cnt = 0
    row_en_ind = []

    row_match = []
    row_short = []
    row_merge = []
    
 
    
    for i in range(1,len(row)):
        
        if super_clumps[i-1]<super_clumps[i]:
            sup_st = i-1
            
        if row[i-1]<row[i]:
            row_st_cnt += 1
            row_st_ind.append(i-1)
        
        elif row[i-1]>row[i]:
            row_en_cnt += 1
            row_en_ind.append(i-1)
            
        if super_clumps[i-1]>super_clumps[i]:
            sup_en = i-1
            
            if (row_st_cnt == 1) & (row_en_cnt == 1):
                if (row_st_ind[0] == sup_st) & (row_en_ind[0] == sup_en):
                    row_match.append([row_st_ind[0],row_en_ind[0]])
            
                else:
                    row_short.append([row_st_ind[0],row_en_ind[0]])
            
            elif (row_st_cnt == 0) & (row_en_cnt == 0):
                
                if str([sup_st,sup_en]) in maybe_lonely.keys():
                    maybe_lonely[str([sup_st,sup_en])] += 1
                                 
                else:
                    maybe_lonely[str([sup_st,sup_en])] = 1
            
            else:
            
                for j in range(len(row_st_ind)):
                    row_merge.append([row_st_ind[j],row_en_ind[j]]) 
                                 
            row_st_cnt = 0
            row_st_ind = []
            
            row_en_cnt = 0
            row_en_ind = []
        
        match[rownum] = row_match
        shorter[rownum] = row_short
        merge[rownum] = row_merge