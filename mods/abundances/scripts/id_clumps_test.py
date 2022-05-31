import pandas as pd
import numpy as np
import pickle

df1=pd.read_csv("/mnt/scratch/f0104093/condensed_pipeline_tests/data/data_AbundanceRow09_C_IV.txt", delim_whitespace=True) ##read in data files
df2=pd.read_csv("/mnt/scratch/f0104093/condensed_pipeline_tests/data/data_AbundanceRow10_C_IV.txt", delim_whitespace=True)

df1_work=df1[df1["lightray_index"]==1] ##filter to only ray1
df2_work=df2[df2["lightray_index"]==1]
df1_clumps = df1_work[["interval_start","interval_end"]].reset_index().drop(columns="index")  ##filter to only indexes
df2_clumps = df2_work[["interval_start","interval_end"]].reset_index().drop(columns="index") 

mx= -np.inf  ##find how long each array should be
rowlist = [df1_clumps, df2_clumps]
for ds in rowlist: #find the cell index of the furthest clump
  row_mx = max(ds["interval_end"])
  if row_mx>mx:
    mx=row_mx
  
super_clumps = np.zeros(int(mx))

clmaps = []
for ds in rowlist: ##make masks for each row and form super_clumps
  ds_clump_loc = np.zeros(int(mx))
  for i in range(ds.shape[0]):
    print(ds["interval_start"][i])
    ds_clump_loc[int(ds["interval_start"][i]):int(ds["interval_end"][i])] = 1
    super_clumps[int(ds["interval_start"][i]):int(ds["interval_end"][i])] = 1
  clmaps.append(ds_clump_loc)
super_clumps=np.append(0, super_clumps)  ##make indexing work  
super_clumps=np.append(super_clumps, 0)
                    
match = {} ##create dictionaries to store indexes of clumps in the row that correspond to one another, keys will be row numbers and values will be indecies except for the lonlies
short = {}
merge = {}
lonely = {}

maybe_lonely = {}

rownum = 0

for row in clmaps:  ##start by iterating over a whole row

    row = np.append(0,row) # adding an extra element to prevent booleans from failing
    row = np.append(row,0)
    rownum += 1 ##define which row we're working on

    ##make all the variables and lists necessary
    row_st_cnt = 0 ##count how many starts of row clumps there have been within a super clump
    row_st_ind = []  ##keep track of start location(s) of a row clump
    
    row_en_cnt = 0 ##how many ends of row clumps there have been within a super clump
    row_en_ind = [] ##keep track of end location(s) of a row clump

    row_match = []
    row_short = []
    row_merge = []
    
 
    
    for i in range(1,len(row)):
        
        if super_clumps[i-1]<super_clumps[i]: ##start of a super clump
            sup_st = i-1 ##keep track of start location of a super clump
            
        if row[i-1]<row[i]:  ##start of a clump in the row
            row_st_cnt += 1
            row_st_ind.append(i-1)
        
        elif row[i-1]>row[i]: ##end of a clump in row
            row_en_cnt += 1
            row_en_ind.append(i-1)
            
        if super_clumps[i-1]>super_clumps[i]: ##end of a super clump
            sup_en = i-1 ##keep track of the location of the end of a super clump
            
            if (row_st_cnt == 1) and (row_en_cnt == 1): ##check for if there is only one row clump in the super clump
    
                if (row_st_ind[0] == sup_st) & (row_en_ind[0] == sup_en): ##if the starts and ends match, the clumps are identical
                    row_match.append([row_st_ind[0],row_en_ind[0]]) ##thus, start and end indecies appended to a list of them
            
                else:
                    row_short.append([row_st_ind[0],row_en_ind[0]]) ##if not, then the row clump must be shorter and the start and end indecies are appended to the appropraite list
            
            elif (row_st_cnt == 0) & (row_en_cnt == 0): ##check if there is nothing in the row that matches the super clump
                
                if str([sup_st,sup_en]) in maybe_lonely.keys(): ##check if we have already seen this superclump, if not make the entry in the dictionary
                    maybe_lonely[str([sup_st,sup_en])] += 1
                                 
                else:
                    maybe_lonely[str([sup_st,sup_en])] = 1
            
            else: ##only other senario is there there was a merge
            
                for j in range(len(row_st_ind)): ##organize the indecies to make the list in order
                    row_merge.append([row_st_ind[j],row_en_ind[j]]) 
                                 
            row_st_cnt = 0
            row_st_ind = []
            
            row_en_cnt = 0
            row_en_ind = []
        
        match[rownum] = row_match
        short[rownum] = row_short
        merge[rownum] = row_merge

                          
#for clump in maybe_lonely: ##later must set a limits on the minimumm number of times something has to appear in maybe lonely for it to actually be considered lonely


pickling_match = open("match.pickle","wb") ##saves the dictonaries so that they can be accesssed later
pickle.dump(match, pickling_match, protocol=3)
pickling_match.close()
                          
pickling_merge = open("merge.pickle","wb")
pickle.dump(merge, pickling_merge, protocol=3)
pickling_merge.close() 
                          
pickling_short = open("short.pickle","wb")
pickle.dump(short, pickling_short, protocol=3)
pickling_short.close()
  
