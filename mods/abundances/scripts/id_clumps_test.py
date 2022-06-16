import pandas as pd
import numpy as np
import pickle

path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/"
datanum = 26 ##number of rows on the abundance table, modify as needed
ndigits= len(str(datanum))
raynum = 4 ##number of rays used, modify as needed
ion  = "C_IV"

for r in range(1,2):
    rowlist = []
    for i in range(datanum):
        m = i+1
        n_len = len(str(m))
        n_zeros = ndigits - n_len
        k = "0" * n_zeros + str(m)
        row_data = pd.read_csv(path+f"data_AbundanceRow{k}_{ion}.txt", delim_whitespace=True) ##read in data files
        row_work = row_data[row_data["lightray_index"]==r] ##filter to only ray1
        df = row_work.reset_index().drop(columns="index") ##make indexing work
        rowlist.append(df)

    mx= -np.inf  ##find how long each array should be
    for ds in rowlist: #find the cell index of the furthest clump
        row_mx = max(ds["interval_end"])
        if row_mx>mx:
            mx=row_mx
    super_clumps = np.zeros(int(mx))
    clmaps = []
    
    problems = []
    hassles = {}
    
    row_tracker = 0
    for ds in rowlist: ##make masks for each row and form super_clumps
        ds_clump_loc = np.zeros(int(mx))
        row_tracker +=1
        hassles_ind = []
        for j in range(ds.shape[0]):
            ds_clump_loc[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1
                #print(f"Row Start: {int(ds['interval_start'][j])}, Row End: {int(ds['interval_end'][j])}")

            if j-1 == -1:
                super_clumps[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1

            elif int(ds["interval_end"][j-1]) == int(ds["interval_start"][j]):

                    #check for wether or not two clumps are within thier standard deviations
                if not (float(ds["delta_v"][j-1]) - 1.5 * (float(ds["vel_dispersion"][j-1]))) <= float(ds["delta_v"][j]) <= 1.5 * (float(ds["delta_v"][j-1]) + (float(ds["vel_dispersion"][j-1]))):
                    ds_clump_loc[int(ds["interval_start"][j])] = 2
                #edge case handling
                    if int(ds["interval_end"][j-1]) not in problems:
                        problems.append(int(ds["interval_end"][j-1])) 
                         #to make sure the bigger clumps stays in superclumps:
                    if (super_clumps[int(ds["interval_end"][j-1])] == 0) or (super_clumps[int(ds["interval_end"][j-1])] == 2):
                        super_clumps[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1
                        super_clumps[int(ds["interval_start"][j])] = 2
                else:
                    hassles_ind.append([int(ds["interval_start"][j-1]), int(ds["interval_end"][j-1])])
                    hassles_ind.append([int(ds["interval_start"][j]), int(ds["interval_end"][j])])                
                    super_clumps[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1
                    print(f"Hassles: {row_tracker, hassles_ind}")

            else:
                super_clumps[int(ds["interval_start"][j]):int(ds["interval_end"][j])] = 1  

        hassles[row_tracker]=hassles_ind
        clmaps.append(ds_clump_loc)

        for index in problems:
            for row in clmaps:
                if ((row[index-1] == 0 and row[index+1] != 0) or (row[index+1] == 0 and row[index-1] != 0)) and (row[index] != 0): ##edge case handling
                    super_clumps[index] = 2
                    
                
        
    
    super_clumps=np.append(0, super_clumps)  ##make indexing work  
    super_clumps=np.append(super_clumps, 0)
    np.save(f'super_clumps_array_{ion}_ray{r}', super_clumps) ##save super_clumps for future reference
    match = {} ##create dictionaries to store indexes of clumps in the row that correspond to one another, keys will be row numbers and values will be indecies except for the lonlies
    short = {}
    split = {}
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
        row_split = []
        sup_st_true =[]
        weird_index = 0

        for i in range(len(row)):
            
            if super_clumps[i-1]<super_clumps[i]: ##start of a super clump
                sup_st = i-1 ##keep track of start location of a super clump
                sup_st_true.append(sup_st)
                
            if row[i-1]<row[i] and (len(row_st_ind) == 0 or len(row_st_ind) ==1):##start of a clump in the row
                print(i-2, i-1, i)
                if super_clumps[i-1] == 2: ##if oops is printed, there is another edge case, debugging must begin again
                    print("oops")
                elif super_clumps[i] != 2 or len(row_st_ind)!=0: ##normal clumps look like this
                    row_st_ind.append(i-1)
                    row_st_cnt += 1
                    print(row_st_ind, rownum)
                else: ##edge case handling ##FIXME
                    weird_index = i-1

            elif row[i-1]>row[i]: ##end of a clump in row
                row_en_cnt += 1
                if row[i-1]==2: ## edge case handling
                    row_en_ind.append(i-2)
                else:
                    row_en_ind.append(i-1)
                
            if super_clumps[i-1]>super_clumps[i]: ##end of a super clump
                if super_clumps[i-1]==2:
                    sup_en = i-2
                    sup_st = sup_st_true[0] ##edge case handling
                    
                if super_clumps[i-1] == 0 or super_clumps[i-1] == 1:
                    sup_en = i-1 ##keep track of the location of the end of a super clump
                
               
               
                if (row_st_cnt == 1) or (row_en_cnt == 1): ##check for if there is only one row clump in the super clump
                    if weird_index != 0: ##edge case handling
                        row_st_ind.append(weird_index)
                        
                    if (len(row_en_ind) == 0 and len(row_st_ind) != 0) or (len(row_st_ind) == 0 and len(row_en_ind) != 0):
                        
                        break
                        
                    hassle_st = []
                    hassle_en = []
                    
                    if rownum in hassles.keys(): ##even more edge case handling
                        for value in hassles.values():
                            for t in value:
                                hassle_st.append(t[0])
                                hassle_en.append(t[1])
                            
                        if row_st_ind[0] in hassle_st:
                            for j in range(len(hassle_st)):
                                if [hassle_st[j], hassle_en[j]] not in row_split:
                                    row_split.append([hassle_st[j], hassle_en[j]])
                        elif row_en_ind[0] in hassle_en:
                            
                            continue
                        
                    if (row_st_ind[0] == sup_st) and (row_en_ind[0] == sup_en) and (row_st_ind[0] not in hassle_st) and (row_en_ind[0] not in hassle_en): ##if the starts and ends match, the clumps are identical
                
                        row_match.append([row_st_ind[0],row_en_ind[0]]) ##thus, start and end indecies appended to a list of them
                    
                    elif (row_st_ind[0] not in hassle_st) and (row_en_ind[0] not in hassle_en):
                        row_short.append([row_st_ind[0],row_en_ind[0]]) ##if not, then the row clump must be shorter and the start and end indecies are appended to the appropraite list
                        
                elif (row_st_cnt == 0) and (row_en_cnt == 0): ##check if there is nothing in the row that matches the super clump

                    if str([sup_st,sup_en]) in maybe_lonely.keys(): ##check if we have already seen this superclump, if not make the entry in the dictionary
                        maybe_lonely[str([sup_st,sup_en])] += 1
                    else:
                        maybe_lonely[str([sup_st,sup_en])] = 1
                        
                else: ##only other senario is there there was a split
                    for j in range(min([len(row_en_ind), len(row_st_ind)])): ##organize the indecies to make the list in order
                        row_split.append([row_st_ind[j],row_en_ind[j]]) 
                        
                        
                if super_clumps[i-1]==2:
                    sup_st=sup_st_true[1]
                    
                ##reset variables      
                row_st_cnt = 0
                row_st_ind = []
                row_en_cnt = 0
                row_en_ind = []
                sup_st_true = []
                match[rownum] = row_match
                short[rownum] = row_short
                split[rownum] = row_split
                
                if super_clumps[i-1] == 2:
                    sup_st = i-2
                    if row[i-1]>row[i]:
                        row_st_ind.append(i-2)
                        row_st_cnt += 1
                          
#for clump in maybe_lonely: ##later must set a limits on the minimumm number of times something has to appear in maybe lonely for it to actually be considered lonely

    pickling_match = open(f"Match_{ion}_Ray{r}.pickle","wb") ##saves the dictonaries so that they can be accesssed later
    pickle.dump(match, pickling_match, protocol=3)	
    pickling_match.close()
                          
    pickling_split = open(f"Split_{ion}_Ray{r}.pickle","wb")
    pickle.dump(split, pickling_split, protocol=3)
    pickling_split.close() 
                      
    pickling_short = open(f"Short_{ion}_Ray{r}.pickle","wb")
    pickle.dump(short, pickling_short, protocol=3)
    pickling_short.close()
    
    pickling_maybe_lonely = open(f"MaybeLonely_{ion}_Ray{r}.pickle","wb")
    pickle.dump(maybe_lonely, pickling_maybe_lonely, protocol=3)
    pickling_maybe_lonely.close()
