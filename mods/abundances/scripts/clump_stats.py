import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

raynum = 4 ##number of rays we have
ray_nums = []
super_cl_nums = []
med_col_dens = []
mad_for_med = []
distances = []
central_v = []
vel_dispersions =[]
densities = []
temperatures = []
num_clumps = []
rows_of_clump_reps = []

for r in range(raynum):
  pickle_match_off = open(f"MatchCIIRay{r}.pickle", 'rb')
  match = pickle.load(pickle_match_off)

  pickle_split_off = open(f"SplitCIIRay{r}.pickle", 'rb')
  split = pickle.load(pickle_split_off)

  pickle_short_off = open(f"ShortCIIRay{r}.pickle", 'rb')
  short = pickle.load(pickle_short_off)

  super_clumps = np.load(f'super_clumps_array_C_II_ray{r}.npy')

  var_rows = []
  path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/"
  datanum = 26 ##number of rows on the abundance table, modify as needed
  ndigits= len(str(datanum))
  for i in range(datanum):
      m = i+1
      n_len = len(str(m))
      n_zeros = ndigits - n_len
      p = "0" * n_zeros + str(m)
      row_data = pd.read_csv(path+f"data_AbundanceRow{p}_C_II.txt", delim_whitespace=True) ##read in data files
      row_work = row_data[row_data["lightray_index"]==r] ##filter to only one ray
      df = row_work.reset_index().drop(columns="index") ##make indexing work
      var_rows.append(df)

  sup_st = [] ##get indexes of superclumps
  sup_en = []
  for i in range(1, len(super_clumps)):
      n= i-1
      if super_clumps[n]<super_clumps[i]: ##start of a super clump
          if super_clumps[n]== 2:
              sup_st.append(n-1)
          else: 
              sup_st.append(n)
      elif super_clumps[n]>super_clumps[i]: ##end of a super clump
          if super_clumps[n] ==2:
              sup_en.append(n-1)
          else:
              sup_en.append(n)

def make_full_list(list_in, list_out):
   for element in list_in:
      list_out.append(element)
   return list_out

  for k in range(len(sup_st)): ##depending on which category each clump belongs to in super_clumps, append its column density to a list
      col_density_match = []
      col_density_split = []
      col_density_short = []
      match_num = 0
      split_num = 0
      short_num = 0

      for row, index in match.items():
          for j in range(len(index)):
            match_num += 1
            if (index[j][0]>=sup_st[k]) and (index[j][1]<=sup_en[k]):
              ds = var_rows[row-1]
              indexq = np.where((index[j][0]) == (ds["interval_start"]))
              if len(list(indexq[0]))>0:
                col_density_match.append(ds["col_dens"][int(indexq[0])])
                
                if match_num == 1:
                    distances.append(ds["radius"][int(indexq[0])])
                    central_v.append(ds["delta_v"][int(indexq[0])])
                    vel_dispersions.append(["vel_dispersion"][int(indexq[0])])
                    densities.append(["density"][int(indexq[0])])
                    temperatures.append(["temperature"][int(indexq[0])])
                    rows_of_rep_clumps.append(row)

      for rows, indexs in short.items():
          for j in range(len(indexs)):
				short_num += 1
              if (indexs[j][0]>=sup_st[k]) and (indexs[j][1]<=sup_en[k]):
                ds = var_rows[rows-1]
                indexq = np.where((indexs[j][0]) == (var_rows[rows-1]["interval_start"]))
				col_density_short.append(ds["col_dens"][int(indexq[0])])
				if match_num == 0 and short_num == 1:
                       			distances.append(ds["radius"][int(indexq[0])])
					central_v.append(ds["delta_v"][int(indexq[0])])
					vel_dispersions.append(["vel_dispersion"][int(indexq[0])])
					densities.append(["density"][int(indexq[0])])
					temperatures.append(["temperature"][int(indexq[0])])
					rows_of_rep_clumps.append(row)
					
      for rowm, indexm in split.items(): ##split is a bit weird so we have to average the densities maybe should sum though?
          temp_col_dens =[]
          for j in range(len(indexm)):
			split_num += 1 
            if (indexm[j][0]>=sup_st[k]) and (indexm[j][1]<=sup_en[k]):
				ds = var_rows[rowm-1]
                indexq = np.where((indexm[j][0]) == (var_rows[rowm-1]["interval_start"]))
                if len(list(indexq[0])) > 0:
			temp_col_dens.append(10 ** ds["col_dens"][int(indexq[0])])
					
			if match_num == 0 and short_num == 0 and 2 <= split_num <= 4:
				distances.append(ds["radius"][int(indexq[0])])
				central_v.append(ds["delta_v"][int(indexq[0])])
				vel_dispersions.append(["vel_dispersion"][int(indexq[0])])
				densities.append(["density"][int(indexq[0])])
				temperatures.append(["temperature"][int(indexq[0])])
			        rows_of_rep_clumps.append(row)
				print("come back and fix this, it happens")
					
          if len(temp_col_dens) != 0:
              log_sum_dens = np.log10(sum(temp_col_dens))
              col_density_split.append(log_sum_dens)
					
      ##make lists to put into dictionaries
      ray_nums.append(r)
      super_cl_nums.append(k)
      
      full_col_density = []
      make_full_list(col_density_match, full_col_density)
      make_full_list(col_density_split, full_col_density)
      make_full_list(col_density_short, full_col_density)
      med_col_dens.append(np.median(full_col_density))
      mad_for_med.append(stats.median_abs_deviation(full_col_density))
      
      num_clumps.append(len(full_col_density))
     
clump_stats = {}
clump_stats["ray_num"] = ray_nums
clump_stats["super_clump_number"] = super_cl_nums 
clump_stats["median_col_desnity"] = med_col_dens 
clump stats["mad_for_col_desnity"] = mad_for_med
clump_stats["distance_from_galaxy"] = distances 
clumps_stats["central_velocity"] = central_v 
clumps_stats["vel_dispersion"] = vel_dispersions 
clumps_stats["density"] = densities 
clumps_stats["temperature"] = temperatures
clumps_stats["num_of_clumps"] = num_clumps
clumps_stats["rep_clumps_row"] = rows_of_clump_reps

df = pd.DataFrame.from_dict(clump_stats)
df.to_csv("halo_rshift_ion_model-type_model-family_clump-type")
	  
    
    
