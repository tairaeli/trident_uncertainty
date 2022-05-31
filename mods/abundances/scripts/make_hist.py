import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

# Takes path to where SALSA abunance catalogs are stored and ion list used in generating them to extract necessary columnn density data into a big boi dictionary
path = "/mnt/scratch/f0104093/condensed_pipeline_tests/data/" # modify as needed
row_dictionary = {}
coldens_dic = {}
ion_list = ["C_II", "O_VI", "C_IV"] # modify as needed
for ion in ion_list:
    row_dictionary[f'{ion}'] = {}
    coldens_dic[f'{ion}'] = {}
    for i in range(26):
        m = i+1
        row_dictionary[f'{ion}'][f'row{m}'] = {}
		n_len = len(str(m))
		n_zeros = ndigits - n_len
		k = "0" * n_zeros + str(m)
        row_data = pd.read_csv(path+f"data_AbundanceRow{k}_{ion}.txt", delim_whitespace=True)
        for ray in range(4):
            row_dictionary[f'{ion}'][f'row{i+1}'][f'ray_index{ray}'] = row_data[row_data['lightray_index']==ray].reset_index()
            # row_dictionary[f'{ion}'][f'row{i+1}'][f'ray_index{ray}'].reset_index()
            coldens_dic[f'{ion}'][f'column_densities{ray}'] = np.zeros(26)

# Somewhat automatically -- when inteerval chunks line up exactly. Use the "ray", "ion", and "interval" variables instead of modifying the for loop
ray = "1"
ion = ion_list[1]
interval = 0
for i in range(26):
    row = row_dictionary[ion][f'row{i+1}']
    coldens_dic[ion][f'column_densities{ray}'][i] += row_dictionary[ion][f'row{i+1}'][f'ray_index{ray}']["col_dens"][interval]
plt.hist(coldens_dic[ion][f'column_densities{ray}'])
plt.title(f"{ion} -- LightRay Index {ray} -- Interval {str(interval)}")
plt.savefig(f"/mnt/scratch/f0104093/condensed_pipeline_tests/data/Hist_{ion}_RayIndex{ray}_Int{interval}.png")

# Manually -- for when interval chunks for rows don't exactly line up
coldens_dic[ion][f'column_densities{ray}'][0] += row_dictionary[ion][f'row1'][f'ray_index{ray}']["col_dens"][0]

plt.hist(coldens_dic[ion][f'column_densities{ray}'])
plt.title(f"{ion} -- LightRay Index {ray} -- Interval {str(interval)}")
plt.savefig(f"/mnt/scratch/f0104093/condensed_pipeline_tests/data/Hist_{ion}}_RayIndex{ray}_Int{str(interval)}.png")
