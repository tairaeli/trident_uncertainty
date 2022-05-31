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
            
##next step is to make histograms for column density for each clump we have data for
