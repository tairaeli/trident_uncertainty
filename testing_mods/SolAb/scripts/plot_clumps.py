import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import argparse
# %matplotlib inline
import sys
import os

parser = argparse.ArgumentParser(description = "args for sorting script")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path to absorber catalogs. The program will create a directory here to store its generated txt files and png files')
parser.add_argument('--fn', nargs='?', action='store', required=True, dest='filename_list', help='txt file of all catalogs corresponding to a given ion')
args = parser.parse_args()
dic_args = vars(args)

directory_path = os.path.expandvars(os.path.expanduser(args.path))

if not os.path.exists(directory_path):
    print("MAKING CLUMP_PLOTS DIRECTORY...")
    os.makedirs(dirctory_path)
    os.mkdir(directory_path+"clump_plots")

with open(args.filename_list) as f:
    files = f.read().splitlines()
# files = args.filename_list
preliminary_dummy_data = pd.read_csv(files[0], delim_whitespace=True)
ion_position = files[0].find('.')-4
ion = files[0][ion_position]
for i in range(ion_position+1, ion_position+4):
    ion += files[0][i]

tick_range = int(round(len(preliminary_dummy_data)/2, 0))
num_ticks = list(range(-1*(tick_range+1), tick_range+1, 1))
tick_labels = [' ']
# print(f"MAX RAY NUM: {preliminary_dummy_data['lightray_index'].max()}")
for i in range(1, len(num_ticks)):
    tick_labels.append(f'row {i}')
for ray_num in range(int(preliminary_dummy_data['lightray_index'].max())):
    for filename in files:
        data = pd.read_csv(filename, delim_whitespace=True)
        ray_data = data[data['lightray_index']==ray_num]
        plt.hlines(np.ones(ray_data.shape[0])*ray_num, ray_data['interval_start'], ray_data['interval_end'])
        plt.title(f"{ion} -- Cell Index {ray_num}")
        plt.xlabel("Cell Index")
        plt.yticks(num_ticks, tick_labels)
    plt.savefig(directory_path+f"clump_plots/ClumpPlot_{ion}_RayIndex{ray_num}.png")