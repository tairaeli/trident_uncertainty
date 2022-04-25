import pandas as pd
import matplotlib.pylot as plt
import os
import sys
import argparse

parser = argparse.ArgumentParser(description = "args for sorting script")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path to absorber catalogs. The program will create a directory here to store its generated txt files and png files')
parser.add_argument('--il', nargs='?', action='store', required=True, dest='interval_list', help='2D list of interval ranges')
parser.add_argument('--fn', nargs='?', action='store', required=True, dest='filename', help='txt file of clump data extacted in plot_clumps.py')
args = parser.parse_args()
dic_args = vars(args)

directory_path = os.path.expandvars(os.path.expanduser(args.path))
clump_dir = directory_path+"clump_plots3"

data = pd.read_csv(args.filename, delim_whitespace=True)
ion_position = args.filename.find('.')-4
ion = args.filename[ion_position]
for i in range(ion_position+1, ion_position+4):
    ion += agrs.filename[i]

for row_num in range(data['row_number'].max()):
    for ray_num in range(data['cell_index'].max()):
        ray_and_row_data = data[data['cell_index']==ray_num & data['row_number']==row_num]
        plt.hist(ray_and_row_data, bins=interval_list, stacked=True)
        Title = f'{ion} -- Cell Index {ray_num}'
    
    plt.title(Title)
    plt.xlabel('interval ranges')
    plt.ylabel('column densities')
        # for  i in range(len(ray_and_row_data)):
        #     for int in args.interval_list:
        #         if (ray_and_row_data['interval_start'][i] >= int[0] & ray_and_row_data['interval_end'] <= int[1]):
        #             # add to bin