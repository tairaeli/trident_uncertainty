import pandas as pd
import argparse
import sys
import os

parser = argparse.ArgumentParser(description = "args for sorting script")
parser.add_argument('--ds', nargs='?', action='store', required=True, dest='path', help='Path to absorber catalogs. The program will create a directory here to store its generated txt files and png files')

args = parser.parse_args()
dic_args = vars(args)

path = os.path.expandvars(os.path.expanduser(args.path)) # path is actually a file, so need to split at filename and make 2 variables
split_path = os.path.split(path)
directory_path = split_path[0]
filename = split_path[1]
if not os.path.exists(directory_path):
    os.makedirs(dirctory_path)
    os.mkdir(directory_path+"/analysis")

data = pd.read_csv(path, delim_whitespace=True)

if filename[18] == "_":
    row_num = int(filename[17])
    element_name = ""
    for i in range(19, 23):
        element_name += filename[i]
else:
    row_num = int(filename[17]+filename[18])
    element_name = ""
    for i in range(20, 24):
       	element_name += filename[i]

if not os.path.exists(directory_path+f"/analysis/master_{element_name}.txt"):
    master = open(directory_path+f"/analysis/master_{element_name}.txt", "w")

master_dic = {}
for i in range(data['lightray_index'].max()):
    master_dic[f'ray{i}'] = []
for i in range(len(data['lightray_index'])):
    master_dic[f'ray{data['lightray_index'][i]}'].append(data['interval_start'][i])

master_df = pd.DataFrame.from_dict(master_dic)

master.write(master_df)

master.close() 
