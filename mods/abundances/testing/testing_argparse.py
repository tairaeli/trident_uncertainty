import argparse
import sys
parser = argparse.ArgumentParser(description = "Preliminary constants for SALSA pipeline.")
parser.add_argument('-ds', nargs='?', action='store', required=True, dest='path', help='Path where rays and output data will be stored. Directory should contain three other directories called "data", "rays", and "visuals" for program to run smoothly.')
parser.add_argument('--nrays', action='store', dest='nrays', default=4, type=int, help='The number of rays to be generated.')
parser.add_argument('--abun', action='store', dest='file_path', default=argparse.SUPPRESS, help='Path to abundance file, if any. Defaults to solar abundances.')
# ask claire about how to use nsc flag without interaction (want to specify just the halo)
# parser.add_argument('-nsc', action='store_true', dest='nsc', help='nsc stands for "not super computer". Use this flag if you would like to manually type out the entire path to your halo data. Otherwise, the interpreter assumes you are using the HPCC FOGGIE simulation data.')
parser.add_argument('--halo', action='store', dest='ds_file', default='/mnt/research/galaxies-REU/sims/FOGGIE/halo_002392/nref11c_nref9f/RD0020/RD0020', help='Path to halo data.')

#args = parser.parse_args()
args = parser.parse_args()
dic_args = vars(args)

# define some variables 
path_list = ['/mnt/research/fuhrmane/test_sal/', args.path]
path = ''.join(path_list)
if 'file_path' in dic_args:
    abundances = args.file_path
else:
    abundances = 'No file given. Using solar abundances.'

cont = input(f"RAY AND OUTPUT DATA WILL BE STORED HERE: {path}\nNUMBER OF RAYS TOO BE GENERATED: {args.nrays}\nUSING HALO DATA CONTAINED HERE: {args.ds_file}\nABUNDANCE DATA TO BE USED: {abundances}\n\ncontinue?[y/n] ")
if cont == 'y':
    print("yeet")

else:
    sys.exit("Pull your shit together, Becky.")

# test whether sys.exit() worked -- this shuoldn't print