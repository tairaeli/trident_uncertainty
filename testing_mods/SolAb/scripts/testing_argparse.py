import argparse
import sys
parser = argparse.ArgumentParser(description = "Preliminary constants for SALSA pipeline.")
parser.add_argument('-ds', nargs='?', action='store', const='/mnt/research/fuhrmane/test_sal/', dest='path', help='Path where rays and output data will be stored. Directory should contain three other directories called "data", "rays", and "visuals" for program to run smoothly.')
parser.add_argument('--nrays', action='store', dest='nrays', default=4, type=int, help='The number of rays to be generated.')
#this is the one where it would be handy if you could interact further with the user append the middle of the string so just providinig a halo
parser.add_argument('--halo', action='store', dest='ds_file', default='/mnt/research/galaxies-REU/sims/FOGGIE/halo_002392/nref11c_nref9f/RD0020/RD0020', help='Path to halo data.')
#would also be handy if could choose between appending to existing directory path or just providing their own...same with data_storage
parser.add_argument('--abun', action='store', dest='file_path', default=argparse.SUPPRESS, help='Path to abundance file, if any. Defaults to solar abundances.')
#args = parser.parse_args()
args = parser.parse_args()
dic_args = vars(args)

#confirm everything is correct with user before proceding
if 'file_path' in dic_args:
    abundances = args.file_path
else:
    abundances = 'None. Using SolAb data.'

cont = input(f"RAY AND OUTPUT DATA WILL BE STORED HERE: {args.path}\nNUMBER OF RAYS TOO BE GENERATED: {args.nrays}\nUSING HALO DATA CONTAINED HERE: {args.ds_file}\nABUNDANCE DATA TO BE USED: {abundances}\n\ncontinue?[y/n] ")
if cont == 'y':
    print("yeet")

else:
    sys.exit("Pull your shit together, Becky.")

# test whether sys.exit() worked -- this shuoldn't print
print('ya fucked up')
