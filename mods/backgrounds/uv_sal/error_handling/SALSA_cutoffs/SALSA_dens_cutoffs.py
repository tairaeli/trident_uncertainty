# importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as sp
import pickle
from mpi4py import MPI
import yt
import trident
import salsa
from salsa.utils import check_rays
import argparse

print(f"Let's do some math, kids")

parser = argparse.ArgumentParser(description = "Select cutoff and UVB for analysis")
parser.add_argument('-cutoff_frac', action='store', default = 0.8,
                    required=False, dest='cutoff_frac', type=float,
                    help='Fraction of mass above which a cutoff is made. \
                        Determines size of iterations of SPICE metod.')

parser.add_argument('-min_dens', action='store', default = 13, 
                    required=False, dest='min_dens', type=float,
                    help='Minimum density at which SALSA method stops\
                        iterating.')

parser.add_argument("-var_arg", action='store',
                    required=True, dest="var_arg",
                    help="argument that is varied between each iteration.\
                        Can either be 'cutoff_frac' or 'min_dens'.")

parser.add_argument('-ion_list', action='store', 
                    required=False, dest='ion_list', 
                    help='Ions to be analyzed')

parser.add_argument('-nrays', action='store', 
                    required=False, dest='nrays', type = int,
                    help='Number of rays to be used for analysis')

parser.add_argument('-uvb_path', action='store', 
                    required=False, dest='uvb', 
                    help='Path to UVB file')

parser.add_argument('-uvb_name', action='store', 
                    required=False, dest='uvb_name', 
                    help='Label to assign to uvb.')

parser.add_argument('-make_plot', action='store', default=False,
                    required=False, dest='make_plot', type=bool,
                    help='Boolean deciding whether to make a plot or not.\
                        Should only be run if the 2 other uvbs have been run.')

args = parser.parse_args()
dic_args = vars(args)

# variable for determining labels of output files
output_val = None

# determining output_val based on input arguments
if args.var_arg == "cutoff_frac":
    # checking to see if input value is valid
    if args.cutoff_frac > 1.0:
        raise Exception("Invalid value of cutoff_frac. Must be value \
                        in between 0 and 1")
    else:
        output_val = str(args.cutoff_frac)

elif args.var_arg == "min_dens":
    output_val = str(args.min_dens)

else:
    raise Exception("Invalid input of var_arg. Input must either be \
                    'cutoff_frac' or 'min_dens'")

# specifying which halo we are looking at and at what redshift
halo = 2392
rs = 20

halo_path = f"/mnt/research/galaxies-REU/sims/FOGGIE/halo_00{halo}/nref11c_nref9f/RD00{rs}/RD00{rs}"

halo_data = yt.load(halo_path)

# EDIT THIS CELL 
foggie_dir = "/mnt/home/tairaeli/astro_libs/foggie/foggie/halo_infos"

center_dat = {}

# loading in FOGGIE data
raw_foggie_dat = pd.read_csv(f"{foggie_dir}/00{halo}/nref11c_nref9f/halo_c_v", sep = '|', 
                             names = ['null','redshift','name','xc','yc','zc','xv','yv','zv','null2'])

# making some fixes specific to these files
raw_foggie_dat = raw_foggie_dat.drop(0)
raw_foggie_dat = raw_foggie_dat.drop(columns = ['null','null2'])

# isolating data to a specific redshift 
raw_foggie_dat = raw_foggie_dat[raw_foggie_dat['name'] == ' RD00'+str(rs)+' ']

# storing the position and velocity data of the galactic center
center_dat['pos'] = [float(raw_foggie_dat["xc"]),float(raw_foggie_dat["yc"]),float(raw_foggie_dat["zc"])]
center_dat['vel'] = [float(raw_foggie_dat["xv"]),float(raw_foggie_dat["yv"]),float(raw_foggie_dat["zv"])]

# input fields as they appear in parameter file
field_items = "density temperature metallicity radius".split(" ")
field_types = "gas NA NA index".split(" ")
units = "g/cm**3 NA Zsun NA".split(" ")

# converting parameters from a string to a dictionary to make SALSA happy
field_units = {}
other_fields = []
for i, item in enumerate(field_items):
    
    if field_types[i] == "NA":
        field_tup = item
    else:
        field_tup = (field_types[i],item)
        
    other_fields.append(field_tup)
    
    if units[i] == "NA":
        continue
    else:
        field_units[f"{item}"] = units[i]


# defining features about calaxy center
center = halo_data.arr(center_dat['pos'], 'kpc')
gal_vel = halo_data.arr(center_dat['vel'], 'km/s')

# defining characteristics of ray behavior
nrays = args.nrays
max_impact = 15

# out_path = "/mnt/home/tairaeli/trident_uncertainty/mods/backgrounds/uv_sal/error_handling/SALSA_cutoffs/cutoff_bin/"

out_path = "/mnt/scratch/tairaeli/cutoff_bin/"

# defining random seed to allow for repitition
np.random.seed(13)

# creating rays and saving results

check = check_rays("./", nrays, [])
if not check:
    print("WARNING: rays not found. Generating new ones.")
    salsa.generate_lrays(halo_data, 
                     center.to('code_length'), 
                     nrays, 
                     max_impact, 
                     length=600, 
                     field_parameters={'bulk_velocity':gal_vel}, 
                     ion_list=['H I'], 
                     fields=other_fields, 
                     out_dir="./")

# initializing ion list
init_ion_list = args.ion_list.split(" ")

# adjusting formating of ion list
ion_list = []
for i, ion in enumerate(init_ion_list):
    
    alt_ion = ion.replace("_"," ")
    ion_list.append(alt_ion)

# iterating through each ray
for ray in range(nrays):
    
    # iterating through each ion
    for ion in ion_list:
        
        atom, istate = ion.split(" ")

        field_name = f"{atom}_p{trident.from_roman(istate)-1}_number_density"

        # defines abundance table used
        # CHANGE THESE TO LOCAL VERSIONS
        abun_table_path = "/mnt/home/tairaeli/trident_uncertainty/mods/abundances/scripts/abun_table.txt"
        abun = pd.read_csv(abun_table_path, delim_whitespace=True)
        abundances = abun.iloc[0].to_dict()

        ray_filename = f"ray{ray}.h5"

        ray_dat = yt.load(ray_filename)

        if not args.make_plot:
            trident.add_ion_number_density_field(atom, trident.from_roman(istate), 
                                                ray_dat, abundance_dict = abundances, 
                                                ionization_table = args.uvb)

            # extracting data from ray analysis
            gas_dens = ray_dat.r[("gas",field_name)].copy()

            save_dens = open(out_path+"density/"+args.uvb_name+"/"+args.uvb_name+"_dens_"+args.var_arg+"_"+output_val+".pickle","wb")
            pickle.dump(gas_dens, save_dens, protocol=3)
            save_dens.close()

            # a bit more prep before the absorber extraction step
            comm = MPI.COMM_WORLD

            ray_list=[ray_filename]

            # CK: Taking a hint from SALSA on how to divvy up the ray list across procs
            # works under assumption we have multiple rays, but in this case we do not
            ray_arr = np.array(ray_list)
            ray_files_split = np.array_split(ray_arr, comm.size)
            my_rays = ray_files_split[comm.rank]

            abs_ext = salsa.AbsorberExtractor(halo_data, 
                                            ray_filename, 
                                            ion_name = ion, 
                                            velocity_res =20, 
                                            abundance_table = abundances, 
                                            calc_missing=True,
                                            frac = args.cutoff_frac,
                                            absorber_min = args.min_dens)

            # mimicing how data is stored in code
            clump_dat = {ion:{args.uvb_name:None}}
            
            clump_dat[ion][args.uvb_name] = salsa.get_absorbers(abs_ext, 
                                            my_rays, 
                                            method='spice', 
                                            fields=other_fields, 
                                            units_dict=field_units)
            
            clump_dat[ion][args.uvb_name] = clump_dat[ion][args.uvb_name].drop(columns='index')

            save_clump = open(out_path+"clump_data/"+args.uvb_name+"/"+args.uvb_name+"_clump_dat_"+args.var_arg+"_"+output_val+".pickle","wb")
            pickle.dump(clump_dat, save_clump, protocol=3)
            save_clump.close()

        else:

            with open(out_path+"density/HM_2012/HM_2012_dens_"+str(args.cutoff_frac)+".pickle", "rb") as dens_dat:
                hm_dens = pickle.load(dens_dat)

            with open(out_path+"clump_data/HM_2012/HM_2012_clump_dat_"+str(args.cutoff_frac)+".pickle", "rb") as salsa_dat:
                hm_clump_dat = pickle.load(salsa_dat)

            with open(out_path+"density/PCW_2019/PCW_2019_dens_"+str(args.cutoff_frac)+".pickle", "rb") as dens_dat:
                pcw_dens = pickle.load(dens_dat)

            with open(out_path+"clump_data/PCW_2019/PCW_2019_clump_dat_"+str(args.cutoff_frac)+".pickle", "rb") as salsa_dat:
                pcw_clump_dat = pickle.load(salsa_dat)

            colors = plt.cm.cool([0.3,0.4,0.9,1])

            ray_pos = ray_dat.r[("gas","l")].to("kpc")

            print("Data loaded, generating plot...")

            # creating density curves for both UVB models
            plt.figure(figsize = [15,8], dpi = 500, facecolor = "white")
            plt.semilogy(ray_pos, hm_dens,label = "HM 2012", color = colors[3])
            plt.semilogy(ray_pos, pcw_dens,label = "PCW 2019", color = colors[1])
            # plt.axhline(1e-13)

            hm_clumps = hm_clump_dat[ion]["HM_2012"][hm_clump_dat[ion]["HM_2012"]["lightray_index"] == str(ray)]
            pcw_clumps = pcw_clump_dat[ion]["PCW_2019"][pcw_clump_dat[ion]["PCW_2019"]["lightray_index"] == str(ray)]

            # showing where SALSA detected clumps in Haart & Madau
            i = 0
            for i in range(len(hm_clumps["interval_start"])):
                
                lb = hm_clumps["interval_start"][i]
                hb = hm_clumps["interval_end"][i]
                rng = [lb,hb]
                
                yb = hm_dens[slice(*rng)]
                xb = ray_pos[lb:hb]
                
                plt.fill_between(xb,yb, color = colors[2], alpha = 0.3)

            # showing where SALSA detected clumps in Puchwein et al
            for i in range(len(pcw_clumps["interval_start"])):
                
                lb = pcw_clumps["interval_start"][i]
                hb = pcw_clumps["interval_end"][i]
                
                rng = [lb,hb]
                yb = pcw_dens[slice(*rng)]
                xb = ray_pos[lb:hb]
                
                plt.fill_between(xb,yb, color = colors[0], alpha = 0.3)    

            plt.grid()
            plt.legend(fontsize = 20)
            # plt.ylim(10**(-15),10**(-6.5))
            # plt.xlim(500,900)
            plt.xlabel("Position Along Ray (kpc)", fontsize = 25)
            plt.ylabel(r"Density ($cm^{-2}$)", fontsize = 25)
            # plt.ylim(1e-12, 2e-6)
            plt.title(f"UVB Number Density Comparison "+ion, fontsize = 30)
            plt.savefig(out_path+args.var_arg+"_plots/UVB_dens_compare_"+args.var_arg+"_"+output_val+".pdf")