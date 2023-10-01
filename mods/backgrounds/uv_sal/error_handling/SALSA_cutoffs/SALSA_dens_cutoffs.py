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
parser.add_argument('-iter_size', action='store', default = 0.8,
                    required=False, dest='cutoff_frac', type=float,
                    help='Fraction of mass above which a cutoff is made. \
                        Determines size of iterations of SPICE metod')

parser.add_argument('-min_dens', action='store', default = 13, 
                    required=False, dest='min_dens', type=float,
                    help='Fraction of mass above which a cutoff is made. \
                        Determines size of iterations of SPICE metod')

parser.add_argument('-uvb_path', action='store', 
                    required=False, dest='uvb', 
                    help='Path to UVB file')

parser.add_argument('-uvb_name', action='store', 
                    required=False, dest='uvb_name', 
                    help='Label to assign to uvb')

parser.add_argument('-make_plot', action='store', default=False,
                    required=False, dest='make_plot', type=bool,
                    help='Boolean deciding whether to make a plot or not.\
                        Should only be run if the 2 other uvbs have been run')

args = parser.parse_args()
dic_args = vars(args)

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
nrays = 1
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

# ray name not especially important here, but useful for when there are multiple rays
ray = 0

# defining data about the particular ascpects of the rays we are looking at
ion = "O VI"
atom,istate = ion.split(" ")

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

    save_dens = open(out_path+"density/"+args.uvb_name+"/"+args.uvb_name+"_dens_"+str(args.cutoff_frac)+".pickle","wb")
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

    """
    Set absorber min individually or globaly (this might not work)
    """

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

    save_clump = open(out_path+"clump_data/"+args.uvb_name+"/"+args.uvb_name+"_clump_dat_"+str(args.cutoff_frac)+".pickle","wb")
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
    plt.savefig(out_path+"min_dens_plots/UVB_dens_compare_"+str(args.min_dens)+".pdf")