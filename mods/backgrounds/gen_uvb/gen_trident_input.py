# Converts CIAOLoop output to a format that can be used in trident
import sys

# setting path to local installation of cloudy grids within cloudy cooling tools
cloudy_grid_path = "/mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/cloudy_grids"

# adding path to system
sys.path.append(cloudy_grid_path)

# importing cloudy_grid functions
from cloudy_grids import convert_ion_balance_tables

# converting '.run' files from useable '.h5' files
convert_ion_balance_tables(
    "/mnt/scratch/tairaeli/pcw_ionization/test1.run", # path to file created from runfile combination
    "/mnt/scratch/tairaeli/trident_inputs/pcw_test_12_4_24.h5", # path to output directory
    ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", # list of ions
     "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", 
     "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"])
    