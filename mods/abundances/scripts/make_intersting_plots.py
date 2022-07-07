import matplotlib.pyplot as plt
import pandas as pd

halo_lis = [] ##list of halos to analyze
rs_lis = [] ##list of rs to analyze, true
ion_lis = [] ##list of ions to analyze

for halo in halo_lis:
  ##get names??
  for rs in rs_lis:
    for ion in ion_lis:
      ds = pd.read_csv(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv)
      med_col_dens = ds["median_col_density"]
      col_dens_spread = ds["mad_for_col_density"]
      density = ds["density"]
      temperature = ds["temperature"]
      distance = ds["distance_from_galaxy"]
