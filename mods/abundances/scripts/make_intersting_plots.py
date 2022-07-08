import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pandas as pd

halo_marker_dict = {'2392' : 'o', '4123' : 'x', '5016' : 'v', '5036' : '*', '8508' : 'D'}
halo_names_dict = {'2392'  :  'Hurricane' ,'2878'  :  'Cyclone' , '4123'  :  'Blizzard' , '5016'  :  'Squall' ,'5036'  :  'Maelstrom' , '8508'  :  'Tempest'}

def get_halo_names(num):
    if str(num) in halo_names_dict.keys():
        halo_name = halo_names_dict[str(num)]
    return halo_name


rs_lis = [2.0, 2.5] ##list of rs to analyze, true


ion_dict = {'C_II':'C4', 'C_IV':'C2', 'O_VI': 'C1'} ##dictionary of ions to colors

legend_elements = [Line2D([0], [0], lw= 0, color = 'k', marker = 'o', label = 'Hurricane'), Line2D([0], [0],  lw= 0, color = 'k', marker = 'x', label = 'Cyclone'), Line2D([0], [0],  lw= 0, color = 'k', marker = 'v', label = 'Squall'), Line2D([0], [0], lw= 0, color = 'k', marker = '*', label = 'Maelstrom'), Line2D([0], [0], lw= 0, color = 'k', marker = 'D', label = 'Tempest'), Patch(fc = 'C4', label = 'C_II'), Patch(fc = 'C2', label = 'C_IV'), Patch(fc = 'C1', label = 'O_IV')]

for rs in rs_lis:
    for halo in halo_marker_dict.keys():
        name = get_halo_names(halo)
        for ion in ion_dict.keys():
            if halo == '2392' and rs == 2.5:
                continue
            else:
                print(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv')
                ds = pd.read_csv(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv', delim_whitespace = True)
                med_col_dens = ds["median_col_desnity"]
                col_dens_spread = ds["mad_for_col_desnity"]
                plt.scatter(med_col_dens, col_dens_spread, c = ion_dict[ion], marker = halo_marker_dict[halo])
      
    plt.title(f"MAD of Column Density vs Column Density, Redshift {rs}")
    plt.legend(handles = legend_elements, ncol = 2)
    plt.xlabel("Median log(Column Density)")
    plt.ylabel("MAD of log(Column Density)")
    plt.savefig(f'/mnt/scratch/f0104093/cgm_abundance_variance/graphs/mad_vs_med_z{rs}.png')
    plt.close()

print("I hope this works.")

for rs in rs_lis:
    for halo in halo_marker_dict.keys():
        name = get_halo_names(halo)
        for ion in ion_dict.keys():
            if halo == '2392' and rs == 2.5:
                continue
            else:
                print(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv')
                ds = pd.read_csv(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv', delim_whitespace = True)
                med_col_dens = ds["median_col_desnity"]
                dens = ds["density"]
                plt.scatter(med_col_dens, dens, c = ion_dict[ion], marker = halo_marker_dict[halo])
      
    plt.title(f"Density vs Column Density, Redshift {rs}")
    plt.yscale('log')
    plt.legend(handles = legend_elements, ncol = 2)
    plt.xlabel("Median log(Column Density)")
    plt.ylabel("log(Density)")
    plt.savefig(f'/mnt/scratch/f0104093/cgm_abundance_variance/graphs/dens_vs_med_z{rs}.png')
    plt.close()
    
    
for rs in rs_lis:
    for halo in halo_marker_dict.keys():
        name = get_halo_names(halo)
        for ion in ion_dict.keys():
            if halo == '2392' and rs == 2.5:
                continue
            else:
                print(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv')
                ds = pd.read_csv(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv', delim_whitespace = True)
                med_col_dens = ds["median_col_desnity"]
                temp = ds["temperature"]
                plt.scatter(med_col_dens, temp, c = ion_dict[ion], marker = halo_marker_dict[halo])
      
    plt.title(f"Temperature vs Column Density, Redshift {rs}")
    plt.yscale('log')
    plt.legend(handles = legend_elements, ncol = 2)
    plt.xlabel("Median log(Column Density)")
    plt.ylabel("log(Temperature)")
    plt.savefig(f'/mnt/scratch/f0104093/cgm_abundance_variance/graphs/temp_vs_med_z{rs}.png')
    plt.close()
    
for rs in rs_lis:
    for halo in halo_marker_dict.keys():
        name = get_halo_names(halo)
        for ion in ion_dict.keys():
            if halo == '2392' and rs == 2.5:
                continue
            else:
                print(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv')
                ds = pd.read_csv(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv', delim_whitespace = True)
                med_col_dens = ds["median_col_desnity"]
                dist = ds["distance_from_galaxy"]
                plt.scatter(med_col_dens, dist, c = ion_dict[ion], marker = halo_marker_dict[halo])
      
    plt.title(f"Distance vs Column Density, Redshift {rs}")
    plt.legend(handles = legend_elements, ncol = 2)
    plt.xlabel("Median log(Column Density)")
    plt.ylabel("Distance from Galaxy")
    plt.savefig(f'/mnt/scratch/f0104093/cgm_abundance_variance/graphs/dist_vs_med_z{rs}.png')
    plt.close()
    
    
for rs in rs_lis:
    for halo in halo_marker_dict.keys():
        name = get_halo_names(halo)
        for ion in ion_dict.keys():
            if halo == '2392' and rs == 2.5:
                continue
            else:
                print(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv')
                ds = pd.read_csv(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv', delim_whitespace = True)
                col_dens_spread = ds["mad_for_col_desnity"]
                dens = ds["density"]
                plt.scatter(dens, col_dens_spread, c = ion_dict[ion], marker = halo_marker_dict[halo])
      
    plt.title(f"Density vs Spread of Column Density , Redshift {rs}")
    plt.xscale('log')
    plt.legend(handles = legend_elements, ncol = 2)
    plt.ylabel("MAD of Median log(Column Density)")
    plt.xlabel("log(Density)")
    plt.savefig(f'/mnt/scratch/f0104093/cgm_abundance_variance/graphs/dens_vs_mad_z{rs}.png')
    plt.close()
    
for rs in rs_lis:
    for halo in halo_marker_dict.keys():
        name = get_halo_names(halo)
        for ion in ion_dict.keys():
            if halo == '2392' and rs == 2.5:
                continue
            else:
                print(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv')
                ds = pd.read_csv(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv', delim_whitespace = True)
                col_dens_spread = ds["mad_for_col_desnity"]
                temp = ds["temperature"]
                plt.scatter(temp, col_dens_spread, c = ion_dict[ion], marker = halo_marker_dict[halo])
      
    plt.title(f"Temperature vs Spread of Column Density, Redshift {rs}")
    plt.xscale('log')
    plt.legend(handles = legend_elements, ncol = 2)
    plt.ylabel("MAD of Median log(Column Density)")
    plt.xlabel("log(Temperature)")
    plt.savefig(f'/mnt/scratch/f0104093/cgm_abundance_variance/graphs/temp_vs_mad_z{rs}.png')
    plt.close()
    
for rs in rs_lis:
    for halo in halo_marker_dict.keys():
        name = get_halo_names(halo)
        for ion in ion_dict.keys():
            if halo == '2392' and rs == 2.5:
                continue
            else:
                print(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv')
                ds = pd.read_csv(f'/mnt/scratch/f0104093/cgm_abundance_variance/halo{halo}/redshift{rs}/stats/{halo}_z{rs}_{ion}_abun_all-model-families_all-clumps.csv', delim_whitespace = True)
                col_dens_spread = ds["mad_for_col_desnity"]
                dist = ds["distance_from_galaxy"]
                plt.scatter(dist, col_dens_spread, c = ion_dict[ion], marker = halo_marker_dict[halo])
      
    plt.title(f"Distance vs Spread of Column Density, Redshift {rs}")
    plt.legend(handles = legend_elements, ncol = 2)
    plt.ylabel("MAD of Median log(Column Density)")
    plt.xlabel("Distance from Galaxy")
    plt.savefig(f'/mnt/scratch/f0104093/cgm_abundance_variance/graphs/dist_vs_mad_z{rs}.png')
    plt.close()
