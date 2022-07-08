import matplotlib.pyplot as plt
import pandas as pd

halo_marker_dict = {'4123' : 'x', '5016' : 'v', '5036' : '*', '8508' : 'D'}
halo_names_dict = {'2392'  :  'Hurricane' ,'2878'  :  'Cyclone' , '4123'  :  'Blizzard' , '5016'  :  'Squall' ,'5036'  :  'Maelstrom' , '8508'  :  'Tempest'}

def get_halo_names(num):
    if str(num) in halo_names_dict.keys():
        halo_name = halo_names_dict[str(num)]
    return halo_name


rs_lis = [2.0, 2.5] ##list of rs to analyze, true


ion_dict = {'C_II':'C4', 'C_IV':'C2', 'O_VI': 'C1'} ##dictionary of ions to colors

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
                plt.scatter(med_col_dens, col_dens_spread, c = ion_dict[ion], marker = halo_marker_dict[halo], label = f"{ion} in {name}")
      
    plt.title(f"MAD of Column Density vs Column Density, Redshift {rs}")
    plt.legend()
    plt.xlabel("Median log(Column Density)")
    plt.ylabel("MAD of log(Column Density)")
    plt.savefig(f'/mnt/scratch/f0104093/cgm_abundance_variance/graphs/mad_vs_med_z{rs}.png')
    plt.close()

print("I hope this works.")
