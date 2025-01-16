# Directory Description

Where the magic happens. Contains all of the scripts necessary for running the analysis
- **sal_params.par**: Parameter file containing all of the arguments that feed into the code
- **ion_dat.txt**: Contains information on atomic mass, ionization energy, and estimated minimum observable CGM densities for each ion in analysis
- **sal_the_super_uvb.py**: Runs main analysis code. Creates SALSA rays and calls on other scripts to perform analysis
- **sal_the_super_uvb_H_I.py**: Performs the same functionality as *sal_the_super_uvb.py*, but only for for H I
- **uvb_abun_pairwise_compare.py**: performs comparison between the SALSA data from the two different UVB models, then makes comparisons between the two sets of data to sort the data into different catagoies based on the differences/similarities of the models to one another 
- **condense_clumps.py**: Called in *sal_the_super_uvb*. Takes the output of the pairwise comparison, then condenses the data into a format that can be used for more comprehensible comparisons
- **plot_total_column.py**: Makes plot of total column densities of selected ions for different ultraviolet backgrounds
- **plot_pair_comp.py**: Plots pairwise comparisons of column densities along with gas densities and temperatures
- **run_salsa_pipeline.sh**: bash script containing the order by which each script os run, along with the arguments included