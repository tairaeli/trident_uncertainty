# Directory Description

Code for analyzing how differences in Ultraviolet Background (UVB) models correspond to differences in the ion column densities found from SALSA.

### pipeline

Where the magic happens. Contains all of the scripts necessary for running the analysis
- **sal_params.par**: Parameter file containing all of the arguments that feed into the code
- **ion_dat.txt**: Contains information on atomic mass, ionization energy, and estimated minimum observable CGM densities for each ion in analysis
- **absorber_extractor.py**: Runs main analysis code. Creates SALSA rays and calls on other scripts to perform analysis
- **absorber_extractor_H_I.py**: Performs the same functionality as *absorber_extractor.py*, but only for for H I
- **uvb_abun_pairwise_compare.py**: performs comparison between the SALSA data from the two different UVB models, then makes comparisons between the two sets of data to sort the data into different catagoies based on the differences/similarities of the models to one another 
- **condense_clumps.py**: Called in *absorber_extractor*. Takes the output of the pairwise comparison, then condenses the data into a format that can be used for more comprehensible comparisons
- **plot_total_column.py**: Makes plot of total column densities of selected ions for different ultraviolet backgrounds
- **plot_pair_comp.py**: Plots pairwise comparisons of column densities along with gas densities and temperatures
- **run_salsa_pipeline.sh**: bash script containing the order by which each script os run, along with the arguments included

### SALSA Cutoffs
Set of scripts for finding an optimal set of absorber cutoff fraction and minimum density settings for use in the absorber extraction methods. For a more thorough explanation of what these two quantities represent: https://salsa.readthedocs.io/en/latest/absorber_extraction.html#detailed-spice-method
- **SALSA_cutoff_finder.py**: Tests a range of settings in either absorber cutoff fraction or minimum density for a selected number of rays for a given UVB, then saves the resulting absorbers from each setting
- **SALSA_cutoff_plotter.py**: Creates histograms for settings tested in *SALSA_cutoff_finder.py*, both absorber cutoff fraction and minimum density
- **cutoff_check.sh**: bash script for running SALSA parameter analysis
