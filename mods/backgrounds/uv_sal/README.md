# Directory Description

Code for analyzing how differences in Ultraviolet Background (UVB) models correspond to differences in the ion column densities found from SALSA.

### pipeline

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

### Current Plan
1. Fix this weird issue with SALSA, need to reach out to Devin about it
2. Afterwards I should reall reach out to the Trident Slack channel about being able to change the UVB model within the code to make my life easier
3. Then comes the plan of expanding this to ALL of the UVB models that are currently avaliable, I will probably need to streamline my process for running my code such that this process can be made easier
4. At some point, I would also like to be able to integrate Alexis' code (abundance table comparisons) into my own and allow us to see how UVB models and abudance table differences interact with one another May involve the scraping of a lot of my own code in favor of using Alexis' (as mine is not designed to handle so many diffferent UVB models)
