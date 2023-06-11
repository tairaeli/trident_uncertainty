# Directory Description

Code for analyzing how differences in Ultraviolet Background (UVB) models correspond to differences in the ion column densities found from SALSA.

### pipeline

Where the magic happens. Contains all of the scripts necessary for running the analysis
- **sal_params.par**: Parameter file containing all of the arguments that feed into the code
- **sal_the_super_uvb.py**: Runs main analysis code. Creates SALSA rays and calls on other scripts to perform analysis
- **uvb_abun_pairwise_compare.py**: performs comparison between the SALSA data from the two different UVB models, then makes comparisons between the two sets of data to sort the data into different catagoies based on the differences/similarities of the models to one another 
- **condense_clumps.py**: Called in *sal_the_super_uvb*. Takes the output of the pairwise comparison, then condenses the data into a format that can be used for more comprehensible comparisons

### notebooks

Contains the notebooks where the majority of the testing of the pipeline is done
- **sal_module_testing.ipynb**: evaluates individual parts on the pipeline to ensure that they are functioning correctly (kind of a mess of at the moment)
- **uvb_col_density_compare.ipynb**: used for creating visualizations of the output from the pipeline data
- **salsa_issues**: directory where I attempt to highlight the issues I've been running into with SALSA

### Current Plan
1. Fix this weird issue with SALSA, need to reach out to Devin about it
2. Afterwards I should reall reach out to the Trident Slack channel about being able to change the UVB model within the code to make my life easier
3. Then comes the plan of expanding this to ALL of the UVB models that are currently avaliable, I will probably need to streamline my process for running my code such that this process can be made easier
4. At some point, I would also like to be able to integrate Alexis' code (abundance table comparisons) into my own and allow us to see how UVB models and abudance table differences interact with one another May involve the scraping of a lot of my own code in favor of using Alexis' (as mine is not designed to handle so many diffferent UVB models)