# In this repository

Contains 2 projects on analysis of a simulated CGM based on the [FOGGIE galaxies](https://foggie.science/). Both of these can be found within the *mods* directory. These projects utilize a Synthetic Absorption Line Surveyor Application ([SALSA](https://salsa.readthedocs.io/en/latest/index.html)), a code build on Trident and yt to analyze dense clouds of gas known as absorbers along many randomly oriented lines of sight across the CGM.

## Element Abundance Analysis
The first project, under the directory *abundances* focuses on the abundance of different elements within the CGM. It tests a suite of different abundance patterns, running them through a SALSA pipeline to analyze the distribution of absorbers and how it varies with different patterns.

## Ultraviolet Background
The next project, under *backgrounds* focuses on the ultraviolet background (UVB) and how different models result in different absorber column densities within the CGM. This project includes code for running pairwise analysis on different UVBs and generating a suite of figures to display the results. However, for the pairwise analysis to work, the UVBs must be included as an ionization table. As such, this directory includes code for running [CLOUDY](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home) to generate these tables. The full paper containing the first analysis performed with this work can be found [here](https://arxiv.org/abs/2503.11775)
