# In This Directory

**gen_cloudy_input.py** - takes in information from UV backgrounds and combines it with UV emissions from a galaxy. Uses results chemical evolution code [Omega+](https://ui.adsabs.harvard.edu/abs/2017ApJ...835..128C/abstract) to determine star formation rate of galaxy, then uses FSPS to determine UVB emissions. Then, the UV data is reformatted to match [CLOUDY]() input requirements for ionization table generation. Takes 5 arguments:
* --ds (str): location of output file
* --om (str): path to omega+ output data
* --uvb (str): location of UVB data
* --rs (str): set range for redshifts to iterate through. Should be formatted with a comma separating the two redshifts (e.g. "1.2,2.7").
* --d (list): integer list of distances (in units of kpc) from the galaxy

**gen_cloudy_input_no_gal.py** - performs the same function as *gen_cloudyInput* without including the contribution of the galactic UV emission. Takes in 3 arguments:
* --ds (str): location of output file
* --uvb (str): location of UVB data
* --rs (str): set range for redshifts to iterate through (e.g. "1.2,2.7")

**full_cloudy.py** - runs cloudy spectral synthesis simulations in parallel. Requires 6 arguments that are defined within the code:
* start_part - first part number in which the code will run from. The cloudy simulation is broken up into a number of parts to be run in parallel.
* end_part - final part number in which the code will run
* total_parts - total number of parts in cloudy run
* number of cores - number of gpu cores the code is run on
* par_file - location of parameter file (more details below)
* CIAOLoop_file - location of CIAOLoop file. Is where the spectral synthesis is run. Should be located within installation of cloudy_cooling_tools

**uvb_params.par** - parameter file for cloudy simulation run. Includes all of the parameters that go into the cloudy run that are specified within the file. Here, there are two files of note that are required for cloudy to run:
* path to cloudy executable
* path to ultraviolet background data

**backup_jobs** - directory containing job scripts for running cloudy simulations in several parts. Should only be used if 1-2 parts take significantly longer to converge than the rest. Otherwise, 'full_cloudy.py' should be run instead.

**combine_runfile_parts.pl** - if job scripts from 'backup_jobs' directory are run, combines output files of each job into a single run file. Note that this file is written in Perl, so it will need to be run as such.

**gen_trident_input.py** - takes in final run file (either from the output of full_cloudy.py or the output of combine_runfile_parts.pl) and outputs an ionization table in an h5 format that can be used in the trident.
