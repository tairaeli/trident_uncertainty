# In This Directory

**full_cloudy.py** - runs cloudy simulations in parallel. Requires 6 arguments:
* start_part - first part number in which the code will run in parallel
* end_part - final part number in which the code will run in parallel
* total_parts - total number of parts in cloudy run
* number of cores - number of gpu cores the code is run on
* par_file - location of parameter file (more details below)
* CIAOLoop_file - location of CIAOLoop file. Should be located within installation of cloudy_cooling_tools

**uvb_params.par** - parameter file for cloudy simulation run. Includes all of the parameters that go into the cloudy run that are specified within the file. Here, there are tow files of note that are required for cloudy to run:
* path to cloudy executable
* path to ultraviolet background data

**backup_jobs** - directory containing job scripts for running cloudy simulations in several parts. Should only be used if 1-2 parts take significantly longer to converge than the rest. Otherwise, 'full_cloudy.py' should be run instead.

**combine_runfile_parts.pl** - if job scripts from 'backup_jobs' directory are run, combines output files of each job into a single run file. Note that this file is written in Perl, so it will need to be run as such.

**gen_trident_input.py** - takes in final run file (either from the output of full_cloudy.py or the output of combine_runfile_parts.pl) and outputs an ionization table in h5 that can be used in the trident pipeline.