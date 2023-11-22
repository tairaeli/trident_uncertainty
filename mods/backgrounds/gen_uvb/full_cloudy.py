"""
Generates a UV Background using Cloudy Cooling Tools
"""
from mpi4py import MPI
import subprocess
import sys
import yt
yt.enable_parallelism()

def print_log(msg):
    """
    Writes a specified message to output file

    args:

        msg (string) - message to print to log file
    
    returns:

        None
    """
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    print ("[P%04d/%04d (%s)]: %s" % (rank, size, name, msg))

def run_command(command, timeout=None):
    """
    Runs a given command in parallel

    args:

        command (string) - path to file to command file that will be run in parallel

        timeout (float) - time before the run times out (seconds). 
                          Set to None by default
    
    returns:

        success (bool) - whether the run succeded or not
    """
    try:
        # where the magic happens, if all goes well, should have
        # a returncode of 0
        proc = subprocess.run(command, shell=True, timeout=timeout)
        if proc.returncode == 0:
            success = True
        else:
            success = False
    except subprocess.TimeoutExpired:
        print ("Process reached timeout of %d s. (%s)" % (timeout, command))
        success = False
    except KeyboardInterrupt:
        print ("Killed by keyboard interrupt!")
        success = False
    return success

if __name__ == "__main__":
    # how many processes to divide the run into
    start_part  = 1
    end_part    = 4
    total_parts = 4

    n_cores     = 128

    # location of parameter file
    par_file    = "./uvb_params.par"

    # CIAOLoop file location. Set to local installation of CIAOLoop
    CIAOLoop_file = "/mnt/home/tairaeli/astro_libs/cloudy_cooling_tools/CIAOLoop"

    for i in yt.parallel_objects(range(start_part, end_part+1), dynamic=True, njobs=-1):
        
        # should be set to local version of CIAO 
        my_command = CIAOLoop_file + " -rx -mp %d %d -np %d %s" % \
                     (i, total_parts, n_cores, par_file)
        print_log("Starting part %d of %d." % (i, total_parts))
        rval = run_command(my_command)
        if not rval:
            print_log("Part %d of %d failed." % (i, total_parts))
            sys.exit(1)

        print_log("Finished part %d of %d." % (i, total_parts))