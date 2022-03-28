# trident_mods
in progress modifications to Trident

A Pipeline designed to generate SALSA absorber catalogs from Trident for the purpose of quantifying uncertainties inherent to synthetic obsrvations. 

# Scripts
* sal_the_snek.py - uses SALSA to generate and save lightray objects and absorber catalogs to a given directory. Accepts a txt or csv file with abundance patterns different from solar abundances. If none are given, the saved data will be from that of solar abundances. The user can specify the number of lightray objects to be generated, the path to the directory where the data and rays will be stored, the path to an optional abundance file, and the path to the halo where the lightray objects will "pass through". This file can be run directly from the command line or though a bash script like sal_bash.txt
* sal_bash.txt - generally used for testing sal_the_snek.py functionality, and can serve as an example of how to use sal_the_snek.py

# Command Line Arguments
* --ds: Path where rays and output data will be stored. Directory should contain two other directories called "data" and "rays",for program to run smoothly (if neither are there, the program will make them using Python's os module)
* --nrays: The number of lightray objects to be generated
* --abun: Path to abundance file, if any. Defaults to solar abundances.
* --halo: Path to halo data.

# Notes for the user
* Make sure the abundance file passed does not have any unnecessary characters (i.e. characters otherr than the element's name and the value associated with it) as this will confuse Trident and result in a runtime error.
* SALSA requires a "center_list" that is the center of the galaxy around which the halo resides. As of right now, this list is a default value in the program and, intuitively, this center_list will need to be changed if a different halo is used. 
* If passing abundance patterns to the program, the user will see the print statement, "Buckle up bc we're about to get wRIGGITY WRIGGITY WRECKED SON". This is intentional; it's how the user will know that Trident is using the abundances that were passed. 