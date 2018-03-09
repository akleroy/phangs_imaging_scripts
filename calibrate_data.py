# This script runs the calibration scripts, first making sure that any
# additional flagging commands are added to the queue.

# Note that this is a pure python script. It *calls* casapy but it
# should not be run inside CASA.

# Edit the "Control Flow" section to use the script.

import os

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Control Flow
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ... a text list. The script will process only these galaxies.

only = []

# ... skip these galaxies

skip = []

# ... set as '12m' or '7m'. Leave it as None to process all data.

just_array = ['7m']

# Steps

overwrite_previous = False
update_flags = False
run_calibration = False
just_print_commands = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Define various routines and functions
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Run the script
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

array_list = ['12m', '7m']

# get the galaxy name and array and directory list

for gal in gals:
    
    if len(only) > 0:
        if only.count(gal) == 0:
            print "Skipping "+gal
            continue
    if len(skip) > 0:
        if skip.count(gal) > 0:
            print "Skipping "+gal
            continue

    for array in array_list:
        
        if len(just_array) > 0:
            if just_array.count(array) == 0:
                print "Skipping "+array
                continue

        # change directory to the relevant scripts/ directory

        # test if the calibrated/ directory is present.

        # if present and requested remove the calibrated data, else fast forward

        # edit the flagging file used for the pipeline run

        # run the calibration script with a non-interactive command call to CASA
            
