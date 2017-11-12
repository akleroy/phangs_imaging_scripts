# This script stages the imaging. It copies the calibrated data into
# the imaging directories, runs any galaxy-specific processing, and
# then extracts line and continuum data sets ready for imaging.

import os
import phangsPipeline as pp
import analysisUtils as au

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Control Flow
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ... process only these galaxies
only = ['ngc0685', 'ngc1097_1', 'ngc1511', 
        'ngc1792_1','ngc1809','ngc2090',
        'ngc2283', 'ngc2775', 'ngc5134',
        'ic5237','ngc7456']

# ... skip these galaxies
skip = []

# ... set this to '12m' or '7m' to stage data only for those
# arrays. Leave it as None to process all data.
just_array = '7m'

# ... set these variables to indicate what steps of the script should
# be performed

do_copy = True
do_custom_scripts = True
do_extract_lines = True
do_extract_cont = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

gals = pp.list_gal_names()

for gal in gals:
    
    if len(only) > 0:
        if only.count(gal) == 0:
            continue
    if len(skip) > 0:
        if skip.count(gal) > 0:
            continue
        
    if do_copy:
        pp.copy_data(
            gal=gal,
            just_array=just_array,
            do_split=True,
            do_statwt=True,
            quiet=False)

    if do_custom_scripts:
        pass

    if do_extract_lines:
        pp.extract_phangs_lines(
            gal=gal,
            just_array=just_array,
            quiet=False)

    if do_extract_cont:
        pp.extract_phangs_continuum(
            gal=gal,
            just_array=just_array,
            quiet=False)


            
