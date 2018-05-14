# This script stages the imaging. It copies the calibrated data into
# the imaging directories, regrids to the desired velocity axis, and
# then extracts line and continuum data sets ready for imaging.

# Edit the "Control Flow" section to use the script.

import os
import phangsPipeline as pp
import analysisUtils as au
import glob

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Control Flow
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ... a text list. The script will process only these galaxies.

only = []

# ... skip these galaxies.

skip = []

# ... set this to '12m' or '7m' to stage data only for those
# arrays. Leave it as None to process all data. If both 12m and 7m
# data are processed, then the script will also create 12m+7m data. So
# you need to rerun the staging when both data sets arrive.

just_array = '7m'

# ... set these variables to indicate what steps of the script should
# carry out. The steps do:

# do_copy - copy the data from the calibrated data directory to the
# working directory. Will create the working directory for this galaxy
# if it doesn't exist already. Uses "ms_file_key.txt" and maps
# multi-part galaxies to directories using "dir_key.txt"

# do_custom_scripts - NOT IMPLEMENTED YET. Will run custom processing
# for each galaxy. For example, additional flagging and uv continuum
# subtraction.

# do_extract_lines - extract lines from the measurement set and regrid
# them onto our working grid. The line set is assumed to be 12co21 and
# c18o21 for PHANGS. The velocity grid is defined in
# "mosaic_definitions.txt" and this step should be rerun after
# changing that grid. Also makes "channel 0" (line integrated) data
# sets.

# do_extract_cont - extract a single-channel continuum data set from
# the measurement set, first flagging lines. The velocity windows used
# for flagging lines is set in "mosaic_definitions.txt"

do_copy = True
do_custom_scripts = True
do_extract_lines = True
do_extract_cont = True
do_only_new = False

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
    
    if do_only_new:        
        this_dir = pp.dir_for_gal(gal)
        has_12m = len(glob.glob(this_dir+gal+'*12m_co21.ms')) > 0
        has_7m = len(glob.glob(this_dir+gal+'*7m_co21.ms')) > 0
        if has_12m or has_7m:
            print ""
            print "... You requested to only stage new data."
            print "... I found an existing file for "+gal+" ."
            print "... I will skip that galaxy."
            print ""
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


            
