# This script post-processes the imaging into science-ready cubes.

import os
import phangsPipeline as pp
import phangsCubePipeline as pcp
import analysisUtils as au
import glob

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Directories and definitions
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

interferometric_array_list = ['12m', '7m', '12m+7m']
full_product_list = ['co21','c18o21','13co21']
gal_part_list = pp.list_gal_names()

inroot_dir = '../'
vstring = 'v3_casa'
outroot_dir = '../release/'+vstring+'/'

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Control Flow
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ... a text list. The script will process only these galaxies.

only = []

# ... skip these galaxies
skip = []

# ... start with this galaxy

first = ""
last = ""

# ... set as '12m', '7m', or '12m+7m' to process only data from that
# array. Leave it as None to process all data.

just_array = ['7m']

# ... set as the products to be handled. Valid choices for the basic
# PHANGS data are 'co21', 'c18o21', 'cont', 'co21_chan0', and
# 'c18o21_chan0'. Note that right now cont and chan0 are not tested.

just_product = ['co21']

# ... set these variables to indicate what steps of the script should
# be performed.

rebuild_directories = False

stage_cubes = True
primary_beam_correct = True
convolve_to_round_beam = True

stage_single_dish = True
write_feather_script = True
copy_feathered_data = True

stage_mosaicking = True
mosaic_data = True
export_and_cleanup = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop #0 - wipe and rebuild if requested
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

if rebuild_directories:
    pcp.rebuild_directories(outroot_dir=outroot_dir)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop #1 - clean up our existing cubes
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

before_first = True
after_last = False

for gal in gal_part_list:
    
    if len(only) > 0:
        if only.count(gal) == 0:
            print "Skipping "+gal
            continue

    if len(skip) > 0:
        if skip.count(gal) > 0:
            print "Skipping "+gal
            continue

    if first != "":
        if gal == first:
            before_first = False
        if before_first:
            continue
    
    if last != "":
        if after_last == True:
            continue
        if gal == last:
            after_last = True
 
    for array in interferometric_array_list:

        if len(just_array) > 0:
            if just_array.count(array) == 0:
                print "Skipping "+array
                continue

        for product in full_product_list:

            if len(just_product) > 0:
                if just_product.count(product) == 0:
                    print "Skipping "+product
                    continue

            print gal, array, product

            if stage_cubes:
                pcp.stage_cubes_in_casa(
                    gal=gal, array=array, product=product, 
                    outroot_dir = outroot_dir, 
                    overwrite=True
                    )

            if primary_beam_correct:
                pcp.primary_beam_correct(
                    gal=gal, array=array, product=product, 
                    root_dir=outroot_dir,
                    overwrite=True
                    )

            if convolve_to_round_beam:
                pcp.convolve_to_round_beam(
                    gal=gal, array=array, product=product, 
                    root_dir=outroot_dir,
                    overwrite=True
                    )                
            
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop #2 - mosaic multi-part cubes
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


