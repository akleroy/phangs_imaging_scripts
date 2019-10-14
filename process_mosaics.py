# This script prepares and then executes the linear mosaicking

import os
import phangsPipeline as pp
import phangsCubePipeline as pcp
import analysisUtils as au
import glob

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Directories and definitions
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

interferometric_array_list = ['12m', '7m', '12m+7m']
full_array_list = ['12m+7m+tp', '12m+7m', '12m', '7m', '7m+tp']
full_product_list = ['co21','c18o21','13co21']
mosaic_key = pp.mosaic_key()

inroot_dir = '../'
vstring = 'v3_casa'
outroot_dir = '../release/'+vstring+'/'

cutoff = 0.25

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

# ... set as '12m', '7m', '7m+tp', '12m+7m', or '12m+7m+tp' to process
# only data from that array. Leave it as None to process all data.

#just_array = []
just_array = ['7m','7m+tp']

# ... set as the products to be handled. Valid choices for the basic
# PHANGS data are 'co21', 'c18o21', 'cont', 'co21_chan0', and
# 'c18o21_chan0'. Note that right now cont and chan0 are not tested.

just_product = ['co21']

# ... set these variables to indicate what steps of the script should
# be performed.

do_common_res = False
do_align = False
do_mosaic = True
do_cleanup_cubes = False

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop over all galaxies to stage, process, mosaic, and cleanup
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

before_first = True
after_last = False

for this_gal in mosaic_key.keys():

    if len(only) > 0:
        not_in_list = (only.count(this_gal) == 0)
        if not_in_list:
            print("Skipping galaxyt "+this_gal)
            continue

    if len(skip) > 0:
        in_list = (skip.count(this_gal) > 0)
        if in_list:
            print("Skipping galaxy part "+this_gal)
            continue

    if first != "":
        if this_gal == first:
            before_first = False
            if before_first:
                continue
                
    if last != "":
        if after_last == True:
            continue
        if this_gal == last:
            after_last = True

    print("Mosaicking for galaxy: "+this_gal)
    print("... which has parts: ")
    print(mosaic_key[this_gal])

    for product in full_product_list:
            
        if len(just_product) > 0:
            if just_product.count(product) == 0:
                print("... skipping line/product: "+product)
                continue
                        
        print("... processing line/product: "+product)
            
        for array in full_array_list:

            if len(just_array) > 0:
                if just_array.count(array) == 0:
                    print("... ... skipping array "+array)
                    continue

            print("... ... processing array: "+array)

            if do_common_res:
                pcp.phangs_common_res_for_mosaic(
                    gal=this_gal, array=array, product=product,
                    root_dir=outroot_dir,
                    overwrite=True)

            if do_align:
                pcp.phangs_align_for_mosaic(
                    gal=this_gal, array=array, product=product,
                    root_dir=outroot_dir,
                    overwrite=True)

            if do_mosaic:
                pcp.phangs_mosaic_data(
                    gal=this_gal, array=array, product=product,
                    root_dir=outroot_dir,
                    overwrite=True)

            if do_cleanup_cubes:
                pcp.phangs_cleanup_cubes(
                    gal=this_gal, array=array, product=product,
                    root_dir=outroot_dir,
                    overwrite=True)
