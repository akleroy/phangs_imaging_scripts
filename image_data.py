# This script carries out the imaging.

import os
import phangsPipeline as pp
import analysisUtils as au

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Control Flow
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ... process only these galaxies
only = []

# ... skip these galaxies
skip = []

# ... set as '12m', '7m', or '12m+7m' to process only data from that
# array. Leave it as None to process all data.
just_array = ['7m']

# ... set as the products to be handled. Valid choices for the basic
# PHANGS data are 'co21', 'c18o21', 'cont', 'co21_chan0',
# 'c18o21_chan0'
just_product = ['co21']

# ... set these variables to indicate what steps of the script should
# be performed
 
make_dirty_image=True
revert_to_dirty=False
read_in_clean_mask=True
run_multiscale_clean=True
revert_to_multiscale=False
make_singlescale_mask=True
run_singlescale_clean=True
export_to_fits=True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

gals = pp.list_gal_names()

array_list = ['12m', '7m', '12m+7m']

product_list = ['co21','c18o21','cont','co21_chan0','c18o21_chan0']

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

        for product in product_list:

            if len(just_product) > 0:
                if just_product.count(product) == 0:
                    print "Skipping "+product
                    continue

            print gal, array, product

            clean_call = \
                pp.buildPhangsCleanCall(
                gal=gal,
                array=array,
                product=product,
                tag='')

            pp.phangsImagingRecipe(
                clean_call=clean_call,
                make_dirty_image=make_dirty_image,
                revert_to_dirty=revert_to_dirty,
                read_in_clean_mask=read_in_clean_mask,
                run_multiscale_clean=run_multiscale_clean,
                revert_to_multiscale=revert_to_multiscale,
                make_singlescale_mask=make_singlescale_mask,
                run_singlescale_clean=run_singlescale_clean,
                do_export_to_fits=export_to_fits,
                )
            
