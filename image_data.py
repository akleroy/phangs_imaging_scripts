# This script carries out the imaging. It makes a dirty image, reads
# and aligns the user-supplied clean mask, carries out multiscale
# imaging down to S/N ~ 4, builds a mask of bright emission, carries
# out a single scale clean down to low S/N within this mask, and then
# exports the output to FITS.

# Right now it should work well for CO21 and C18O21. Imaging of the
# two-d products is still TBD.

# Edit the "Control Flow" section to use the script.

# WARNING! Right now a bug prevents the script from stopping and
# starting effectively, so it needs to be run end-to-end to work. This
# shouldn't be a big problem for the 7m data.

import os
import phangsPipeline as pp
import analysisUtils as au
import glob

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Control Flow
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ... a text list. The script will process only these galaxies.

only = []

if False:
    only = ['ic1954',
            'ic5273',
            'ic5332',
            'ngc0628',
            'ngc0685',
            'ngc1087',
            'ngc1097_1',
            'ngc1097_2',
            'ngc1300',
            'ngc1365',
            'ngc1385',
            'ngc1433',
            'ngc1511',
            'ngc1512',
            'ngc1546',
            'ngc1559',
            'ngc1566',
            'ngc1637',
            'ngc1672',
            'ngc1792_1',
            'ngc1792_2',
            'ngc1809',
            'ngc2090',
            'ngc2283',
            'ngc2566',
            'ngc2775',
            'ngc2835',
            'ngc2903_1',
            'ngc2903_2',
            'ngc2903_3',
            'ngc2997_1',
            'ngc2997_2',
            'ngc2997_3',
            'ngc3059',
            'ngc3137',
            'ngc3239',
            'ngc3351',
            'ngc3507',
            'ngc3511',
            'ngc3521_1', 
            'ngc3521_2', 
            'ngc3596', 
            'ngc3621_1', 
            'ngc3621_2',
            'ngc3626']

if True:
    only = ['ngc3627north',
            'ngc3627south',
            'ngc4207', 
            'ngc4254north', 
            'ngc4254south', 
            'ngc4293', 
            'ngc4298', 
            'ngc4303', 
            'ngc4321north',
            'ngc4321south', 
            'ngc4424', 
            'ngc4457', 
            'ngc4496a', 
            'ngc4535', 
            'ngc4536_1', 
            'ngc4536_2',
            'ngc4540', 
            'ngc4548', 
            'ngc4569', 
            'ngc4571', 
            'ngc4579', 
            'ngc4654', 
            'ngc4689', 
            'ngc4694',
            'ngc4731', 
            'ngc4781', 
            'ngc4826', 
            'ngc4941',  
            'ngc4951', 
            'ngc5042', 
            'ngc5068north', 
            'ngc5068south',
            'ngc5128', 
            'ngc5134', 
            'ngc5248_1', 
            'ngc5248_2', 
            'ngc5530',
            'ngc5643_1',
            'ngc5643_2', 
            'ngc6300', 
            'ngc6744north', 
            'ngc6744south', 
            'ngc7456', 
            'ngc7496']
    
# ... skip these galaxies

skip = []

# ... set as '12m', '7m', or '12m+7m' to process only data from that
# array. Leave it as None to process all data.

just_array = ['7m']

# ... set as the products to be handled. Valid choices for the basic
# PHANGS data are 'co21', 'c18o21', 'cont', 'co21_chan0', and
# 'c18o21_chan0'. Note that right now cont and chan0 are not tested.

just_product = ['co21']

# ... set these variables to indicate what steps of the script should
# be performed. The steps do:

# make_dirty_image - make a niter=0 image cube. Useful for checking
# mosaic parameters and that sort of thing. Used as a template for the
# clean mask, so this needs to be done first.

# revert_to_dirty - (DOES NOT WORK RIGHT NOW) reset the whole process
# so that the cube is now the dirty cube. Unfortunately there's some
# 'memory' related to the clean call in CASA that is causing a bug so
# that things don't cleanly resume from here.

# read_in_clean_mask - if a clean mask is found in the ../clean_masks/
# directory, read it in and align it to the astrometric grid of the
# dirty image. This is now the .mask file and will be used in future
# imaging.

# run_multiscale_clean - run a multiscale clean using a set of scales
# selected for that array combination. By default, clean down to a
# signal to noise of 4 in the residuals or stop when successive
# iterations of clean change the flux in the model by less than 2%.

# revert_to_multiscale - (DOES NOT WORK RIGHT NOW) reset the process
# to just after the multiscale clean. See revert_to_dirty above.

# make_singlescale_mask - construct a signal-to-noise based mask for
# use in single scale clean. Uses the .image to do this, so that
# bright regions in the image so far get more cleaning. This maks is
# joined with the clean mask, so no regions are included outside the
# clean mask.

# run_singlescale_clean - run a single scale hogbom clean inside the
# bright regions defined above. Clean with a very low (S/N ~ 1)
# threshold and stop based on convergence in flux between successive
# model images.

# export_to_fits - export the dirt, multiscale, and final images and
# associated products to FITS files for subsequent processing.
 
# Right now I recommend to not change the flow, and to rerun imaging
# from beginning to end when using the script.

make_dirty_image=True
revert_to_dirty=False
read_in_clean_mask=True
run_multiscale_clean=True
revert_to_multiscale=False
make_singlescale_mask=True # True
run_singlescale_clean=True # True
export_to_fits=True

do_only_new = False # True

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

            if do_only_new:
                this_dir = pp.dir_for_gal(gal)
                out_image_name = this_dir+gal+'_'+array+'_'+product+'.image'
                has_image = len(glob.glob(out_image_name)) > 0
                if has_image:
                    print ""
                    print "... You requested to only image new data."
                    print "... I found an existing image named "+out_image_name+" ."
                    print "... I will skip this combination of galaxy, array, and product."
                    print ""
                    continue

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
            
