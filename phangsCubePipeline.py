# NB: CASA doesn't always include the pwd in the module search path. I
# had to modify my init.py file to get this to import. See the README.

import os
import numpy as np
import scipy.ndimage as ndimage
import glob

# Other PHANGS scripts
from phangsPipelinePython import *
import line_list

# Analysis utilities
import analysisUtils as au

# CASA imports
from taskinit import *

# Import specific CASA tasks
from concat import concat
from exportfits import exportfits
from flagdata import flagdata
from imhead import imhead
from immath import immath
from imstat import imstat
from imregrid import imregrid
from importfits import importfits
from impbcor import impbcor
from imsmooth import imsmooth
from makemask import makemask
from mstransform import mstransform
from split import split
from statwt import statwt
from tclean import tclean
from uvcontsub import uvcontsub
from visstat import visstat

# Physical constants
sol_kms = 2.9979246e5

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# DIRECTORY STRUCTURE AND FILE MANAGEMENT
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def rebuild_directories(outroot_dir=None):
    """
    Reset and rebuild the directory structure.
    """
    
    if outroot_dir is None:
        print("Specify a root directory.")
        return

    check_string = raw_input("Reseting directory structure in "+outroot_dir+". Type 'Yes' to confirm.")
    if check_string != "Yes":
        print("Aborting")
        return

    if os.is_dir(outroot_dir) == False:
        print("Got make "+outroot_dir+" manually. Then I will make the subdirectories.")

    for subdir in ['raw/','process/','products/','feather/','delivery/']:
        os.system('rm -rf '+outroot_dir+'/'+subdir)       
        os.system('mkdir '+outroot_dir+'/'+subdir)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# SET UP THE CUBES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def stage_cubes_in_casa(
    gal=None, array=None, product=None,
    outroot_dir=None, 
    overwrite=False,
    ):
    """
    Stage cubes for further processing in CASA. This imports them back
    into CASA as .image files.
    """
    if gal is None or array is None or product is None or \
            outroot_dir is None:
        print("Missing required input.")
        return
    
    input_dir = dir_for_gal(gal)
    in_cube_name = input_dir+gal+'_'+array+'_'+product+'.fits'
    in_pb_name = input_dir+gal+'_'+array+'_'+product+'_pb.fits'

    out_dir = outroot_dir+'raw/'
    out_cube_name = out_dir+gal+'_'+array+'_'+product+'.image'
    out_pb_name = out_dir+gal+'_'+array+'_'+product+'.pb'
    
    importfits(fitsimage=in_cube_name, imagename=out_cube_name,
               zeroblanks=False, overwrite=overwrite)
    importfits(fitsimage=in_pb_name, imagename=out_pb_name,
               zeroblanks=True, overwrite=overwrite)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# IMAGE MANIPULATION AND CLEANUP
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def primary_beam_correct(
    gal=None, array=None, product=None,
    root_dir=None, cutoff=0.5, overwrite=False):
    """
    Construct primary-beam corrected images.
    """

    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return

    input_dir = root_dir+'raw/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'.image'
    input_pb_name = input_dir+gal+'_'+array+'_'+product+'.pb'

    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_pbcorr.image'

    if overwrite:
        os.system('rm -rf '+output_cube_name)

    impbcor(imagename=input_cube_name,
            pbimage=input_pb_name,
            outfile=output_cube_name,
            cutoff=cutoff)

def convolve_to_round_beam(
    gal=None, array=None, product=None,
    root_dir=None,
    overwrite=False):
    """
    Convolve supplied image to have a round beam.
    """
    
    if gal is None or array is None or product is None or \
            root_dir is None:
        print("Missing required input.")
        return

    input_dir = root_dir+'process/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'_pbcorr.image'
    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_pbcorr_round.image'
    
    hdr = imhead(input_cube_name)
    if (hdr['axisunits'][0] != 'rad'):
        print("ERROR: Based on CASA experience. I expected units of radians.")
        print("I did not find this. Returning. Adjust code or investigate file "+input_cube_name)
        return
    pixel_as = abs(hdr['incr'][0]/np.pi*180.0*3600.)

    if (hdr['restoringbeam']['major']['unit'] != 'arcsec'):
        print("ERROR: Based on CASA experience. I expected units of arcseconds for the beam.")
        print("I did not find this. Returning. Adjust code or investigate file "+input_cube_name)
        return    
    bmaj = hdr['restoringbeam']['major']['value']
    
    target_bmaj = np.sqrt((bmaj)**2+(2.0*pixel_as)**2)

    imsmooth(imagename=input_cube_name,
             outfile=output_cube_name,
             targetres=True,
             major=str(target_bmaj)+'arcsec',
             minor=str(target_bmaj)+'arcsec',
             pa='0.0deg',
             overwrite=overwrite
             )

    input_dir = root_dir+'raw/'
    input_cube_name = input_dir+gal+'_'+array+'_'+product+'.image'

    output_dir = root_dir+'process/'
    output_cube_name = output_dir+gal+'_'+array+'_'+product+'_flat_round.image'    

    imsmooth(imagename=input_cube_name,
             outfile=output_cube_name,
             targetres=True,
             major=str(target_bmaj)+'arcsec',
             minor=str(target_bmaj)+'arcsec',
             pa='0.0deg',
             overwrite=overwrite
             )    

    return

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# ALIGN AND CONVERT SINGLE DISH DATA
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    
