"""
Routines to carry out basic noise estimation, masking, and mask
manipulation steps in CASA.
"""

#region Imports and definitions

import os
import numpy as np
from scipy.special import erfc
import pyfits # CASA has pyfits, not astropy
import glob

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Analysis utilities
import analysisUtils as au

# CASA stuff
import casaStuff as casa

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

#region Noise estimation

def estimate_noise(
    data=None,
    mask=None,
    method='mad',
    niter=3,
    ):
    """
    Return a noise estimate given a vector and associated mask.
    """
    
    if data is None:
        logger.error("No data supplied.")
        return(None)

    if mask is not None:
        if len(mask) != len(data):
            logger.error("Mask and data have mismatched sizes.")
            return(None)
    
    valid_methods = ['std','mad','chauv']
    if method not in valid_methods:
        logger.error("Invalid method - "+method+" valid methods are "+str(valid_methods))
        return(None)
    
    if mask is None:
        use_mask = np.isfinite(data)
    else:
        use_mask = mask*np.isfinite(data)

    if np.sum(use_mask) == 0:
        logger.error("No valid data. Returning NaN.")
        return(np.nan)

    use_data = data[use_mask]

    if method == 'std':
        this_noise = np.std(use_data)
        return(this_noise)

    if method == 'mad':
        this_med = np.median(use_data)
        this_dev = np.abs(use_data - this_med)
        this_mad = np.median(this_dev)
        this_noise = this_mad / 0.6745
        return(this_noise)

    if method == 'chauv':
        for ii in range(niter):
            this_mean = np.mean(use_data)
            this_std = np.std(use_data)
            this_dev = np.abs((use_data-this_mean)/this_std)/2.0**0.5
            this_prob = erfc(this_dev)
                        
            chauv_crit = 1.0/(2.0*len(use_data))
            keep = this_prob > chauv_crit
            if np.sum(keep) == 0:
                logger.error("Rejected all data. Returning NaN.")
                return(np.nan)
            use_data = use_data[keep]
            
        this_noise = np.std(use_data)
        return(this_noise)

    return(None)
    
def test_noise(
    ):
    """
    Test the noise estimation routine.
    """
    
    tol = 1e-2

    vec = np.random.randn(1e5)
    mad_est = estimate_noise(vec, method='mad')
    std_est = estimate_noise(vec, method='std')
    chauv_est = estimate_noise(vec, method='chauv')
    
    logger.info("mad estimate accuracy: "+str(np.abs(mad_est-1.0)))
    if np.abs(mad_est - 1.0) > tol:
        logger.error("mad estimate exceeds tolerance.")

    logger.info("std estimate accuracy: "+str(np.abs(std_est-1.0)))
    if np.abs(std_est - 1.0) > tol:
        logger.error("std estimate exceeds tolerance.")

    logger.info("chauv estimate accuracy: "+str(np.abs(chauv_est-1.0)))
    if np.abs(chauv_est - 1.0) > tol:
        logger.error("chauv estimate exceeds tolerance.")

def noise_for_cube(
    infile=None,
    mask=None,
    method='mad',
    ):
    """
    Get a single noise estimate for an image cube.
    """

    if infile is None:
        logger.error('No infile specified.')
        return(None)
    
    if not os.path.isdir(infile) and not os.path.isfile(infile):
        logger.error('Infile specified but not found - '+infile)
        return(None)
        
    imstat_dict = imstat(cube_file)
    
    return imstat_dict

#endregion

#region Mask creation and manipulation

def signal_mask(
    cube_root=None,
    out_file=None,
    operation='AND',
    high_snr = 4.0,
    low_snr = 2.0,
    absolute = False,
    ):
    """
    A simple signal mask creation routine used to make masks on the
    fly during imaging. Leverages CASA statistics and scipy.
    """
    
    if os.path.isdir(cube_root+'.image') == False:
        print 'Need CUBE_ROOT.image to be an image file.'
        print 'Returning. Generalize the code if you want different syntax.'
        return

    myia = au.createCasaTool(iatool)
    if operation == 'AND' or operation == 'OR':
        if os.path.isdir(cube_root+'.mask') == True:
            myia.open(cube_root+'.mask')
            old_mask = myia.getchunk()
            myia.close()
        else:
            print "Operation AND/OR requested but no previous mask found."
            print "... will set operation=NEW."
            operation = 'NEW'    

    if os.path.isdir(cube_root+'.residual') == True:
        stats = stat_clean_cube(cube_root+'.residual')
    else:
        stats = stat_clean_cube(cube_root+'.image')
    rms = stats['medabsdevmed'][0]/0.6745
    hi_thresh = high_snr*rms
    low_thresh = low_snr*rms

    header = imhead(cube_root+'.image')
    if header['axisnames'][2] == 'Frequency':
        spec_axis = 2
    else:
        spec_axis = 3

    myia.open(cube_root+'.image')
    cube = myia.getchunk()
    myia.close()

    if absolute:
        hi_mask = (np.abs(cube) > hi_thresh)
    else:
        hi_mask = (cube > hi_thresh)
    mask = \
        (hi_mask + np.roll(hi_mask,1,axis=spec_axis) + \
             np.roll(hi_mask,-1,axis=spec_axis)) >= 1

    if high_snr > low_snr:
        if absolute:
            low_mask = (np.abs(cube) > low_thresh)
        else:
            low_mask = (cube > low_thresh)
        rolled_low_mask = \
            (low_mask + np.roll(low_mask,1,axis=spec_axis) + \
                 np.roll(low_mask,-1,axis=spec_axis)) >= 1
        mask = ndimage.binary_dilation(hi_mask, 
                                       mask=rolled_low_mask, 
                                       iterations=-1)

    if operation == 'AND':
        mask = mask*old_mask
    if operation == 'OR':
        mask = (mask + old_mask) > 0
    if operation == 'NEW':
        mask = mask

    os.system('rm -rf '+cube_root+'.mask')
    os.system('cp -r '+cube_root+'.image '+cube_root+'.mask')
    myia.open(cube_root+'.mask')
    myia.putchunk(mask)
    myia.close()

def apply_additional_mask(
    old_mask_file=None,
    new_mask_file=None,
    new_thresh=0.0,
    operation='AND'
    ):
    """
    Combine a mask with another mask on the same grid and some
    threshold. Can run AND/OR operations. Can be used to apply primary
    beam based masks by setting the PB file to new_mask_file and the
    pb_limit as new_thresh.
    """
    if root_mask == None:
        print "Specify a cube root file name."
        return

    myia = au.createCasaTool(iatool)    
    myia.open(new_mask_file)
    new_mask = myia.getchunk()
    myia.close()

    myia.open(old_mask_file)
    mask = myia.getchunk()
    if operation == "AND":
        mask *= (new_mask > new_thresh)
    else:
        mask = (mask + (new_mask > new_thresh)) >= 1.0
    myia.putchunk(mask)
    myia.close()

    return

#endregion

#region Mask alignment

def import_and_align_mask(  
    in_file=None,
    out_file=None,
    template=None,
    ):
    """
    Align a mask to a target astrometry. Some klugy steps to make
    things work most of the time.
    """

    # Import from FITS (could make optional)
    os.system('rm -rf '+out_file+'.temp_copy')
    importfits(fitsimage=in_file, 
               imagename=out_file+'.temp_copy'
               , overwrite=True)

    # Align to the template grid
    os.system('rm -rf '+out_file+'.temp_aligned')
    imregrid(imagename=out_file+'.temp_copy', 
             template=template, 
             output=out_file+'.temp_aligned', 
             asvelocity=True,
             interpolation='nearest',         
             replicate=False,
             overwrite=True)

    # Make an EXACT copy of the template, avoids various annoying edge cases
    os.system('rm -rf '+out_file)
    myim = au.createCasaTool(imtool)
    myim.mask(image=template, mask=out_file)

    hdr = imhead(template)

    # Pull the data out of the aligned mask and place it in the output file
    myia = au.createCasaTool(iatool)
    myia.open(out_file+'.temp_aligned')
    mask = myia.getchunk(dropdeg=True)
    myia.close()

    # Need to make sure this works for two dimensional cases, too.
    if (hdr['axisnames'][3] == 'Frequency') and \
            (hdr['ndim'] == 4):
        myia.open(out_file)
        data = myia.getchunk(dropdeg=False)
        data[:,:,0,:] = mask
        myia.putchunk(data)
        myia.close()
    elif (hdr['axisnames'][2] == 'Frequency') and \
            (hdr['ndim'] == 4):
        myia.open(mask_root+'.mask')
        data = myia.getchunk(dropdeg=False)
        data[:,:,:,0] = mask
        myia.putchunk(data)
        myia.close()
    else:
        print "ALERT! Did not find a case."

    os.system('rm -rf '+out_file+'.temp_copy')
    os.system('rm -rf '+out_file+'.temp_aligned')
    return

#endregion
