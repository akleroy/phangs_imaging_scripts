"""
Stand alone routines to carry out basic noise estimation, masking, and
mask manipulation steps in CASA.
"""

# 
# 20200210 dzliu: moved "stat_clean_cube()" to here, as it is required by "signal_mask()"
# 20200210 dzliu: changed "casa." to "casaStuff.", as "casa" is a dict used by CASA itself.
# 20200210 dzliu: changed "print +(.*)$" to "logger.info(\1)"
# 

#region Imports and definitions

import os
import numpy as np
from scipy.special import erfc
import scipy.ndimage as ndimage
import pyfits # CASA has pyfits, not astropy
import glob

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Analysis utilities
import analysisUtils as au

# CASA stuff
import casaStuff

# Pipeline versionining
from pipelineVersion import version as pipeVer

#endregion

#region Noise estimation

def mad(
    data=None,
    as_sigma=True
    ):
    """
    Helper routine to calculate median absolute deviation (MAD). This
    is present already in scipy.stats but not in the version of scipy
    that CASA ships with. The MAD is a use a fast, useful robust noise
    estimator.

    data : the vector of data used to calculate the MAD. The routine
    flattens the array, so no along-axis operations.

    as_sigma (default True) : scale the output so that the returned
    value represents the RMS or 1-sigma value for a normal
    distribution. For Gaussian noise, this implies that the result can
    just be used as a standard noise estimate.  
    """

    if data is None:
        logger.error("No data supplied.")

    this_med = np.median(data)
    this_dev = np.abs(data - this_med)
    this_mad = np.median(this_dev)
    if as_sigma:
        return(this_mad / 0.6745)
    else:
        return(this_mad)

def estimate_noise(
    data=None,
    mask=None,
    method='mad',
    niter=None,
    ):
    """
    Return a noise estimate given a vector and associated mask.

    data : the data used to calculate the noise.

    mask : a mask used to indicate which subset of data to
    consider. In this routine mask values of True will be included in
    the calculation and mask values of False will be excluded. This
    matches the CASA syntax, but might require "inverting" a signal
    mask when the intention is to avoid bright signal.

    method (default "mad") : Method to use. Either "std" for standard
    deviation, "mad" for median absolute deviation, "chauvstd" for
    standard deviation with outlier rejection, or "chauvmad" for mad
    with outlier rejection.

    niter : number of iterations used in outlier rejection.

    Method "mad" is preferred for fast calculation and "chauvmad" for
    accurate calculation. Both should be reasonably robust.
    """
    
    if data is None:
        logger.error("No data supplied.")
        return(None)

    if mask is not None:
        if len(mask) != len(data):
            logger.error("Mask and data have mismatched sizes.")
            return(None)
    
    if niter is None:
        niter = 5

    valid_methods = ['std','mad','chauvstd','chauvmad']
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
        this_noise = mad(use_data, as_sigma=True)
        return(this_noise)

    if method == 'chauvstd' or method == 'chauvmad':
        for ii in range(niter):
            this_mean = np.mean(use_data)
            if method == 'chauvstd':
                this_std = np.std(use_data)
            elif method == 'chauvmad':
                this_std = mad(use_data, as_sigma=True)

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

def noise_for_cube(
    infile=None,
    maskfile=None,
    exclude_mask=True,
    method='mad',
    niter=None,
    ):
    """
    Get a single noise estimate for an image cube.
    """

    if infile is None:
        logger.error('No infile specified.')
        return(None)
    
    if not os.path.isdir(infile) and not os.path.isfile(infile):
        logger.error('infile specified but not found - '+infile)
        return(None)

    if maskfile is not None:
        if not os.path.isdir(maskfile) and not os.path.isfile(maskfile):
            logger.error('maskfile specified but not found - '+maskfile)
            return(None)
            
    myia = au.createCasaTool(casaStuff.iatool)
    myia.open(infile)
    data = myia.getchunk()
    mask = myia.getchunk(getmask=True)
    myia.close()

    if maskfile is not None:
        myia.open(maskfile)
        user_mask = myia.getchunk()
        user_mask_mask = myia.getchunk(getmask=True)
        myia.close()
        if exclude_mask:
            mask = mask*user_mask_mask*(user_mask < 0.5)
        else:
            mask = mask*user_mask_mask*(user_mask >= 0.5)

    this_noise = estimate_noise(
        data=data, mask=mask, method=method, niter=niter)
    
    return(this_noise)

def stat_cube(
    cube_file=None, 
    ):
    """
    Calculate statistics for an image cube. Right now this is a thin
    wrapper to imstat.
    """
    if cube_file == None:
        logger.info("No cube file specified. Returning")
        return
    
    imstat_dict = casaStuff.imstat(cube_file)
    
    return imstat_dict

#endregion

#region Mask creation and manipulation

def signal_mask(
    cube_root=None,
    out_file=None,
    suffix_in='',
    suffix_out='',
    operation='AND',
    high_snr = 4.0,
    low_snr = 2.0,
    absolute = False,
    ):
    """
    A simple signal mask creation routine used to make masks on the
    fly during imaging. Leverages CASA statistics and scipy.
    """
    
    if os.path.isdir(cube_root+'.image'+suffix_in) == False:
        logger.error('Data file not found: "'+cube_root+'.image'+suffix_in+'"')
        logger.info('Need CUBE_ROOT.image to be an image file.')
        logger.info('Returning. Generalize the code if you want different syntax.')
        return

    myia = au.createCasaTool(casaStuff.iatool)
    if operation == 'AND' or operation == 'OR':
        if os.path.isdir(cube_root+'.mask'+suffix_in) == True:
            myia.open(cube_root+'.mask'+suffix_in)
            old_mask = myia.getchunk()
            myia.close()
        else:
            logger.info("Operation AND/OR requested but no previous mask found.")
            logger.info("... will set operation=NEW.")
            operation = 'NEW'    

    if os.path.isdir(cube_root+'.residual'+suffix_in) == True:
        stats = stat_cube(cube_root+'.residual'+suffix_in)
    else:
        stats = stat_cube(cube_root+'.image'+suffix_in)
    rms = stats['medabsdevmed'][0]/0.6745
    hi_thresh = high_snr*rms
    low_thresh = low_snr*rms

    header = casaStuff.imhead(cube_root+'.image'+suffix_in)
    if header['axisnames'][2] == 'Frequency':
        spec_axis = 2
    else:
        spec_axis = 3

    myia.open(cube_root+'.image'+suffix_in)
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

    os.system('rm -rf '+cube_root+'.mask'+suffix_out)
    os.system('cp -r '+cube_root+'.image'+suffix_in+' '+cube_root+'.mask'+suffix_out)
    myia.open(cube_root+'.mask'+suffix_out)
    myia.putchunk(mask.astype(int)) #<TODO><DL># modified mask --> mask.astype(int)
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
        logger.info("Specify a cube root file name.")
        return

    myia = au.createCasaTool(casaStuff.iatool)    
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

def import_and_align_mask(  
    in_file=None,
    out_file=None,
    template=None,
    ):
    """
    Align a mask to a target astrometry. This includes some klugy
    steps (especially related to axes and interpolation) to make this
    work, e.g., for clean masks, most of the time.
    """

    # Import from FITS (could make optional)
    os.system('rm -rf '+out_file+'.temp_copy'+' 2>/dev/null')
    logger.debug('Importing mask file: "'+in_file+'"')
    casaStuff.importfits(fitsimage=in_file, 
               imagename=out_file+'.temp_copy', 
               overwrite=True)
    
    # Prepare analysis utility tool
    myia = au.createCasaTool(casaStuff.iatool)
    myim = au.createCasaTool(casaStuff.imtool)
    
    # Read mask data
    #myia.open(out_file+'.temp_copy')
    #mask = myia.getchunk(dropdeg=True)
    #myia.close()
    #print('**********************')
    #print('type(mask)', type(mask), 'mask.dtype', mask.dtype, 'mask.shape', mask.shape) # note that here the mask is in F dimension order, not Pythonic.
    #print(np.max(mask), np.min(mask))
    #print('**********************')
    
    # Read template image header
    hdr = casaStuff.imhead(template)
    #print('hdr', hdr)
    
    maskhdr = casaStuff.imhead(out_file+'.temp_copy')
    #print('maskhdr', maskhdr)
    
    # Check if 2D or 3D
    logger.debug('Template data axis names: '+str(hdr['axisnames'])+', shape: '+str(hdr['shape']))
    logger.debug('Mask data axis names: '+str(maskhdr['axisnames'])+', shape: '+str(maskhdr['shape']))
    is_template_2D = (np.prod(list(hdr['shape'])) == np.prod(list(hdr['shape'])[0:2]))
    is_mask_2D = (np.prod(list(maskhdr['shape'])) == np.prod(list(maskhdr['shape'])[0:2]))
    if is_template_2D and not is_mask_2D:
        logger.debug('Template image is 2D but mask is 3D, collapsing the mask over channel axes: '+str(np.arange(maskhdr['ndim']-1, 2-1, -1)))
        # read mask array
        myia.open(out_file+'.temp_copy')
        mask = myia.getchunk(dropdeg=True)
        myia.close()
        #print('**********************')
        #print('type(mask)', type(mask), 'mask.dtype', mask.dtype, 'mask.shape', mask.shape) # Note that here array shapes are in F dimension order, i.e., axis 0 is RA, axis 1 is Dec, axis 2 is Frequency, etc.
        #print('**********************')
        # 
        # collapse channel and higher axes
        #mask = np.any(mask.astype(int).astype(bool), axis=np.arange(maskhdr['ndim']-1, 2-1, -1)) # Note that here array shapes are in F dimension order, i.e., axis 0 is RA, axis 1 is Dec, axis 2 is Frequency, etc.
        #mask = mask.astype(int)
        #while len(mask.shape) < len(hdr['shape']):
        #    mask = np.expand_dims(mask, axis=len(mask.shape)) # Note that here array shapes are in F dimension order, i.e., axis 0 is RA, axis 1 is Dec, axis 2 is Frequency, etc.
        #print('**********************')
        #print('type(mask)', type(mask), 'mask.dtype', mask.dtype, 'mask.shape', mask.shape) # Note that here array shapes are in F dimension order, i.e., axis 0 is RA, axis 1 is Dec, axis 2 is Frequency, etc.
        #print('**********************')
        #os.system('rm -rf '+out_file+'.temp_collapsed'+' 2>/dev/null')
        ##myia.open(out_file+'.temp_copy')
        #newimage = myia.newimagefromarray(outfile=out_file+'.temp_collapsed', pixels=mask.astype(int), overwrite=True)
        #newimage.done()
        #myia.close()
        # 
        # collapse channel and higher axes
        os.system('rm -rf '+out_file+'.temp_collapsed'+' 2>/dev/null')
        myia.open(out_file+'.temp_copy')
        collapsed = myia.collapse(outfile=out_file+'.temp_collapsed', function='max', axes=np.arange(maskhdr['ndim']-1, 2-1, -1))
        collapsed.done()
        myia.close()
        # ia tools -- https://casa.nrao.edu/docs/CasaRef/image-Tool.html
        os.system('rm -rf '+out_file+'.temp_copy'+' 2>/dev/null')
        os.system('cp -r '+out_file+'.temp_collapsed'+' '+out_file+'.temp_copy'+' 2>/dev/null')
        os.system('rm -rf '+out_file+'.temp_collapsed'+' 2>/dev/null')
    
    # Align to the template grid
    os.system('rm -rf '+out_file+'.temp_aligned'+' 2>/dev/null')
    casaStuff.imregrid(imagename=out_file+'.temp_copy', 
                template=template, 
                output=out_file+'.temp_aligned', 
                asvelocity=True,
                interpolation='nearest',         
                replicate=False,
                overwrite=True)

    # Make an EXACT copy of the template, avoids various annoying edge cases
    os.system('rm -rf '+out_file+' 2>/dev/null')
    myim.mask(image=template, mask=out_file)
    
    # Pull the data out of the aligned mask and place it in the output file
    myia.open(out_file+'.temp_aligned')
    mask = myia.getchunk(dropdeg=True)
    myia.close()
    
    # 
    if is_template_2D and not is_mask_2D:
        while len(mask.shape) < len(hdr['shape']):
            mask = np.expand_dims(mask, axis=len(mask.shape))
        #print('**********************')
        #print('type(mask)', type(mask), 'mask.dtype', mask.dtype, 'mask.shape', mask.shape) # Note that here array shapes are in F dimension order, i.e., axis 0 is RA, axis 1 is Dec, axis 2 is Frequency, etc.
        #print('**********************')
        myia.open(out_file)
        data = myia.getchunk(dropdeg=False)
        data = mask
        myia.putchunk(data)
        myia.close()
    else:
        # Need to make sure this works for two dimensional cases, too.
        if (hdr['axisnames'][3] == 'Frequency') and (hdr['ndim'] == 4):
            myia.open(out_file)
            data = myia.getchunk(dropdeg=False)
            data[:,:,0,:] = mask
            myia.putchunk(data)
            myia.close()
        elif (hdr['axisnames'][2] == 'Frequency') and (hdr['ndim'] == 4):
            myia.open(out_file)
            data = myia.getchunk(dropdeg=False)
            data[:,:,:,0] = mask
            myia.putchunk(data)
            myia.close()
        else:
            logger.info("ALERT! Did not find a case.")

    os.system('rm -rf '+out_file+'.temp_copy'+' 2>/dev/null')
    os.system('rm -rf '+out_file+'.temp_aligned'+' 2>/dev/null')
    return()

#endregion
