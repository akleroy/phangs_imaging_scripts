import logging
import scipy.ndimage.morphology as morph
import scipy.ndimage as nd
from scipy.signal import savgol_coeffs
import numpy as np
from astropy.stats import mad_std
from astropy.convolution import convolve, Gaussian2DKernel
import scipy.stats as ss

from spectral_cube import SpectralCube
import astropy.wcs as wcs
import astropy.units as u
# from pipelineVersion import version as pipeVer
from astropy.io import fits

from scNoiseRoutines import mad_zero_centered

np.seterr(divide='ignore', invalid='ignore')

mad_to_std_fac = 1.482602218505602

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def nchan_thresh_mask(cube, thresh=5., nchan=2):
    """
    Make mask from a cube applying a threshold across some number of
    consecutive channels.

    Parameters:

    -----------

    cube : np.array

        A cube in the same units as threshold.

    Keywords:
    ---------

    thresh : float

        Threshold for inclusion in mask if.

     nchan : int

        Number of consecutive channels at thresh required for a detection. Default: 2

    """

    # TBD Error checking on types, dimensionality, etc.

    mask = np.greater_equal(cube, thresh,
                            where=(~np.isnan(cube)),
                            out=np.full(cube.shape, False, dtype=bool))
    
    kernel = np.ones(nchan, dtype=np.bool)
    kernel = kernel[:,np.newaxis,np.newaxis]

    mask = morph.binary_opening(mask, kernel)

    return(mask)
    
def reject_small_regions(mask, min_volume=0, min_area=0):
    """
    Remove small regions from a mask. Small can be defined in either
    volume or area.

    Parameters:

    -----------

    mask : np.array
        A mask.

    Keywords:
    ---------

    minvolume : int    
        Minimum volume in pixels. Default 0.

    minarea : int
        Minimum area in pixels. Default 0.

    """

    # TBD Error checking on types, dimensionality, etc.

    # Blob color the cube and loop over regions
    
    regions, regct = nd.label(mask)
    objslices = nd.find_objects(regions)

    for ii, thisslice in enumerate(objslices):
        subcube = regions[thisslice]

        if min_volume > 0:
            volume =  np.sum(subcube == (ii+1))
            if pixct < min_volume:
                mask[regions == (ii+1)] = False
            
        if min_area > 0:
            if mask.ndim == 3:
                area = np.sum(np.any(subcube == (ii+1), axis=0))
            if mask.ndim == 2:
                area = np.sum(np.any(subcube == (ii+1)))

            if area < min_area:
                mask[regions == (ii+1)] = False

    return(mask)

def grow_mask(mask, iters_xy=0, iters_v=0, constraint=None):
    """
    Grow a boolean mask via dilation in the spectral (v) dimension,
    spatial (xy) dimension, or into a constraint.

    Logic:

    (Case I) If iters_xy and/or iters_v suppplied:

    1. Dilate original mask by iters_xy

    2. Dilate original mask by iters_v

    3. Take the union of the two new masks.

    4. Apply the constraint.

    (Case II) Only a constraint is supplied:

    1. Include all regions of the constraint that include an element
    of the original mask.

    Parameters:

    -----------

    mask : np.array
        A mask.

    Keywords:
    ---------

    iters_xy : int    
        Number of iterations of expansion in spatial dimensions.

    iters_v : int    
        Number of iterations of expansion in spectral dimension.

    constraint : np.array that can be broadcast to mask
        Another mask to use as a constraint.

    """

    # TBD Error checking on types, dimensionality, etc.

    if iters_xy > 0:

        # Generate the structure to dilate by
        struct = morph.iterate_structure(
            morph.generate_binary_structure(2, 1), iters_xy)

        # Fill in a third dimension if needed
        if mask.ndim == 3:
            struct = struct[np.newaxis, :, :]

        if iters_v > 0:
            mask_xy = morph.binary_dilation(mask, struct)
        else:
            mask = morph.binary_dilation(mask, struct)

    if iters_v > 0:

        struct = np.ones(grow_v, dtype=np.bool)
        struct = struct[:, np.newaxis, np.newaxis]

        if iters_xy > 0:
            mask_v = morph.binary_dilation(mask, struct)
        else:
            mask = morph.binary_dilation(mask, struct)

    if iters_v > 0 and iters_xy > 0:
        mask = np.logical_or(mask_v, mask_xy)

    if (iters_v > 0 or iters_xy > 0) and (constraint is not None):
        mask = np.logical_and(mask, constraint)

    if (iters_v == 0 and iters_xy == 0) and (constraint is not None):

        # blob color regions in the constraint
        regions, regct = nd.label(constraint)

        # get a list of all region assignments that have a True value
        # in the original mask
        good_regions = np.unique(regions[mask])

        # create a new mask that includes only these good new regions
        mask = np.zeros_like(mask, dtype=np.bool)
        for hit in good_regions:
            mask[regions == hit] = True

    return(mask)

def cprops_mask(data, noise=None, 
                hi_thresh=5, hi_nchan=2,
                lo_thresh=None, lo_nchan=None,
                min_pix=None, min_area=None,
                grow_xy=None, grow_v=None, 
                invert=False):
    """
    Standard CPROPS masking recipe.

    Parameters:

    -----------

    data : np.array

        Original data

    noise : np.array

        Estimate of the amplitude of the noise at every position in the data
        (or an array that will broadcast to that under division).

    Keywords:
    ---------

    hi_thresh : float
        Threshold for detection (in units of sigma).  Default: 5

    hi_nchan : int
        Number of consecutive channels needed for detection.  Default: 2

    lo_thresh : float
        Threshold for inclusion in mask if connected to a hi_thresh
        detection. Default: None

    lo_nchan : int
        Number of consecutive channels at lo_thresh required for a detection
        if connected to a hi_thresh detection. Default: None

    min_pix : int
        Number of pixels required for a detection.  Default: None

    min_area : int
        Minimum number of pixels required in area projection for a detection.
        Default: None

    grow_xy : int
        Not implemented

    grow_v : int
        Not implemented

    invert : bool
        If True, invert the data before applying masking.

        Used for assessing the number of false positives given masking
        criteria. Default: False.

    """

    # TBD error checking, dimensions, types, etc.

    if noise is None:
        logger.warning("Need a noise estimate. Making a simple one.")
        noise = mad_zero_centered(data)

    # Recase the cube into a signal-to-noise cube
    signif = data / noise

    # If requested, invert the data
    if invert:
        signif *= -1

    # Create a the core mask

    hi_mask = nchan_thresh_mask(
        signif, thresh=hi_thresh, nchan=hi_nchan)

    # If requested, reject small regions from the mask
    if (min_pix is not None) or (min_area is not None):
        
        if min_pix is None:
            min_pix = 0

        if min_area is None:
            min_area = 0
        
        hi_mask = reject_small_regions(
            hi_mask, min_volume=min_pix, min_area=min_area)

    # If supplied, make a lower significance mask and expand into it

    if (lo_thresh is not None) and (lo_nchan is not None):
        lo_mask = nchan_thresh_mask(
            signif, thresh=lo_thresh, nchan=lo_nchan)

        # Now expand the original mask into the lower mask
        mask = grow_mask(hi_mask, constraint=lo_mask)
        
    # If requested, grow the mask in xy and v directions. Sequential
    # calls mean that the xy is applied then the v.

    if grow_xy is not None:
        mask = grow_mask(mask, xy_iters=grow_xy)
    
    if grow_v is not None:
        mask = grow_mask(mask, v_iters=grow_v)
    
    return(mask)

def recipe_hybridize_mask(hires_in, lores_in, order='bilinear',
                          return_cube=False):

    if type(hires_in) is str:
        hires_hdulist = fits.open(hires_in)
        hires = SpectralCube(np.array(hires_hdulist[0].data,
                                      dtype=np.bool),
                             wcs=wcs.WCS(hires_hdulist[0].header),
                             header=hires_hdulist[0].header)
    elif type(hires_in) is SpectralCube:
        hires = hires_in
    else:
        logging.error('Unrecognized input type')
        raise NotImplementedError

    if type(lores_in) is str:
        lores_hdulist = fits.open(lores_in)
        lores = SpectralCube(np.array(lores_hdulist[0].data,
                                      dtype=np.bool),
                             wcs=wcs.WCS(lores_hdulist[0].header),
                             header=lores_hdulist[0].header)
    elif type(lores_in) is SpectralCube:
        lores = lores_in
    else:
        logging.error('Unrecognized input type')
        raise NotImplementedError

    lores = lores.reproject(hires.header, order=order)
    lores = lores.spectral_interpolate(hires.spectral_axis)
    mask = np.logical_or(np.array(hires.filled_data[:].value, dtype=np.bool),
                         np.array(lores.filled_data[:].value, dtype=np.bool))
    if return_cube:
        mask = SpectralCube(mask, wcs=wcs.WCS(hires_hdulist[0].header),
                            header=hires_hdulist[0].header,
                            meta={'BUNIT': ' ', 'BTYPE': 'Mask'})
    return(mask)

def recipe_phangs_strict_mask(
    incube, innoise, outfile=None, 
    mask_kwargs=None, 
    return_spectral_cube=False,
    overwrite=False):
    """
    Task to create the PHANGS-style "strict" masks.

    Parameters:

    -----------

    cube : string or SpectralCube

        The cube to be masked.

    noise : string or SpectralCube

        The noise estimate.

    Keywords:
    ---------

    outfile : string

        Filename where the mask will be written. The mask is also
        returned.

    masks_kwargs : dictionary

        Parameters to be passed to the masking routine.

    """

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Error checking and work out inputs
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # TBD error checking, dimensions, types, etc.

    if type(incube) is SpectralCube:
        cube = incube
    elif type(incube) == type("hello"):
        cube = SpectralCube.read(incube)
    else:
        logger.error("Input cube must be a SpectralCube object or a filename.")

    if type(innoise) is SpectralCube:
        rms = innoise
    elif type(innoise) == type("hello"):
        rms = SpectralCube.read(innoise)
    else:
        logger.error("Input noise must be a SpectralCube object or a filename.")

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Set up the masking kwargs
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # Initialize an empty kwargs dictionary
    if mask_kwargs is None:
        mask_kwargs = {}

    # Fill in strict mask defaults
    if 'hi_thresh' not in mask_kwargs:
        mask_kwargs['hi_thresh'] = 4.0

    if 'hi_nchan' not in mask_kwargs:
        mask_kwargs['hi_nchan'] = 2

    if 'lo_thresh' not in mask_kwargs:
        mask_kwargs['lo_thresh'] = 2.0

    if 'lo_nchan' not in mask_kwargs:
        mask_kwargs['lo_nchan'] = 2

    if 'min_pix' not in mask_kwargs:
        mask_kwargs['min_pix'] = None

    if 'min_area' not in mask_kwargs:
        mask_kwargs['min_area'] = None

    if 'grow_xy' not in mask_kwargs:
        mask_kwargs['grow_xy'] = 0

    if 'grow_v' not in mask_kwargs:
        mask_kwargs['grow_v'] = 0

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Create the mask
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    mask = cprops_mask(cube.filled_data[:].value, rms.filled_data[:].value, 
                       **mask_kwargs)

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Write to disk and return
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # In this case can avoid a recast
    if not return_spectral_cube and (outfile is None):
        return(mask)

    # Recast from numpy array to spectral cube    
    mask = SpectralCube(mask*1.0, wcs=cube.wcs, header=cube.header)
    
    # Write to disk, if desired
    if outfile is not None:        
        mask.write(outfile, overwrite=overwrite)
        
    if return_spectral_cube:
        return(mask)
    else:
        return(mask.filled_data[:].value)

def recipe_phangs_broad_mask(
    template_mask, outfile=None, list_of_masks = [],
    grow_xy = None, grow_v = None):
    """
    Task to create the PHANGS-style "broad" masks from the combination
    of a set of other masks. Optionally also grow the mask at the end.

    Parameters:

    -----------

    template_mask : string or SpectralCube

        The original mask that holds the target WCS. The other masks
        will be reprojected onto this one. This mask is included in
        the final output.

    Keywords:
    ---------

    outfile : string

        Filename where the mask will be written. The mask is also
        returned.

    list_of_masks : list

        List of masks or spectral cubes. These will be reprojected
        onto the template mask and then combined via logical or to
        form the final mask.

    grow_xy : int    
        
        Number of spatial dilations of the mask.

    grow_v : int

        Number of spectralwise dilations of the mask.

    """
    
    pass
