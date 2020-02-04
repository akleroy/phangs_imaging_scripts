from scipy.stats import norm
from spectral_cube import SpectralCube
import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
from scipy.special import erfcinv
from astropy.convolution import convolve, Gaussian2DKernel
import matplotlib.pyplot as plt  # for debug only
import scipy.ndimage.morphology as morph
import scipy.ndimage as nd
from scipy.signal import savgol_filter


def make_mask(data, noise, hi_thresh=5, hi_nchan=2,
              lo_thresh=5, lo_nchan=2,
              min_pix=None, min_area=None,
              grow_xy=None, grow_v=None, invert=False):
    signif = data / noise
    if invert:
        signif *= -1
    hi_mask = signif >= hi_thresh
    lo_mask = signif >= lo_thresh
    histr = np.ones(hi_nchan, dtype=np.bool)
    lostr = np.ones(lo_nchan, dtype=np.bool)
    hi_mask = morph.binary_opening(hi_mask, histr[:,
                                                  np.newaxis,
                                                  np.newaxis])
    lo_mask = morph.binary_opening(lo_mask, lostr[:,
                                                  np.newaxis,
                                                  np.newaxis])
    # Prune small volume regions
    if (min_pix is not None):
        regions, regct = nd.label(hi_mask)
        objslices = nd.find_objects(regions)
        for ii, thisslice in enumerate(objslices):
            subcube = regions[thisslice]
            pixct = np.sum(subcube == (ii+1))
            if pixct < min_pix:
                hi_mask[regions == (ii+1)] = False
    if (min_area is not None):
        regions, regct = nd.label(hi_mask)
        objslices = nd.find_objects(regions)
        for ii, thisslice in enumerate(objslices):
            subcube = regions[thisslice]
            area = np.sum(np.any(subcube == (ii+1), axis=0))
            if area < min_area:
                hi_mask[regions == (ii+1)] = False

    # Expand in to lo_mask
    regions, regct = nd.label(lo_mask)
    good_regions = np.unique(regions[hi_mask])
    mask = np.zeros_like(hi_mask, dtype=np.bool)
    # revserve_indices like functinoality here
    for hit in good_regions:
        mask[regions == hit] = True
    return(mask)


def noise_cube(data, mask=None, box=None, spec_box=None,
               nThresh=30, iterations=1):
    """
    Makes an empirical estimate of the noise in a cube assuming that it 
    is normally distributed.
    
    Parameters:
    -----------
    
    data : np.array
        Array of data (floats)
    
    Keywords:
    ---------
    
    mask : np.bool
        Boolean array with True indicating where data can be 
        used in the noise estimate.
    
    box : int
        Spatial size of the box over which noise is calculated (correlation 
        scale).  Default: no box
    
    spec_box : int
        Spectral size of the box overwhich the noise is calculated.  Default:
        no box
    
    nThresh : int
        Minimum number of date to be used in a noise estimate.
    
    iterations : int
        Number of times to iterate the noise solution to force Gaussian 
        statistics.  Default: no iterations.
    
    """
    if mask is None:
        mask = np.isfinite(data)
    step = 1
    boxr = step // 2
    if box is not None:
        step = np.floor(box/2.5).astype(np.int)
        boxr = step // 2
    if spec_box is not None:
        spec_step = np.floor(box/2).astype(np.int)
        boxv = spec_step // 2
    else:
        boxv = 0
    noise_cube_out = np.ones_like(data)
    for i in np.arange(iterations):
        noise_map = np.zeros(data.shape[1:]) + np.nan
        noise_spec = np.zeros(data.shape[0]) + np.nan
        sig_false = erfcinv(0.5/data.shape[0])
        sigma = 1.4826 * np.abs(np.median(data[data < 0]))
        xx = np.arange(data.shape[2])
        yy = np.arange(data.shape[1])
        zz = np.arange(data.shape[0])
        for x in xx[boxr::step]:
            for y in yy[boxr::step]:
                spec = data[:, (y-boxr):(y+boxr+1),
                            (x-boxr):(x+boxr+1)]
                spec_mask = mask[:, (y-boxr):(y+boxr+1),
                                 (x-boxr):(x+boxr+1)]
                if np.sum(np.isfinite(spec)) > nThresh:
                    sigma1 = 1.4826 * np.abs(np.median(spec[(spec < 0)
                                                            * (spec_mask)]))
                    noise_map[y, x] = mad_std(spec[(spec < (sigma1 * sig_false))
                                                   * (spec_mask)])
        if boxr > 0:
            data_footprint = np.any(np.isfinite(data), axis=0)
            kernel = Gaussian2DKernel(box / np.sqrt(8 * np.log(2)))
            wt_map = np.isfinite(noise_map).astype(np.float)
            # wt_map[~data_footprint] = np.nan
            noise_map[np.isnan(noise_map)] = 0.0
            # noise_map[~data_footprint] = np.nan
            noise_map = convolve(noise_map, kernel)
            wt_map = convolve(wt_map, kernel)
            noise_map /= wt_map
            noise_map[~data_footprint] = np.nan

        for z in zz:
            lowz = np.clip(z - boxv, 0, data.shape[0])
            hiz = np.clip(z + boxv + 1, 0, data.shape[0])
            plane = data[lowz:hiz, :, :] / noise_map[np.newaxis, :, :]
            plane_mask = mask[lowz:hiz, :, :]
            sig_false = erfcinv(0.5 / np.sum(np.isfinite(plane)))
            sigma1 = 1.4826 * np.abs(np.median(plane[(plane < 0)
                                                     * plane_mask]))
            noise_spec[z] = mad_std(plane[(plane < (sig_false * sigma1))
                                          * plane_mask])
        # Smooth spectral shape
        if boxv > 0:
            noise_spec = savgol_filter(noise_spec, 4 * boxv + 1, 3)
        noise_spec /= np.nanmedian(noise_spec)
        noise_cube = np.ones_like(data)
        noise_cube *= (noise_map[np.newaxis, :]
                       * noise_spec[:, np.newaxis, np.newaxis])
        if iterations == 1:
            return(noise_cube)
        else:
            data = data / noise_cube
            noise_cube_out *= noise_cube
            print('Iteration {0} complete'.format(i))
    return(noise_cube_out)
