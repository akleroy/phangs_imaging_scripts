import scipy.ndimage.morphology as morph
import scipy.ndimage as nd
from scipy.signal import savgol_filter
import numpy as np
from astropy.stats import mad_std
from astropy.convolution import convolve, Gaussian2DKernel
import scipy.stats as ss
# from pipelineVersion import version as pipeVer

mad_to_std_fac = 1.482602218505602

def mad_zero_centered(data, mask=None):
    if mask is None:
        sig_false = ss.norm.isf(0.5 / data.size)
        mad1 = mad_to_std_fac * np.abs(np.median(data[data < 0]))
        mad2 = mad_to_std_fac * np.abs(np.median(np.abs(data[data <
                                                        (sig_false * mad1)])))
    else:
        nData = mask.sum()
        sig_false = ss.norm.isf(0.5 / nData)
        mad1 = mad_to_std_fac * np.abs(np.median(data[(data < 0) * mask]))
        mad2 = mad_to_std_fac * np.abs(np.median(np.abs(data[(data <
                                                        (sig_false * mad1)) 
                                                        * mask])))
    return(mad2)


def simple_mask(data, noise, hi_thresh=5, hi_nchan=2,
                lo_thresh=5, lo_nchan=2,
                min_pix=None, min_area=None,
                grow_xy=None, grow_v=None, invert=False):
    """
    Mask generation given data and estimate of the local noise
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
        detection. Default: 5
    lo_nchan : int
        Number of consecutive channels at lo_thresh required for a detection
        if connected to a hi_thresh detection. Default: 2
    min_pix : int
        Number of pixels required for a detection.  Default: None
    min_area : int
        Minimum number of pixels required in area projection for a detection.
        Default: None
    grow_xy : int
        NotImplemented
    grow_v : int
        NotImplemented
    invert : bool
        If True, invert the data before applying masking.
        Used for assessing the number of false positives given masking
        criteria. Default: False.
    """

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
    # revserve_indices like functionality here
    for hit in good_regions:
        mask[regions == hit] = True
    if grow_xy is not None:
        struct = morph.iterate_structure(morph.generate_binary_structure(2, 1),
                                         grow_xy)
        mask = morph.binary_dilation(mask, struct[np.newaxis, :, :])
    if grow_v is not None:
        struct = np.ones(grow_v, dtype=np.bool)
        mask = morph.binary_dilation(mask, struct[:, np.newaxis, np.newaxis])
    return(mask)


def noise_cube(data, mask=None, box=None, spec_box=None,
               nThresh=30, iterations=1,
               bandpass_smooth_window=None,
               bandpass_smooth_order=3):
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
        Boolean array with False indicating where data can be 
        used in the noise estimate. (i.e., True is Signal)
    
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
    
    bandpass_smooth_window : int
        Number of channels used in bandpass smoothing kernel.  Defaults to 
        nChan / 4 where nChan number of channels.  Set to zero to suppress 
        smoothing. Uses Savitzky-Golay smoothing
        
    bandpass_smooth_order : int
        Polynomial order used in smoothing kernel.  Defaults to 3.
    
    """
    noisemask = np.isfinite(data)
    if mask is not None:
        noisemask[mask] = False
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

    if bandpass_smooth_window is None:
        bandpass_smooth_window = 2 * (data.shape[0] // 8) + 1
        
    for i in np.arange(iterations):
        noise_map = np.zeros(data.shape[1:]) + np.nan
        noise_spec = np.zeros(data.shape[0]) + np.nan
        xx = np.arange(data.shape[2])
        yy = np.arange(data.shape[1])
        zz = np.arange(data.shape[0])
        for x in xx[boxr::step]:
            for y in yy[boxr::step]:
                spec = data[:, (y-boxr):(y+boxr+1),
                            (x-boxr):(x+boxr+1)]
                spec_mask = noisemask[:, (y-boxr):(y+boxr+1),
                                      (x-boxr):(x+boxr+1)]
                if np.sum(spec_mask) > nThresh:
                    noise_map[y, x] = mad_zero_centered(spec, mask=spec_mask)

        if boxr > 0:
            data_footprint = np.any(np.isfinite(data), axis=0)
            kernel = Gaussian2DKernel(box / np.sqrt(8 * np.log(2)))
            wt_map = np.isfinite(noise_map).astype(np.float)
            noise_map[np.isnan(noise_map)] = 0.0
            noise_map = convolve(noise_map, kernel)
            wt_map = convolve(wt_map, kernel)
            noise_map /= wt_map
            noise_map[~data_footprint] = np.nan

        for z in zz:
            lowz = np.clip(z - boxv, 0, data.shape[0])
            hiz = np.clip(z + boxv + 1, 0, data.shape[0])
            plane = data[lowz:hiz, :, :] / noise_map[np.newaxis, :, :]
            plane_mask = noisemask[lowz:hiz, :, :]
            noise_spec[z] = mad_zero_centered(plane, mask=plane_mask)
        # Smooth spectral shape
        if bandpass_smooth_window > 0:
            noise_spec = savgol_filter(noise_spec,
                                       bandpass_smooth_window,
                                       bandpass_smooth_order)
        noise_spec /= np.nanmedian(noise_spec)
        noise_cube = np.ones_like(data)
        noise_cube *= (noise_map[np.newaxis, :]
                       * noise_spec[:, np.newaxis, np.newaxis])
        if iterations == 1:
            return(noise_cube)
        else:
            data = data / noise_cube
            noise_cube_out *= noise_cube
    return(noise_cube_out)
