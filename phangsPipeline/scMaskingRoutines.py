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

np.seterr(divide='ignore', invalid='ignore')

mad_to_std_fac = 1.482602218505602

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def mad_zero_centered(data, mask=None):
    where_data_valid = np.logical_and(~np.isnan(data), np.isfinite(data))
    if mask is None:
        sig_false = ss.norm.isf(0.5 / data.size)
        data_lt_zero = np.less(data, 0,
                               where=where_data_valid,
                               out=np.full(where_data_valid.shape,
                                           False, dtype=bool))
        mad1 = mad_to_std_fac * np.abs(np.median(data[data_lt_zero]))
    
        data_lt_mad1 = np.less(data, sig_false * mad1,
                               where=where_data_valid,
                               out=np.full(where_data_valid.shape,
                                           False, dtype=bool))
        mad2 = mad_to_std_fac * np.abs(np.median(np.abs(data[data_lt_mad1])))
    else:
        nData = mask.sum()
        if nData == 0:
            logger.info('No data in mask. Returning NaN, which will now break things')
            return(np.nan)
        sig_false = ss.norm.isf(0.5 / nData)
        data_lt_zero = np.logical_and(np.less(data, 0, 
                                      where=where_data_valid,
                                      out=np.full(where_data_valid.shape,
                                                  False, dtype=bool)), mask)
        mad1 = mad_to_std_fac * np.abs(np.median(data[data_lt_zero]))
        data_lt_mad1 = np.logical_and(np.less(data, (sig_false * mad1),
                                      where=where_data_valid,
                                      out=np.full(where_data_valid.shape,
                                                  False, dtype=bool)), mask)
        mad2 = mad_to_std_fac * np.abs(np.median(np.abs(data[data_lt_mad1])))
    return(mad2)


def simple_mask(data, noise, hi_thresh=5, hi_nchan=2,
                lo_thresh=None, lo_nchan=2,
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
    if lo_thresh is None:
        lo_thresh = hi_thresh
    
    hi_mask = np.greater_equal(signif, hi_thresh,
                               where=(~np.isnan(signif)),
                               out=np.full(signif.shape, False, dtype=bool))
    lo_mask = np.greater_equal(signif, lo_thresh,
                               where=(~np.isnan(signif)),
                               out=np.full(signif.shape, False, dtype=bool))
    
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
        boxr = int(box // 2)
    if spec_box is not None:
        spec_step = np.floor(box/2).astype(np.int)
        boxv = int(spec_box // 2)
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
            kernel = savgol_coeffs(int(bandpass_smooth_window),
                                   int(bandpass_smooth_order))
            noise_spec = convolve(noise_spec, kernel, boundary='extend')
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


def recipe_phangs_noise(cube, noise_kwargs=None,
                        return_spectral_cube=False):
    if noise_kwargs is None:
        pixels_per_beam = cube.pixels_per_beam
        # pixels_per_beam = (cube.beam.sr
        #                    / wcs.utils.proj_plane_pixel_area(cube.wcs)
        #                    / u.deg**2).to(u.dimensionless_unscaled).value
        box = np.ceil(2.5 * pixels_per_beam**0.5)
        spectral_smooth = np.ceil(cube.shape[0] / 5) // 2 * 2 + 1
        # This looks for a non-trivial signal mask.
        if (np.sum(cube.mask.include())
                < np.sum(np.isfinite(cube.filled_data[:].value))):
            m = cube.mask.include()
        else:
            m = None
        rms = noise_cube(cube.filled_data[:].value,
                         mask=m,
                         box=box,
                         bandpass_smooth_window=spectral_smooth,
                         spec_box=5,
                         iterations=3)
    else:
        rms = noise_cube(cube.filled_data[:].value,
                         mask=cube.mask.include(),
                         **noise_kwargs)
    if return_spectral_cube:
        rms = SpectralCube(rms, wcs=cube.wcs, header=cube.header)
    return(rms)


def recipe_phangs_mask(cube,
                       mask_kwargs=None,
                       noise_kwargs=None,
                       return_rms=False):

    rms = recipe_phangs_noise(cube, noise_kwargs=noise_kwargs)

    if mask_kwargs is None:
        mask = simple_mask(cube.filled_data[:].value,
                           rms, hi_thresh=4, hi_nchan=2,
                           lo_thresh=2, lo_nchan=2)
    else:
        mask = simple_mask(cube.filled_data[:].value,
                           rms, **mask_kwargs)

    if return_rms:
        return(cube.with_mask(mask, inherit_mask=False),
               SpectralCube(rms, wcs=cube.wcs, header=cube.header))
    else:
        return(cube.with_mask(mask, inherit_mask=False))
