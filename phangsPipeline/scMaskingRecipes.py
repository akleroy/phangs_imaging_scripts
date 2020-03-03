from astropy.io import fits
from spectral_cube import SpectralCube
from astropy import wcs
import numpy as np
import warnings
from scMaskingRoutines import simple_mask, noise_cube
import logging
import astropy.units as u

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def hybridize_mask(hires_in, lores_in, order='bilinear',
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
        lores_hdulist = fits.open(lores_filename)
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
    mask = np.logical_or(np.array(hires_hdulist[0].data, dtype=np.bool),
                         np.array(lores.filled_data[:].value, dtype=np.bool))
    if return_cube:
        mask = SpectralCube(mask, wcs=wcs.WCS(hires_hdulist[0].header),
                            header=hires_hdulist[0].header,
                            meta={'BUNIT':' ', 'BTYPE':'Mask'})
    return(mask)


def phangs_noise(cube, noise_kwargs=None,
                 return_spectral_cube=False):
    if noise_kwargs is None:
        pixels_per_beam = (cube.beam.sr
                        / wcs.utils.proj_plane_pixel_area(cube.wcs)
                        / u.deg**2).to(u.dimensionless_unscaled).value
        box = np.ceil(2.5 * pixels_per_beam**0.5)
        spectral_smooth = np.ceil(cube.shape[0] / 5)
        # This looks for a non-trivial signal mask.
        if (np.sum(cube.mask.include()) 
            < np.sum(np.isfinite(cube.unmasked_data[:].value))):
            m = cube.mask.include()
        else:   
            m = None
        rms = noise_cube(cube.unmasked_data[:].value,
                         mask=m,
                         box=box,
                         bandpass_smooth_window=spectral_smooth,
                         spec_box=5,
                         iterations=3)
    else:
        rms = noise_cube(cube.unmasked_data[:].value,
                         mask=cube.mask.include(),
                         **noise_kwargs)
    if return_spectral_cube:
        rms = SpectralCube(rms, wcs=cube.wcs, header=cube.header)
    return(rms)


def phangs_mask(cube,
                mask_kwargs=None,
                noise_kwargs=None, 
                return_rms=False):
    
    rms = phangs_noise(cube, noise_kwargs=noise_kwargs)
    
    if mask_kwargs is None:
        mask = simple_mask(cube.unmasked_data[:].value,
                           rms, hi_thresh=4, hi_nchan=2,
                           lo_thresh=2, lo_nchan=2)
    else:
        mask = simple_mask(cube.unmasked_data[:].value,
                           rms, **mask_kwargs)
        
    if return_rms:
        return(cube.with_mask(mask, inherit_mask=False),
               SpectralCube(rms, wcs=cube.wcs, header=cube.header))
    else:
        return(cube.with_mask(mask, inherit_mask=False))
