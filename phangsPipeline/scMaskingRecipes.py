from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.wcs import wcs
import numpy as np
import warnings

import logging
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
        logging.log('Unrecognized input type')
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
        logging.log('Unrecognized input type')
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

