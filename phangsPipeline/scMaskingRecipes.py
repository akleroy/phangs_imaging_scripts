from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.wcs import wcs
import numpy as np

def hybridize_mask(hires_filename, lores_filename, order='bilinear',
                   return_cube=False):
    hires_hdulist = fits.open(hires_filename)
    lores_hdulist = fits.open(lores_filename)
    
    hires = SpectralCube(np.array(hires_hdulist[0].data, dtype=np.bool),
                         wcs=wcs.WCS(hires_hdulist[0].header),
                         header=hires_hdulist[0].header)
    lores = SpectralCube(np.array(lores_hdulist[0].data, dtype=np.bool),
                         wcs=wcs.WCS(lores_hdulist[0].header),
                         header=lores_hdulist[0].header)
    lores = lores.reproject(hires.header, order=order)
    lores = lores.spectral_interpolate(hires.spectral_axis)
    mask = np.logical_or(np.array(hires_hdulist[0].data, dtype=np.bool),
                         np.array(lores.filled_data[:].value, dtype=np.bool))
    if return_cube:
        mask = SpectralCube(mask, wcs=wcs.WCS(hires_hdulist[0].header),
                            header=hires_hdulist[0].header,
                            meta={'BUNIT':' ', 'BTYPE':'Mask'})
    return(mask)
