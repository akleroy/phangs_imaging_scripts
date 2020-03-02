import scProductRoutines as scpr
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.io import fits
import numpy as np


product_dict = {'mom0':
    }


def moment_generator(cubefile,
                     mask_file=None,
                     rms_file=None,
                     products = ['mom0','mom1','mom2',
                                 'ew','vquad',
                                 'tpeak', 'tpeak_12p5kms',
                                 'vpeak','vpeak_12p5kms',
                                 'mom1xs']):
    cube = SpectralCube.read(cubefile)
    cube = cube.to(u.K)

    if mask_file is not None:
        mask_hdu = fits.open(mask_file)
        mask = np.array(mask_hdu[0].data, dtype=np.bool)
        cube = cube.with_mask(mask, inherit_mask=False)

    if rms_file is not None:
        rms = SpectralCube.read(rms_file)

    
    

