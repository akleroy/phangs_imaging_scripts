import scProductRoutines as scpr
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.io import fits
import numpy as np
from radio_beam import Beam


product_dict = {'mom0': None}


def moment_generator(cubefile,
                     mask_file=None,
                     rms_file=None,
                     products=['mom0','mom1','mom2',
                               'ew','vquad',
                               'tpeak', 'tpeak',
                               'vpeak','vpeak',
                               'mom1'],
                     angular_resolution=None,
                     velocity_resolution=None,
                     generate_mask=False,
                     masking_keywords=None):
    cube = SpectralCube.read(cubefile)
    cube = cube.to(u.K)
    
    if resolution is not None:
        beam = Beam(major=resolution,
                    minor=resolution,
                    pa=0 * u.deg)
        cube = cube.convolve_to(beam)

    if mask_file is not None:
        mask_hdu = fits.open(mask_file)
        mask = np.array(mask_hdu[0].data, dtype=np.bool)
        cube = cube.with_mask(mask, inherit_mask=False)
    elif generate_mask:
        mask = 
    if rms_file is not None:
        rms = SpectralCube.read(rms_file)

    raise NotImplementedError    
    

