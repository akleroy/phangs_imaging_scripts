import scProductRoutines as scpr
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.io import fits
import numpy as np
from radio_beam import Beam
from scMaskingRecipes import phangs_mask

product_dict = {'mom0': None}


def moment_generator(cubefile,
                     root_name='',
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
                     masking_kwargs=None):

    cube = SpectralCube.read(cubefile)
    cube = cube.to(u.K)
    
    angres_name = None
    if angular_resolution is not None:
        beam = Beam(major=angular_resolution,
                    minor=angular_resolution,
                    pa=0 * u.deg)
        cube = cube.convolve_to(beam)
        angres_name = '_' + str(angular_resolution).replace(
            ' ', '').replace('.', 'p').replace('/', '')

    velres_name = None
    if velocity_resolution is not None:
        from astropy.convolution import Box1DKernel
        dv = scpr.channel_width(cube)
        nChan = (velocity_resolution / dv).to(u.dimensionless_unscaled).value
        if nChan > 1:
            cube = cube.spectral_smooth(Box1DKernel(nChan))
        velres_name = '_' + str(velocity_resolution).replace(
            ' ', '').replace('.', 'p').replace('/', '')
        
    if mask_file is not None:
        mask_hdu = fits.open(mask_file)
        mask = np.array(mask_hdu[0].data, dtype=np.bool)
        cube = cube.with_mask(mask, inherit_mask=False)
    elif generate_mask:
        cube, rms = phangs_mask(cube, return_rms=True)

    if rms_file is not None:
        rms = SpectralCube.read(rms_file)

    raise NotImplementedError    
    

