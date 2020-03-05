import scDerivativeRoutines as scdr
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.io import fits
import numpy as np
from radio_beam import Beam
from scMaskingRoutines import recipe_phangs_mask as phangs_mask
from scMaskingRoutines import recipe_phangs_noise as phangs_noise


def _nicestr(quantity):
    if quantity.value == int(quantity.value):
        return(str(int(quantity.value))+' '+str(quantity.unit))
    else:
        return(str(quantity))

def _func_n_kwargs(derivative_name):
    if derivative_name == 'mom0':
        func = scdr.write_moment0
        kwargs ={'unit': u.K * u.km / u.s}
        return(func, kwargs)
    elif derivative_name == 'mom1':
        func = scdr.write_moment1
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif derivative_name == 'mom2':
        func = scdr.write_moment2
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif derivative_name == 'ew':
        func = scdr.write_ew
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif derivative_name == 'vquad':
        func = scdr.write_vquad
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif derivative_name == 'vpeak':
        func = scdr.write_vmax
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif derivative_name == 'tpeak':
        func = scdr.write_tmax
        kwargs = {'unit': u.K}
        return(func, kwargs)
    elif derivative_name == 'mom1hybrid':
        func = scdr.write_moment1_hybrid
        kwargs = {'unit': u.K}
        return(func, kwargs)
    
    
def moment_generator(cubefile,
                     root_name='',
                     mask=None,
                     rms=None,
                     rms_name='noise',
                     mask_name='signalmask',
                     derivatives=['mom0','mom1','mom2',
                                  'ew','vquad',
                                  'tpeak', 'vpeak'],
                     angular_resolution=None,
                     velocity_resolution=None,
                     linear_resolution=None,
                     distance=None,
                     generate_mask=True,
                     mask_kwargs=None,
                     noise_kwargs=None,
                     generate_noise=False,
                     channel_correlation=None):

    if type(cubefile) is SpectralCube:
        cube = cubefile
    else:
        cube = SpectralCube.read(cubefile)

    # We will be unit agnostic later
    cube = cube.to(u.K)
    
    if angular_resolution is not None and linear_resolution is not None:
        logger.error('Only one of angular_resolution or ',
                     'linear_resolution can be set')
    
    angres_name = ''
    if angular_resolution is not None:
        if type(angular_resolution) is str:
            angular_resolution = u.Quantity(angular_resolution)
        beam = Beam(major=angular_resolution,
                    minor=angular_resolution,
                    pa=0 * u.deg)
        cube = cube.convolve_to(beam)
        angres_name = '_' + _nicestr(angular_resolution).replace(
            ' ', '').replace('.', 'p').replace('/', '')

    linres_name = ''
    if linear_resolution is not None and distance is not None:
        if type(distance) is str:
            distance = u.Quantity(distance)
        if type(linear_resolution) is str:
            linear_resolution = u.Quantity(linear_resolution)
        angular_resolution = (linear_resolution / distance * u.rad).to(u.arcsec)
        beam = Beam(major=angular_resolution,
                    minor=angular_resolution,
                    pa=0 * u.deg)
        cube = cube.convolve_to(beam)
        linres_name = '_' + _nicestr(linear_resolution).replace(
            ' ', '').replace('.', 'p').replace('/', '')

    velres_name = ''
    if velocity_resolution is not None:
        if type(velocity_resolution) is str:
            velocity_resolution = u.Quantity(velocity_resolution)
        from astropy.convolution import Box1DKernel
        dv = scdr.channel_width(cube)
        nChan = (velocity_resolution / dv).to(u.dimensionless_unscaled).value
        if nChan > 1:
            cube = cube.spectral_smooth(Box1DKernel(nChan))
        velres_name = '_' + _nicestr(velocity_resolution).replace(
            ' ', '').replace('.', 'p').replace('/', '')

    if angres_name != '' or linres_name != '' or velres_name != '':
        cube.write(root_name + angres_name 
                   + linres_name + velres_name + '.fits',
                   overwrite=True)
        generate_mask = True

    if mask is not None:
        if type(mask) is str:
            mask_hdu = fits.open(mask)
            mask = np.array(mask_hdu[0].data, dtype=np.bool)
        cube = cube.with_mask(mask, inherit_mask=False)
    elif generate_mask:
        cube, rms = phangs_mask(cube, mask_kwargs=mask_kwargs,
                                noise_kwargs=noise_kwargs,
                                return_rms=True)
        m = SpectralCube(cube.mask.include().astype(np.uint8),
                         wcs=cube.wcs,
                         header=cube.header)
        m.write(root_name + '_' + mask_name
                + angres_name + linres_name + velres_name
                + '.fits', overwrite=True)
        rms.write(root_name + '_' + rms_name +
                  angres_name + linres_name + velres_name +
                  '.fits', overwrite=True)

    if rms is not None:
        if type(rms) is str:
            rms = SpectralCube.read(rms)

    elif generate_noise:
        rms = phangs_noise(cube, noise_kwargs=noise_kwargs,
                           return_spectral_cube=True)
        rms.write(root_name + '_noise'
                  + angres_name
                  + linres_name
                  + velres_name 
                  + '.fits', overwrite=True)
    for thisderivative in derivatives:
        func, prodkwargs = _func_n_kwargs(thisderivative)
        derivativefile = (root_name + '_' + thisderivative 
                       + angres_name 
                       + linres_name
                       + velres_name
                       + '.fits')
        if rms is not None:
            errorfile = derivativefile.replace('.fits', '_error.fits')
        else:
            errorfile = None
        func(cube, outfile=derivativefile,
             errorfile=errorfile,
             rms=rms,
             channel_correlation=channel_correlation,
             **prodkwargs)
    

