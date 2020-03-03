import scProductRoutines as scpr
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.io import fits
import numpy as np
from radio_beam import Beam
from scMaskingRecipes import phangs_mask, phangs_noise

def _func_n_kwargs(productname):
    if productname == 'mom0':
        func = scpr.write_moment0
        kwargs ={'unit': u.K * u.km / u.s}
        return(func, kwargs)
    elif productname == 'mom1':
        func = scpr.write_moment1
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif productname == 'mom2':
        func = scpr.write_moment2
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif productname == 'ew':
        func = scpr.write_ew
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif productname == 'vquad':
        func = scpr.write_vquad
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif productname == 'vpeak':
        func = scpr.write_vmax
        kwargs = {'unit': u.km / u.s}
        return(func, kwargs)
    elif productname == 'tpeak':
        func = scpr.write_tmax
        kwargs = {'unit': u.K}
        return(func, kwargs)
    
    
def moment_generator(cubefile,
                     root_name='',
                     mask_file=None,
                     rms_file=None,
                     products=['mom0','mom1','mom2',
                               'ew','vquad',
                               'tpeak', 'vpeak'],
                     angular_resolution=None,
                     velocity_resolution=None,
                     generate_mask=True,
                     mask_kwargs=None,
                     noise_kwargs=None,
                     generate_noise=False,
                     channel_correlation=None):

    cube = SpectralCube.read(cubefile)
    cube = cube.to(u.K)
    
    angres_name = ''
    if angular_resolution is not None:
        if type(angular_resolution) is str:
            angular_resolution = u.Quantity(angular_resolution)
        beam = Beam(major=angular_resolution,
                    minor=angular_resolution,
                    pa=0 * u.deg)
        cube = cube.convolve_to(beam)
        angres_name = '_' + str(angular_resolution).replace(
            ' ', '').replace('.', 'p').replace('/', '')

    velres_name = ''
    if velocity_resolution is not None:
        if type(velocity_resolution) is str:
            velocity_resolution = u.Quantity(velocity_resolution)
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
        cube, rms = phangs_mask(cube,mask_kwargs=mask_kwargs,
                                noise_kwargs=noise_kwargs,
                                return_rms=True)
        m = SpectralCube(cube.mask.include().astype(np.uint8),
                         wcs=cube.wcs,
                         header=cube.header)
        m.write(root_name + '_signalmask' 
                + angres_name + velres_name 
                + '.fits', overwrite=True)
        rms.write(root_name + '_noise' +
                  angres_name + velres_name +
                  '.fits', overwrite=True)

    if rms_file is not None:
        rms = SpectralCube.read(rms_file)
    elif generate_noise:
        rms = phangs_noise(cube, noise_kwargs=noise_kwargs)
        rms.write(root_name + '_noise' + 
                  angres_name + velres_name + 
                  '.fits', overwrite=True)
    for thisproduct in products:
        func, prodkwargs = _func_n_kwargs(thisproduct)
        productfile = (root_name + '_' + thisproduct 
                       + angres_name 
                       + velres_name
                       + '.fits')
        if rms is not None:
            errorfile = productfile.replace('.fits', '_error.fits')
        else:
            errorfile = None
        func(cube, outfile=productfile,
             errorfile=errorfile,
             rms=rms,
             channel_correlation=channel_correlation,
             **prodkwargs)
    

