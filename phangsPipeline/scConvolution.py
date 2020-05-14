from spectral_cube import SpectralCube
from radio_beam import Beam

import astropy.units as u
from astropy.io import fits
from astropy.convolution import Box1DKernel

import numpy as np

def smoothCube(
    cubefile,
    outfile=None,
    angular_resolution=None,
    linear_resolution=None,
    distance=None,
    velocity_resolution=None,
    ):
    """
    Smooth an input cube to coarser angular or spectral
    resolution. This only lightly wraps spectral cube and some of the
    error checking is left to that.
    """

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Error checking
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # Require a valid cube or map input
    if type(cubefile) is SpectralCube:
        cube = cubefile
    elif type(cubefile) == type("hello"):
        cube = SpectralCube.read(cubefile)
    else:
        logger.error("Input must be a SpectralCube object or a filename.")

    # Check that only one target scale is set
    if (angular_resolution is not None) and (linear_resolution is not None):
        logger.error('Only one of angular_resolution or ',
                     'linear_resolution can be set')
        return(None)

    # Work out the target angular resolution
    if angular_resolution is not None:
        if type(angular_resolution) is str:
            angular_resolution = u.Quantity(angular_resolution)

    if linear_resolution is not None:
        if distance is None:
            logger.error('Convolution to linear resolution requires a distance.')
            return(None)

        if type(distance) is str:
            distance = u.Quantity(distance)
        if type(linear_resolution) is str:
            linear_resolution = u.Quantity(linear_resolution)
        angular_resolution = (linear_resolution / distance * u.rad).to(u.arcsec)

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Convolution to coarser beam
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
    if angular_resolution is not None:
        beam = Beam(major=angular_resolution,
                    minor=angular_resolution,
                    pa=0 * u.deg)
        cube = cube.convolve_to(beam)

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Spectral convolution
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # This is only a boxcar smooth right now and does not downsample
    # or update the header.

    if velocity_resolution is not None:
        if type(velocity_resolution) is str:
            velocity_resolution = u.Quantity(velocity_resolution)

        dv = scdr.channel_width(cube)
        nChan = (velocity_resolution / dv).to(u.dimensionless_unscaled).value
        if nChan > 1:
            cube = cube.spectral_smooth(Box1DKernel(nChan))

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Write or return as requested
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    if outfile is not None:
        cube.write(outfile, overwrite=True)

    return(cube)
