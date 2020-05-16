from spectral_cube import SpectralCube
from radio_beam import Beam

import astropy.units as u
from astropy.io import fits
from astropy.convolution import Box1DKernel
from astropy.convolution import convolve, convolve_fft

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

import numpy as np

def smooth_cube(
    incube=None,
    outfile=None,
    angular_resolution=None,
    linear_resolution=None,
    distance=None,
    velocity_resolution=None,
    nan_treatment='interpolate', # can also be 'fill'
    tol=None,
    overwrite=True
    ):
    """
    Smooth an input cube to coarser angular or spectral
    resolution. This lightly wraps spectral cube and some of the error
    checking is left to that.

    tol is a fraction. When the target beam is within tol of the
    original beam, we just copy.
    """

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Error checking
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # Require a valid cube or map input
    if type(incube) is SpectralCube:
        cube = incube
    elif type(incube) == type("hello"):
        cube = SpectralCube.read(incube)
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
        dist_mpc_val = float(distance.to(u.pc).value) / 1e6
        cube._header.append(('DIST_MPC',dist_mpc_val,'Used in convolution'))

    if tol is None:
        tol = 0.0

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Convolution to coarser beam
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
    if angular_resolution is not None:
        logger.info("... convolving from beam: "+str(cube.beam))
        target_beam = Beam(major=angular_resolution,
                           minor=angular_resolution,
                           pa=0 * u.deg)
        logger.info("... convolving to beam: "+str(target_beam))

        new_major = float(target_beam.major.to(u.arcsec).value)
        old_major = float(cube.beam.major.to(u.arcsec).value)        
        delta = (new_major-old_major)/old_major

        logger.info("... fractional change: "+str(delta))
        
        if delta > tol:
            logger.info("... proceeding with convolution.")
            cube = cube.convolve_to(target_beam,
                                    nan_treatment=nan_treatment)
        if np.abs(delta) < tol:
            logger.info("... current resolution meets tolerance.")
        if delta < -1.0*tol:
            logger.info("... resolution cannot be matched. Returning")
            return(None)

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
        cube.write(outfile, overwrite=overwrite)

    return(cube)
