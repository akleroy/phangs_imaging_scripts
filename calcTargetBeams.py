from spectral_cube import SpectralCube, VaryingResolutionSpectralCube
from radio_beam import Beam
import numpy as np
import astropy.units as u
import astropy.utils.console as console
import copy
import warnings

def calc_target_beam(incube):
    if isinstance(incube,str):
        cube = SpectralCube.read(incube)

    if isinstance(incube, VaryingResolutionSpectralCube):
        cube = incube

    if not isinstance(cube, VaryingResolutionSpectralCube):
        warnings.warn("No information about multiple beams")
        return(None)
    
    beams = cube.beams
    major_axes = np.array([bm.major.to(u.deg).value for bm in beams])
    target_beamsize = np.array(major_axes.max())
    target_beam = Beam(major=target_beamsize*u.deg,
                       minor=target_beamsize*u.deg,
                       pa=0.0*u.deg)
    print ("Target beam is :", target_beam)
    return target_beam

gals = {'ngc0628':'arcsec' \
            , 'ngc1672':'arcsec' \
            , 'ngc3351':'arcsec' \
            , 'ngc3627':'arcsec' \
            , 'ngc4254':'arcsec' \
            , 'ngc4303':'arcsec' \
            , 'ngc4321':'arcsec' \
            , 'ngc4535':'arcsec' \
            , 'ngc5068':'arcsec' \
            , 'ngc6744north':'arcsec' \
            , 'ngc6744south':'arcsec' \
            }

for this_gal in gals.keys():
    print this_gal
    target_beam = calc_target_beam('../release/v0p4/raw/'+this_gal+'_co21.fits')
    print target_beam
