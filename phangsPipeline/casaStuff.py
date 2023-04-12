# CASA imports

# This is a huge pain. Check that it works correctly by running

# casapy-XYZ -c casaStuff.py --nologger

# with XYZ each relevant version.

# AKL - checked with 6.5, 6.4, 6.3, 6.2.1, 5.8, 5.7, 5.6.1, 5.4, 5.3, 5.1.1, 5.0, 4.7.2, 4.5.3

# Obtain a version tuple (note the syntax change from < 6 to > 6)

try:
    from taskinit import *
except ModuleNotFoundError:
    pass

if ('casa' in locals()) or ('casa' in globals()):
    casa_version = tuple(map(int, casa['build']['version'].replace('-','.').split('.')[0:3])) # tested CASA 4, 5
else:
    # This works in CASA 6 where the casatools has a version attribute
    import casatools
    casa_version = (casatools.version()[0], casatools.version()[1], casatools.version()[2])

print("CASA version: ", casa_version)

# Import specific CASA tasks. Not all of these are used by this
# package, so this could be pared in the future.
        
if casa_version[0] < 6:
    from taskinit import *
    
    from concat import concat
    from exportfits import exportfits
    from feather import feather
    from flagdata import flagdata
    from imhead import imhead
    from immath import immath
    from impbcor import impbcor
    from importfits import importfits
    from imrebin import imrebin
    from imregrid import imregrid
    from imsmooth import imsmooth
    from imstat import imstat
    from imsubimage import imsubimage
    from imtrans import imtrans
    from imval import imval
    from makemask import makemask
    from mstransform import mstransform
    from split import split
    from statwt import statwt
    from tclean import tclean
    from uvcontsub import uvcontsub
    from visstat import visstat

    from importasdm import importasdm
    from listobs import listobs
    from flagcmd import flagcmd
    from flagdata import flagdata
    from plotms import plotms
    from viewer import viewer
    from gencal import gencal
    from plotbandpass import plotbandpass
    from sdbaseline import sdbaseline
    from sdimaging import sdimaging
    from sdcal import sdcal
    from taskinit import msmdtool
    from taskinit import tbtool
        
# imports for singledish processing when CASA version < 5
    
if casa_version[0] < 5: 
    from sdsave import sdsave
    from sdlist import sdlist
    from sdcal2 import sdcal2
    from sdscale import sdscale
    from sdplot import sdplot

# Imports for CASA versions above 6
    
if casa_version[0] >= 6:

    import casatools
    from casatools import (table, image, imager, msmetadata, synthesisimager, synthesisutils, regionmanager)

    import casatasks
    from casatasks import (casalog,
                           concat,
                           exportfits,
                           feather,
                           flagdata,
                           imhead,
                           immath,
                           impbcor,
                           importfits,
                           imrebin,
                           imregrid,
                           imsmooth,
                           imstat,
                           imsubimage,
                           imtrans,
                           imval,
                           makemask,
                           mstransform,
                           split,
                           statwt,
                           tclean,
                           uvcontsub,
                           visstat)
    
    # sdintimaging imports
    
    from casatasks.private import sdint_helper

    #from .taskSDIntImaging import sdintimaging
    from casatasks import sdintimaging

    # singledish processing imports
    #   see some documents at 
    #   - https://casadocs.readthedocs.io/en/stable/api/casatasks.html?highlight=sdcal#single-dish
    #   - https://casadocs.readthedocs.io/en/stable/notebooks/synthesis_calibration.html?highlight=recipes
        
    iatool = image
    rgtool = regionmanager
    imtool = imager
    msmdtool = msmetadata
    tbtool = table
    
    from casatasks import (importasdm, listobs, flagcmd, flagdata)

    import almatasks
    
    import casaplotms
    plotms = casaplotms.gotasks.plotms.plotms
    import casaviewer
    viewer = casaviewer.gotasks.imview.imview
    import casashell
    gencal = casashell.private.gencal.gencal
    plotbandpass = casashell.private.plotbandpass.plotbandpass
    sdbaseline = casashell.private.sdbaseline.sdbaseline
    sdimaging = casashell.private.sdimaging.sdimaging
    sdcal = casashell.private.sdcal.sdcal

# sdintimaging import

if (casa_version[0] >= 5) and (casa_version[1] >= 7):
    from sdintimaging import sdintimaging

# tsysspwmap import
    
if casa_version[0] <= 5:
    from recipes.almahelpers import tsysspwmap
else:
    #from almatasks.private.almahelpers import tsysspwmap
    from almahelpers_localcopy import tsysspwmap
