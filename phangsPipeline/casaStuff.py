# CASA imports
try:
    from taskinit import *

    # Import specific CASA tasks. Not all of these are used by this
    # package, so this could be pared in the future.
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

    # version tuple

    casa_version = tuple(map(int, casa['build']['version'].replace('-','.').split('.')[0:3])) # tested CASA 4, 5
    
    # singledish processing imports
    
    if casa_version < (5, 0): # for singledish processing with precasa5
        from sdsave import sdsave
        from sdlist import sdlist
        from sdcal2 import sdcal2
        from sdscale import sdscale
        from sdplot import sdplot

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
    from recipes.almahelpers import tsysspwmap
    
    # sdintimaging imports

    if casa_version >= (5, 7):
        from sdintimaging import sdintimaging

except (ImportError, ModuleNotFoundError):

    # This is for CASA6

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
    from .taskSDIntImaging import sdintimaging

    # singledish processing imports
    # see some documents at 
    # - https://casadocs.readthedocs.io/en/stable/api/casatasks.html?highlight=sdcal#single-dish
    # - https://casadocs.readthedocs.io/en/stable/notebooks/synthesis_calibration.html?highlight=recipes
    
    iatool = image
    rgtool = regionmanager
    imtool = imager
    msmdtool = msmetadata
    tbtool = table

    import almatasks
    from almatasks.private.almahelpers import tsysspwmap
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

    # version tuple
    
    casa_version = casatools.version()




