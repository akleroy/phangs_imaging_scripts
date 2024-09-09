# CASA imports

# This is a huge pain. Check that it works correctly by running

# casapy-XYZ -c casaStuff.py --nologger

# with XYZ each relevant version.

# AKL - checked with 6.5, 6.4, 6.3, 6.2.1, 5.8, 5.7, 5.6.1, 5.4, 5.3, 5.1.1, 5.0, 4.7.2, 4.5.3

# Obtain a version tuple (note the syntax change from < 6 to > 6)

from packaging import version

try:
    from taskinit import *
except ModuleNotFoundError:
    pass

if ('casa' in locals()) or ('casa' in globals()):
    casa_version = tuple(map(int, casa['build']['version'].replace('-','.').split('.')[0:3])) # tested CASA 4, 5
    casa_version_str = '.'.join([str(casa_version_no) for casa_version_no in casa_version])
else:
    # This works in CASA 6 where the casatools has a version attribute
    import casatools
    casa_version = (casatools.version()[0], casatools.version()[1], casatools.version()[2])
    casa_version_str = '.'.join([str(casa_version_no) for casa_version_no in casa_version])

print("CASA version: ", casa_version_str)

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
    from gencal import gencal
    from sdcal import sdcal
    from sdbaseline import sdbaseline
    from sdimaging import sdimaging

    from recipes.almahelpers import tsysspwmap


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
    from taskinit import metool
    from taskinit import qatool

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
    from casatools import (table, image, imager, msmetadata,
                           synthesisimager, synthesisutils, regionmanager,
                           measures, quanta)

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
                           visstat,
                           importasdm,
                           listobs,
                           flagcmd,
                           gencal,
                           sdcal,
                           sdbaseline,
                           sdimaging,
                           tsdimaging,
                           visstat,
                           rmtables)

    from casatasks.private import sdint_helper

    # from recipes.almahelpers import tsysspwmap

    # TODO: For now, uvcontsub doesn't work as we want it in newer CASA versions, fall back to old version
    if version.parse(casa_version_str) >= version.parse('6.5.2'):
        from casatasks import uvcontsub_old as uvcontsub

    # sdintimaging imports

    from casatasks.private import sdint_helper

    from .taskSDIntImaging import sdintimaging
    # from casatasks import sdintimaging

    # singledish processing imports
    #   see some documents at
    #   - https://casadocs.readthedocs.io/en/stable/api/casatasks.html?highlight=sdcal#single-dish
    #   - https://casadocs.readthedocs.io/en/stable/notebooks/synthesis_calibration.html?highlight=recipes

    iatool = image
    rgtool = regionmanager
    imtool = imager
    msmdtool = msmetadata
    tbtool = table
    metool = measures
    qatool = quanta

    from casatasks import (importasdm, listobs, flagcmd, flagdata)

    import almatasks

    import casaplotms

    # Depending on version, imports can be different
    try:
        plotms = casaplotms.gotasks.plotms.plotms
    except AttributeError:
        plotms = casaplotms.plotms

    import casaviewer

    try:
        viewer = casaviewer.gotasks.imview.imview
    except AttributeError:
        viewer = casaviewer.imview

    import casashell

    try:
        gencal = casashell.private.gencal.gencal
    except AttributeError:
        from casatasks import gencal

    try:
        plotbandpass = casashell.private.plotbandpass.plotbandpass
    except AttributeError:
        from casatasks import plotbandpass

    try:
        sdbaseline = casashell.private.sdbaseline.sdbaseline
    except AttributeError:
        from casatasks import sdbaseline

    try:
        sdimaging = casashell.private.sdimaging.sdimaging
    except AttributeError:
        from casatasks import sdimaging

    try:
        sdcal = casashell.private.sdcal.sdcal
    except AttributeError:
        from casatasks import sdcal

# sdintimaging import

if (casa_version[0] >= 5) and (casa_version[1] >= 7):
    from sdintimaging import sdintimaging

# tsysspwmap import

if casa_version[0] <= 5:
    from recipes.almahelpers import tsysspwmap
else:
    #from almatasks.private.almahelpers import tsysspwmap
    from almahelpers_localcopy import tsysspwmap
