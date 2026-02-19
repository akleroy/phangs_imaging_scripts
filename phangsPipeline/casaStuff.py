# CASA imports

# This is a huge pain. Check that it works correctly by running

# casapy-XYZ -c casaStuff.py --nologger

# with XYZ each relevant version.

# AKL - checked with 6.5, 6.4, 6.3, 6.2.1, 5.8, 5.7, 5.6.1, 5.4, 5.3, 5.1.1, 5.0, 4.7.2, 4.5.3
# TGW - checked with 6.7
try:
    from packaging import version
except ImportError:
    class version:
        def parse(self, vstr):
            return tuple(map(int, vstr.replace('-','.').split('.')[0:3]))

# Obtain a version tuple
if ("casa" in locals()) or ("casa" in globals()):
    casa_version = tuple(
        map(int, casa["build"]["version"].replace("-", ".").split(".")[0:3])
    )  # tested CASA 4, 5
    casa_version_str = ".".join(
        [str(casa_version_no) for casa_version_no in casa_version]
    )
else:
    # This works in CASA 6 where the casatools has a version attribute
    import casatools

    casa_version = (
        casatools.version()[0],
        casatools.version()[1],
        casatools.version()[2],
    )
    casa_version_str = ".".join(
        [str(casa_version_no) for casa_version_no in casa_version]
    )

print("CASA version: ", casa_version_str)

# Import specific CASA tasks

# Imports for CASA versions above 6
if casa_version[0] >= 6:

    import casatools
    from casatools import (
        table,
        image,
        imager,
        msmetadata,
        synthesisimager,
        synthesisutils,
        regionmanager,
        measures,
        quanta,
    )

    import casatasks
    from casatasks import (
        casalog,
        concat,
        exportfits,
        feather,
        flagcmd,
        flagdata,
        gencal,
        imhead,
        immath,
        impbcor,
        importasdm,
        importfits,
        imrebin,
        imregrid,
        imsmooth,
        imstat,
        imsubimage,
        imtrans,
        imval,
        listobs,
        makemask,
        mstransform,
        rmtables,
        sdbaseline,
        sdcal,
        split,
        statwt,
        tclean,
        tsdimaging,
        uvcontsub,
        visstat,
    )

    # sdintimaging imports
    from casatasks.private import sdint_helper

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

    import casaplotms

    # Depending on version, imports can be different
    try:
        plotms = casaplotms.gotasks.plotms.plotms
    except AttributeError:
        plotms = casaplotms.plotms

    try:
        import casaviewer
    except (ImportError, ModuleNotFoundError):
        casaviewer = None
        print("Could not import casaviewer")

    if casaviewer is not None:
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
        sdcal = casashell.private.sdcal.sdcal
    except AttributeError:
        from casatasks import sdcal

    # tsysspwmap import
    from almahelpers_localcopy import tsysspwmap

if casa_version[0] < 6:
    from taskinit import *

    from concat import concat
    from exportfits import exportfits
    from feather import feather
    from flagcmd import flagcmd
    from flagdata import flagdata
    from gencal import gencal
    from imhead import imhead
    from immath import immath
    from impbcor import impbcor
    from importasdm import importasdm
    from importfits import importfits
    from imrebin import imrebin
    from imregrid import imregrid
    from imsmooth import imsmooth
    from imstat import imstat
    from imsubimage import imsubimage
    from imtrans import imtrans
    from imval import imval
    from listobs import listobs
    from makemask import makemask
    from mstransform import mstransform
    from plotbandpass import plotbandpass
    from plotms import plotms
    from sdbaseline import sdbaseline
    from sdcal import sdcal
    from sdimaging import sdimaging
    from split import split
    from statwt import statwt
    from taskinit import metool
    from taskinit import msmdtool
    from taskinit import qatool
    from taskinit import tbtool
    from tclean import tclean
    from uvcontsub import uvcontsub
    from viewer import viewer
    from visstat import visstat

    from recipes.almahelpers import tsysspwmap

# Imports for singledish processing when CASA version < 5
if casa_version[0] < 5:
    from sdsave import sdsave
    from sdlist import sdlist
    from sdcal2 import sdcal2
    from sdscale import sdscale
    from sdplot import sdplot

# sdintimaging import
if (casa_version[0] >= 5) and (casa_version[1] >= 7):
    from .taskSDIntImaging import sdintimaging
