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

    # Version

    simple_version = '.'.join((casa['version'].split('.'))[0:2])

except (ImportError, ModuleNotFoundError):

    # This is for CASA6

    from casatools import table, image, imager

    tb = table()

    iatool = image
    imtool = imager

    import casatools
    simple_version = casatools.version()

    from casatasks import (concat,
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
