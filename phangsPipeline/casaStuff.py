
import sys
import numpy as np
from collections import namedtuple
VersionInfo = namedtuple('VersionInfo', ['major', 'minor', 'micro', 'patch']) # 'releaselevel', serial'
str2versioninfo = lambda x: VersionInfo(*(np.pad(np.array(str(x).replace('-','.').split('.')).astype(int), (0,4), mode='constant', constant_values=0)[0:4]))

def get_casa_version(
        return_python_version = False,
    ):
    """
    Get CASA version if we are in CASA environment. 
    
    If `return_python_version` is set to True, then we will also return the python version major number. 
    """
    casa_enabled = False
    casa_version = None
    python_version = None
    # 
    if not casa_enabled:
        try:
            import taskinit
            casa_enabled = True
            import casadef
            casa_version = str2versioninfo(casadef.casa_version)
            #python_version = [2]
        except ImportError:
            pass
    # 
    if not casa_enabled:
        try:
            from casatasks import tclean
            casa_enabled = True
            from casatools import utils as casatools_utils
            casa_version = str2versioninfo(casatools_utils.utils().version_string()) # see https://casadocs.readthedocs.io/en/stable/api/tt/casatools.utils.html#casatools.utils.utils.version
            #python_version = [3]
        except ImportError:
            pass
    # 
    if return_python_version:
        python_version = sys.version_info
        return casa_version, python_version
    else:
        return casa_version


# Get CASA VersionInfo

casa_version = get_casa_version()


# Import CASA modules

if casa_version.major <= 5:
    
    # For CASA <= 5.X.X
    
    # Import taskinit tools
    #   'aftool', 'at', 'attool', 'casa', 'casac', 'casaglobals', 'casalog', 'cbtool', 'cltool', 'coordsystool', 
    #   'cptool', 'cstool', 'cu', 'dctool', 'find_casa', 'fitool', 'fntool', 'gentools', 'iatool', 'imdtool', 
    #   'imtool', 'inspect', 'lmtool', 'me', 'metool', 'mptool', 'ms', 'msmdtool', 'mstool', 'mttool', 'os', 
    #   'pm', 'pmtool', 'potool', 'qa', 'qatool', 'rgtool', 'sbstool', 'sdms', 'sdmstool', 'sltool', 'smtool', 
    #   'stack_find', 'static_var', 'string', 'sys', 'tb', 'tbtool', 'tptool', 'utilstool', 'viewertool', 
    #   'vptool', 'write_history', 'xmlpath'
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
    from importasdm import importasdm
    from listobs import listobs
    from flagcmd import flagcmd
    from flagdata import flagdata
    
    try:
        # useful for singledish processing with precasa5
        from sdsave import sdsave
        from sdlist import sdlist
        from sdcal2 import sdcal2
        from sdscale import sdscale
        from sdplot import sdplot
    except:
        pass
    
    from plotms import plotms
    from viewer import viewer
    from gencal import gencal
    from plotbandpass import plotbandpass
    from sdbaseline import sdbaseline
    from sdimaging import sdimaging
    from sdcal import sdcal
    #from taskinit import msmdtool
    #from taskinit import tbtool
    from recipes.almahelpers import tsysspwmap
    
    # Also get a simple two-digit version string
    #simple_version = '.'.join((casa['version'].split('.'))[0:2])
    simple_version = '{}.{}'.format(casa_version.major, casa_version.minor)

else:
    
    # For CASA >= 6.X.X
    
    # Import taskinit tools
    from casatools import table as casatools_table
    tb = casatools_table()
    from casatools import image as casatools_image
    iatool = casatools_image
    from casatools import imager as casatools_imager
    imtool = casatools_imager
    
    from casatasks import casalog
    
    # Import specific CASA tasks. Not all of these are used by this
    # package, so this could be pared in the future.
    from casatasks import (concat, exportfits, feather, flagdata, flagcmd, imhead, immath, impbcor, importfits, 
        imrebin, imregrid, imsmooth, imstat, imsubimage, imtrans, imval, makemask, mstransform, 
        split, statwt, tclean, uvcontsub, visstat)
    
    # Import stuff for single dish: plotms, viewer, gencal, plotbandpass, sdbaseline, sdimaging, sdcal
    # some documentation can be found at https://casadocs.readthedocs.io/en/stable/api/casatasks.html?highlight=sdcal#single-dish
    #from casatasks.private import sdutil
    # some documentation can be found at https://casadocs.readthedocs.io/en/stable/notebooks/synthesis_calibration.html?highlight=recipes
    import almatasks
    from almatasks.private.almahelpers import tsysspwmap
    # searched "tsysspwmap" in casa subdirectories
    import casaplotms
    plotms = casaplotms.gotasks.plotms.plotms
    # checked type(plotms) in casa6
    import casaviewer
    viewer = casaviewer.gotasks.imview.imview
    # checked type(viewer) in casa6
    import casashell
    gencal = casashell.private.gencal.gencal
    plotbandpass = casashell.private.plotbandpass.plotbandpass
    sdbaseline = casashell.private.sdbaseline.sdbaseline
    sdimaging = casashell.private.sdimaging.sdimaging
    sdcal = casashell.private.sdcal.sdcal
    
    # Also get a simple two-digit version string
    simple_version = '{}.{}'.format(casa_version.major, casa_version.minor)


