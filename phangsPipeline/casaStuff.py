
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

    # Also get a simple two-digit version string
    #simple_version = '.'.join((casa['version'].split('.'))[0:2])
    simple_version = '{}.{}'.format(casa_version.major, casa_version.minor)

else:
    
    # Import taskinit tools
    from casatools import table as casatools_table
    tb = casatools_table()
    from casatools import image as casatools_image
    iatool = casatools_image()
    from casatools import imager as casatools_imager
    imtool = casatools_imager()

    from casatasks import casalog

    # Import specific CASA tasks. Not all of these are used by this
    # package, so this could be pared in the future.
    from casatasks import (concat, exportfits, feather, flagdata, imhead, immath, impbcor, importfits, 
        imrebin, imregrid, imsmooth, imstat, imsubimage, imtrans, imval, makemask, mstransform, 
        split, statwt, tclean, uvcontsub, visstat)

    # Also get a simple two-digit version string
    simple_version = '{}.{}'.format(casa_version.major, casa_version.minor)


