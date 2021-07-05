# CASA imports
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

# Version

#try:
#    casa
#    if 'version' in casa:
#        simple_version = '.'.join((casa['version'].split('.'))[0:2])
