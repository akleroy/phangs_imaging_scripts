# Licensed under a MIT license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

from .casa_check import is_casa_installed
casa_enabled = is_casa_installed()

from .phangsLogger import setup_logger
from .handlerKeys import KeyHandler
from .handlerSingleDish import SingleDishHandler
from .handlerAlmaDownload import AlmaDownloadHandler
from .handlerVis import VisHandler
from .handlerPostprocess import PostProcessHandler
from .handlerDerived import DerivedHandler
from .handlerRelease import ReleaseHandler

if casa_enabled:
    from .handlerImaging import ImagingHandler

__all__ = ["setup_logger", "KeyHandler", "SingleDishHandler", "VisHandler", "PostProcessHandler", "DerivedHandler",
           "ReleaseHandler"]

if casa_enabled:
    __all__.append("ImagingHandler")

try:
    from .handlerAlmaDownload import AlmaDownloadHandler
    __all__.append("AlmaDownloadHandler")
except ImportError:
    pass
