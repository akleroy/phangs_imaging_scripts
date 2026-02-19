from importlib.metadata import version

# Ensure CASA is installed
from .casa_check import is_casa_installed

casa_enabled = is_casa_installed()
if not casa_enabled:
    raise ImportError("CASA not detected; pipeline is incorrectly installed")

from .handlerAlmaDownload import AlmaDownloadHandler
from .handlerDerived import DerivedHandler
from .handlerImaging import ImagingHandler
from .handlerImagingChunked import ImagingChunkedHandler
from .handlerKeys import KeyHandler
from .handlerPostprocess import PostProcessHandler
from .handlerRelease import ReleaseHandler
from .handlerSingleDish import SingleDishHandler
from .handlerVis import VisHandler
from .phangsLogger import setup_logger

__version__ = version(__name__)

__all__ = [
    "AlmaDownloadHandler",
    "DerivedHandler",
    "ImagingChunkedHandler",
    "ImagingHandler",
    "KeyHandler",
    "PostProcessHandler",
    "ReleaseHandler",
    "SingleDishHandler",
    "VisHandler",
    "setup_logger",
]
