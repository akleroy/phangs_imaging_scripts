from importlib.metadata import PackageNotFoundError, version

# Ensure CASA is installed
from .casa_check import is_casa_installed

casa_enabled = is_casa_installed()

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "dev"

from .handlerAlmaDownload import AlmaDownloadHandler
from .handlerDerived import DerivedHandler
from .handlerKeys import KeyHandler
from .handlerRelease import ReleaseHandler
from .phangsLogger import setup_logger

__all__ = [
    "AlmaDownloadHandler",
    "DerivedHandler",
    "KeyHandler",
    "ReleaseHandler",
    "setup_logger",
]

# Modules that require CASA to be installed
if casa_enabled:
    from .handlerImagingChunked import ImagingChunkedHandler
    from .handlerImaging import ImagingHandler
    from .handlerPostprocess import PostProcessHandler
    from .handlerSingleDish import SingleDishHandler
    from .handlerTestImaging import TestImagingHandler
    from .handlerVis import VisHandler

    __all__.extend([
        "ImagingChunkedHandler",
        "ImagingHandler",
        "PostProcessHandler",
        "SingleDishHandler",
        "TestImagingHandler",
        "VisHandler",
    ])
