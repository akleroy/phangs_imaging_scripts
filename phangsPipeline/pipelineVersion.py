from __future__ import unicode_literals

# Update this as versions increase
tableversion = '1.6'

try:
    from .version import version
except ImportError:
    # NOTE: this is here to match with previous versions and when
    # the version.py file has not been generated on package installation

    # For now this needs to be updated manually.
    # Setting to last tagged version on 11/11/2021.
    version = "v2.1"
