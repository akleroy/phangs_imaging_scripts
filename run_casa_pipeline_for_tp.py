"""
To run the Total Power pipeline, make sure you have set up the key files like the following example:

In "keys/master_key.txt":

    singledish_key          singledish_key.txt

In "keys/singledish_key.txt":

    ngc3351      co10      ALMA_TP.NGC3351.CO10.image.VLSRK.fits

In "keys/ms_file_key.txt" (point to the TP observation's member uid directory):

    ngc3351  2022.1.00360  all  tp   1  2022.1.00360.S/science_goal.uid___A001_X2d20_X2d03/group.uid___A001_X2d20_X2d04/member.uid___A001_X2d20_X2d09

In "keys/config_definitions.txt" (no white space inside '{}'):

    singledish_config  tp       {'bl_order':1,'chan_dv_kms':2.5,'doplots':True}

Then, we also need to install the 'astropy' package, because the TP pipeline needs the 'pyfits' package, but casa5?6? does not have that anymore.
If your casa's path is like /software/casa/casa-6.6.4-34-py3.8.el7, then run this to install 'astropy' into the directory './scripts/local/':

    /software/casa/casa-6.6.4-34-py3.8.el7/bin/python3 -m pip install astropy==5.2.2 --prefix=./scripts/local/ --ignore-installed

In this script below, you can see that we will add this path into sys.path so that the pipeline can import the 'astropy' package.

To run the TP pipeline, we will also need the CASA analysisUtil script, which can be obtained from {}. 
We put this analysisUtil script under the directory './scripts/analysis_scripts/'. 
We also put the 'phangsPipeline' under the directory './scripts/phangs_imaging_scripts/'. 
If you have different directories, please change them in this code below. 

"""
import os, sys
sys.path.insert(0, os.getcwd()+'/scripts/analysis_scripts')
sys.path.insert(0, os.getcwd()+'/scripts/phangs_imaging_scripts/phangsPipeline')
sys.path.insert(0, os.getcwd()+'/scripts/phangs_imaging_scripts')
sys.path.insert(0, os.getcwd()+'/scripts/local/lib/python3.8/site-packages')
"""Run these in advance:
    mkdir scripts/local
    /software/casa/casa-6.6.4-34-py3.8.el7/bin/python3 -m pip install astropy==5.2.2 --prefix=./scripts/local/ --ignore-installed
"""
if 'importlib' in sys.modules:
    del sys.modules['importlib']
if 'importlib.metadata' in sys.modules:
    del sys.modules['importlib.metadata']
import packaging as importlib
from importlib import metadata
print('importlib.__path__', importlib.__path__)
#print('metadata.__path__', metadata.__path__)
import astropy

key_file = os.getcwd()+'/keys/master_key.txt'

from phangsPipeline import phangsLogger as pl
pl.setup_logger(level='DEBUG', logfile=None)

from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerSingleDish as sdh

this_kh = kh.KeyHandler(master_key=key_file)
this_sdh = sdh.SingleDishHandler(key_handler=this_kh)

this_sdh.loop_singledish()

