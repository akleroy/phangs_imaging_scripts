
# Tested under:
#   CASA 5.4.0

# Imports
from __future__ import print_function
import os, sys
sys.path.append('/home/liu.8691/data_tycho_liu/Scripts/phangs_imaging_scripts')
sys.path.append("/home/maury/leroy.42/casapy/analysis_scripts")
#sys.path.insert(1,"/home/maury/leroy.42/casapy/casa-release-5.4.0-68.el6/lib/python2.7/site-packages")
from phangsPipeline import keyHandler as kh
from phangsPipeline import uvdataHandler as uvh

# Instantiate handlers
reload(kh)
reload(uvh)
this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
this_uvh = uvh.UVDataHandler(key_handler = this_kh, dry_run = False) # dry_run = True

# Set which data to process
this_uvh.set_line_products(only=['co21'])
this_uvh.set_interf_configs(only=['12m+7m'])
#this_uvh.set_no_cont_products(True)
this_uvh.set_targets(only=['ngc3621'])

# Run end-to-end loop
this_uvh.loop_stage_uvdata(\
        do_copy=True, 
        do_extract_line=True,
        do_extract_cont=True,
        do_concat_line=True,
        do_concat_cont=True,
        extra_ext = '', 
        )
        # extra_ext = '_test_2', 

# Run each of the steps individually

