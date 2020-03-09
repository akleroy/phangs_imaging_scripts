#!/usr/bin/env python
# 
# Run this script inside CASA!
# 

import os, sys
sys.path.append(os.getcwd())
casa_enabled = (sys.argv[0].endswith('start_casa.py'))
if not casa_enabled:
    print('Please run this script inside CASA!')
    sys.exit()

# Set the logging
from phangsPipeline import phangsLogger as pl
reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)

# Imports

#sys.path.insert(1, )
from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerVis as uvh
from phangsPipeline import handlerImaging as imh
from phangsPipeline import handlerPostprocess as pph

# Reloads for debugging
reload(kh)
reload(uvh)
reload(imh)
reload(pph)

# Initialize key handler

this_kh = kh.KeyHandler(master_key = 'phangshi_keys/master_key.txt')
this_uvh = uvh.VisHandler(key_handler = this_kh)
this_imh = imh.ImagingHandler(key_handler = this_kh)
this_pph = pph.PostProcessHandler(key_handler= this_kh)

this_uvh.set_dry_run(False)
#this_uvh.set_targets(skip=['ngc1512','ngc1637','ngc2775','ngc3239','ngc4303','ngc4540','ngc5248'])
#this_uvh.set_targets(only=['ngc1300','ngc5134','ngc3626','ngc3596','ngc2283','ngc4571','ngc5042','ngc5068','ngc4207','ngc3511']
this_uvh.set_targets()

##############################################################################
# Stage the uv data
##############################################################################

# This step loses some archival data due to channel width

this_uvh.loop_stage_uvdata(do_copy=True, do_concat=False, do_remove_staging=False,
                           do_contsub=False, do_extract_line=False, do_extract_cont=False,
                           overwrite=True, timebin='20s')

this_uvh.loop_stage_uvdata(do_copy=False, do_concat=True, do_remove_staging=False,
                           do_contsub=False, do_extract_line=False, do_extract_cont=False,
                           overwrite=True)

this_uvh.loop_stage_uvdata(do_copy=False, do_concat=False, do_remove_staging=True,
                           do_contsub = False, do_extract_line=False, do_extract_cont=False,
                           overwrite=True)

##############################################################################
# Process the staged uv data.
##############################################################################

# Had to exclude some archival data here - NGC 1512

this_uvh.loop_stage_uvdata(do_copy=False, do_concat=False, do_remove_staging=False,
                           do_contsub = True, do_extract_line=False, do_extract_cont=False,
                           overwrite=True)

this_uvh.loop_stage_uvdata(do_copy=False, do_concat=False, do_remove_staging=False,
                           do_contsub = False, do_extract_line=True, do_extract_cont=False,
                           overwrite=True)

##############################################################################
# Step through imaging
##############################################################################

##############################################################################
# Step through postprocessing
##############################################################################

#this_imh.loop_imaging(\
        #do_dirty_image = True, 
        #do_revert_to_dirty = True, 
        #do_read_clean_mask = True, 
        #do_multiscale_clean = True, 
        #do_revert_to_multiscale = True, 
        #do_singlescale_mask = True, 
        #do_singlescale_clean = True, 
        #do_export_to_fits = True, 
        #extra_ext_in = '', 
        #extra_ext_out = '', 
        #)

#this_pph.loop_postprocess(\
        #do_prep=True,
        #do_feather=True, 
        #do_mosaic=True,
        #do_cleanup=True,
        #do_convolve=True,
        #feather_apod=True, 
        #feather_noapod=True,
        #feather_before_mosaic=True,
        #feather_after_mosaic=True,
        #extra_ext_out NOT IMPLEMENTED
        #)

