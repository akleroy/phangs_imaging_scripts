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

this_kh = kh.KeyHandler(master_key = 'key_templates/master_key.txt')
this_uvh = uvh.UVDataHandler(key_handler = this_kh)
this_imh = imh.ImagingHandler(key_handler = this_kh)
this_pph = pph.PostProcessHandler(key_handler= this_kh)

# Set which data to process

#for this_handler in [this_uvh, this_imh, this_pph, this_prh]:
#    this_handler.set_targets(only=['ngc3627'])
#    this_handler.set_interf_configs(only=['7m'])
#    this_handler.set_no_feather_configs(True)
#    #this_handler.set_feather_configs(only=['7m+tp'])
#    this_handler.set_line_products(only=['co21'])
#    this_handler.set_no_cont_products(True)

# Run all steps

this_uvh.loop_stage_uvdata(\
        )
        #do_copy=True, 
        #do_extract_line=True,
        #do_extract_cont=True,
        #do_concat_line=True,
        #do_concat_cont=True,
        #extra_ext = '', 


this_imh.loop_imaging(\
        )
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

this_pph.loop_postprocess(\
        )
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


