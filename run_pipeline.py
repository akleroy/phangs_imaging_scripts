# Run this script inside CASA!
# 

import os, sys
casa_enabled = (sys.argv[0].endswith('start_casa.py'))
if not casa_enabled:
    print('Please run this script inside CASA!')
    sys.exit()


# Imports

#sys.path.insert(1, )
from phangsPipeline import keyHandler as kh
from phangsPipeline import uvdataHandler as uvh
from phangsPipeline import imagingHandler as imh
from phangsPipeline import postprocessHandler as pph
from phangsPipeline import productHandler as prh
from phangsPipeline import phangsLogger as pl
reload(kh)
reload(uvh)
reload(imh)
reload(pph)
reload(prh)
reload(pl)


# Set the logging level
pl.setup_logger(level='DEBUG', logfile=None)


# Initialize key handler

this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
this_uvh = uvh.UVDataHandler(key_handler = this_kh)
this_imh = imh.ImagingHandler(key_handler = this_kh)
this_pph = pph.PostProcessHandler(key_handler= this_kh)
this_prh = prh.ProductHandler(key_handler = this_kh)


# Set which data to process

for this_handler in [this_uvh, this_imh, this_pph, this_prh]:
    this_handler.set_targets(only=['ngc3627'])
    this_handler.set_interf_configs(only=['7m'])
    this_handler.set_no_feather_configs(True)
    #this_handler.set_feather_configs(only=['7m+tp'])
    this_handler.set_line_products(only=['co21'])
    this_handler.set_no_cont_products(True)


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


this_prh.loop_products_making(\
        )
        # 


#this_rrh.loop_data_release()
# 











