#!/usr/bin/env python
# 
# Run this script inside CASA!
# 

import os, sys

# Set the logging
from phangsPipeline import phangsLogger as pl
reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)

# Imports

#sys.path.insert(1, )
from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerProducts as prh
#release handler will go here
reload(kh)
reload(prh)
reload(pl)

# Initialize key handler

this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
this_prh = prh.ProductHandler(key_handler = this_kh)
# release handler

# Set which data to process

#for this_handler in [this_uvh, this_imh, this_pph, this_prh]:
#    this_handler.set_targets(only=['ngc3627'])
#    this_handler.set_interf_configs(only=['7m'])
#    this_handler.set_no_feather_configs(True)
#    #this_handler.set_feather_configs(only=['7m+tp'])
#    this_handler.set_line_products(only=['co21'])
#    this_handler.set_no_cont_products(True)


# Run all steps

this_prh.loop_products_making(\
        )
        # 


#this_rrh.loop_data_release()
# 











