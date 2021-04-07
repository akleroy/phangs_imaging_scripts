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
from phangsPipeline import handlerDerived as dh
#release handler will go here
reload(kh)
reload(dh)
reload(pl)

# Initialize key handler

this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
this_dh = dh.DerivedHandler(key_handler=this_kh)
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

this_dh.loop_derive_products(
    # do_all=False,
    do_convolve=True,
    do_noise=True,
    do_strictmask=True,
    do_broadmask=True,
    do_moments=True,
    do_secondary=True,
    make_directories=True,
    # extra_ext_in='',
    # extra_ext_out='',
    overwrite=True
    )


#this_rrh.loop_data_release()
# 











