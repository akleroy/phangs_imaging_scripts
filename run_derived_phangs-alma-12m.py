#!/usr/bin/env python
# 
# Run this script inside CASA!
# 

##############################################################################
# Load routines, initialize handlers
##############################################################################

import os, sys
import importlib

pipedir = '/data/tycho/0/leroy.42/reduction/alma/phangs_imaging_scripts/'
key_file = '/data/tycho/0/leroy.42/reduction/alma/phangs_pipeline_configs/phangs-alma/phangs-alma_keys/master_key.txt'

os.chdir(pipedir)

# Set the logging
from phangsPipeline import phangsLogger as pl
importlib.reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)

# Imports

#sys.path.insert(1, )
from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerDerived as der

# Reloads for debugging
importlib.reload(kh)
importlib.reload(der)

# Initialize key handler

this_kh = kh.KeyHandler(master_key = key_file)
this_der = der.DerivedHandler(key_handler= this_kh)

# Make missing directories

this_kh.make_missing_directories(imaging=True,derived=True,postprocess=True,release=True)

# Set arrays and targets

#this_der.set_targets()
#this_der.set_targets(only=['ngc3627'])
this_der.set_targets(first='ngc1809')
this_der.set_interf_configs(only=['12m+7m'])
this_der.set_feather_configs(only=['12m+7m+tp'])
this_der.set_line_products(only=['co21'])
this_der.set_no_cont_products(True)

# Set steps

do_convolve = False #True
do_noise = False #True
do_strictmask = False #True
do_broadmask = False #True
do_moments = False #True
do_secondary = True #True

##############################################################################
# Step through derived product creation
##############################################################################

if do_convolve:
    this_der.loop_derive_products(do_convolve = True, do_noise = False, 
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = False, do_secondary = False)

if do_noise:
    this_der.loop_derive_products(do_convolve = False, do_noise = True, 
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = False, do_secondary = False)

if do_strictmask:
    this_der.loop_derive_products(do_convolve = False, do_noise = False, 
                                  do_strictmask = True, do_broadmask = False,
                                  do_moments = False, do_secondary = False)

if do_broadmask:
    this_der.loop_derive_products(do_convolve = False, do_noise = False, 
                                  do_strictmask = False, do_broadmask = True,
                                  do_moments = False, do_secondary = False)

if do_moments:
    this_der.loop_derive_products(do_convolve = False, do_noise = False, 
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = True, do_secondary = False)

if do_secondary:
    this_der.loop_derive_products(do_convolve = False, do_noise = False, 
                                  do_strictmask = False, do_broadmask = False,
                                  do_moments = False, do_secondary = True)
