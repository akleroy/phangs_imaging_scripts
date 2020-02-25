
# Tested under:
#   CASA 5.4.0

# Imports
from __future__ import print_function
import os, sys
from phangsPipeline import keyHandler as kh
from phangsPipeline import imagingHandler as imh

# Instantiate handlers
reload(kh)
reload(imh)
this_kh = kh.KeyHandler()
this_imh = imh.ImagingHandler(key_handler = this_kh, dry_run = False)

# Set which data to process
this_imh.set_line_products(only=['co21'])
this_imh.set_interf_configs(only=['7m'])
this_imh.set_no_cont_products(True)
this_imh.set_targets(only=['ngc3627'])

# Run end-to-end loop
this_imh.loop_imaging(\
        #make_dirty_image = False, 
        #revert_to_dirty = False, 
        #read_in_clean_mask = False, 
        #run_multiscale_clean = False, 
        #revert_to_multiscale = True, 
        #make_singlescale_mask = True, 
        #run_singlescale_clean = True, 
        do_export_to_fits = True, 
        )
# Error found when using revert_*!

# Run each of the steps individually

