# Imports
from phangsPipeline import keyHandler as kh
from phangsPipeline import imagingHandler as imh

# Instantiate handlers
this_kh = kh.KeyHandler()
this_imh = imh.ImagingHandler(key_handler = this_kh, dry_run = True)

# Set which data to process
this_imh.set_line_products(only=['co21'])
this_imh.set_no_cont(True)
this_imh.set_targets(only=['ngc3627_1'])

# Run end-to-end loop
this_imh.loop_imaging()

# Run each of the steps individually

