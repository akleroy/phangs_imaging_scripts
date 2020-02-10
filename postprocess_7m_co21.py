# Imports
from phangsPipeline import phangsLogger as pl
from phangsPipeline import keyHandler as kh
from phangsPipeline import postprocessHandler as pph

reload(pph)
reload(kh)
reload(pl)

# Set the logging level
pl.setup_logger(level='DEBUG', logfile=None)

# Instantiate handlers
this_kh = kh.KeyHandler()
this_kh.make_missing_directories(postprocess=True)
this_pph = pph.PostProcessHandler(key_handler=this_kh, dry_run=True)
this_pph.set_dry_run(False)

# Set which data to process
this_pph.set_targets(only=[])
this_pph.set_targets(only=['ngc4321','ngc4321_1','ngc4321_2'])
#this_pph.set_targets(only=['circinus_1','circinus_2'])
this_pph.set_interf_configs(only=['7m'])
this_pph.set_feather_configs(only=['7m+tp'])
this_pph.set_line_products(only=['co21'])
this_pph.set_no_cont(True)

# Run a recipe

this_pph.loop_postprocess(
    do_prep=False,
    do_feather=False, 
    do_mosaic=False,
    do_cleanup=True,
    do_convolve=False,
    feather_apod=True, 
    feather_noapod=True,
    feather_before_mosaic=True,
    feather_after_mosaic=True,
    )

