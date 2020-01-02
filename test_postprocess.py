# Imports
from phangsPipeline import keyHandler as kh
from phangsPipeline import postprocessHandler as pph

# Instantiate handlers
this_kh = kh.KeyHandler()
this_pph = pph.PostProcessHandler(key_handler=this_kh, dry_run=True)
this_pph.set_dry_run(False)

# Set which data to process
this_pph.set_interf_configs(only=['7m'])
this_pph.set_line_products(only=['co21'])
this_pph.set_no_cont(True)
this_pph.set_targets(only=['ngc3521_1'])

# Run each of the steps individually
this_pph.stage_interferometer_data()
this_pph.primary_beam_correct()
this_pph.convolve_to_round_beam()
this_pph.stage_singledish_data()
this_pph.feather()
this_pph.compress()
this_pph.convert()
this_pph.export()

