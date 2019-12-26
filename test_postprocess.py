from phangsPipeline import keyHandler as kh
from phangsPipeline import postprocessHandler as pph

this_kh = kh.KeyHandler()
this_pph = pph.PostProcessHandler(key_handler=this_kh, dry_run=True)

this_pph.set_interf_configs(only=['7m'])
this_pph.set_line_products(only=['co21'])
this_pph.set_no_cont(True)
this_pph.set_targets(only=['ngc3521'])

this_pph.stage_interferometer_data()
