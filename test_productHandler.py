# Imports
from phangsPipeline import keyHandler as kh
from phangsPipeline import productHandler as prh

# Instantiate handlers
this_kh = kh.KeyHandler(master_key = 'phangsalma_keys/master_key.txt')
this_prh = prh.ProductHandler(key_handler = this_kh)

# Set which data to process
this_prh.set_interf_configs(only=['7m'])
this_prh.set_feather_configs(only=['7m+tp'])
this_prh.set_line_products(only=['co21'])
this_prh.set_no_cont_products(True)
this_prh.set_targets(only=['ngc4321'])

# Run each of the steps individually
this_prh.loop_product()

