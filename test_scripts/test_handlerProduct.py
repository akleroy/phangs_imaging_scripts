sys.path.insert(1,'/home/saito.50/phangs/phangs_imaging_scripts/')
# Imports
from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerProducts as prh


# Instantiate handlers
this_kh = kh.KeyHandler(master_key = '../key_templates/master_key.txt')
this_prh = prh.ProductHandler(key_handler = this_kh)

# Set which data to process
this_prh.set_interf_configs(only=['7m'])
this_prh.set_feather_configs(only=['7m+tp'])
this_prh.set_line_products(only=['co21'])
this_prh.set_no_cont_products(True)
this_prh.set_targets(only=['ngc4321_1'])

# Run each of the steps individually
this_prh.loop_make_products()

