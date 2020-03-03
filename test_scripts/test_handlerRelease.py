sys.path.insert(1,"/home/saito.50/phangs/phangs_imaging_scripts/")
sys.path.insert(1,"/home/saito.50/phangs/phangs_imaging_scripts/phangsPipeline/")
# Imports
from phangsPipeline import handlerKeys as hk
from phangsPipeline import handlerRelease as hr

# Instantiate handlers
this_hk = hk.KeyHandler(master_key = '../phangsalma_keys/master_key.txt')
this_hr = hr.ReleaseHandler(key_handler = this_hk)

# Set which data to process
this_hr.set_interf_configs(only=['7m'])
this_hr.set_line_products(only=['co21'])
this_hr.set_no_cont_products(True)
this_hr.set_targets(only=['ngc4321'])

# Run each of the steps individually
fname_dict = this_hr.loop_build_release()

