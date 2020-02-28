#!/usr/bin/env python
# 

from __future__ import print_function
import os, sys, re, copy, shutil
import numpy as np
import astropy
from astropy.coordinates import SkyCoord, FK5
from astropy.table import Table, Column
import astropy.io.ascii as asciitable
import astroquery
try:
    import gnureadline as readline
except ImportError:
    import readline

if sys.version_info.major <= 2:
    from ConfigParser import ConfigParser
else:
    raw_input = input
    from configparser import ConfigParser


# Get template phangsalma_keys dir, we will copy documentation and header from there. 
template_key_dir = os.path.dirname(os.path.abspath(__file__))+os.sep+'phangsalma_key_templates'+os.sep
if not os.path.isdir(template_key_dir):
    print('Error! Directory "%s" was not found! Please make sure you have downloaded our code completely.'%(template_key_dir))
    sys.exit()



# Check last input records
global confparfile
global confpar
global confparsection
confparfile = 'keyprep.last'
confpar = ConfigParser()
confparsection = None
if os.path.isfile(confparfile):
    print('Reading "%s"'%(confparfile))
    with open(confparfile, 'r') as fp:
        confpar.read_file(fp)
        confparsection = confpar.sections()[0]

#confpar_template_dict = {}
#confpar_template_dict['galaxy_name'] = None
#confpar_template_dict['has_multipart'] = None
#confpar_template_dict['multipart_num'] = None
#confpar_template_dict['phase_center'] = None
#confpar_template_dict['ms_data_path'] = None
#confpar_template_dict['array_tag'] = None
#confpar_template_dict['project_tag'] = None
#confpar_template_dict['more_ms_data'] = None
#confpar_template_dict['galaxy_vsys'] = None
#confpar_template_dict['galaxy_vwidth'] = None
#confpar_template_dict['line_products'] = None
#confpar_template_dict['channel_kms'] = None
#confpar_template_dict['do_continuum'] = None
#confpar_template_dict['key_dir_ok'] = None
#confpar_template_dict['key_dir'] = None
#confpar_template_dict['work_dir_ok'] = None
#confpar_template_dict['work_dir'] = None
    


# Prepare needed variables
galaxy_name_var_dict = {    'name': 'galaxy_name', 
                            'prompt': 'Please input a galaxy name (no whitespace)', 
                            'prompt_fields': None, 
                            'regex': r'^([0-9a-zA-Z_]+)$', 
                            'func': lambda x: x.strip(), 
                            'type': str, 
                            'value': None, 
                            }

has_multipart_var_dict = {  'name': 'has_multipart', 
                            'prompt': 'Does this galaxy have multiple mosaic observations (different phase centers)? [Y/N]', 
                            'prompt_fields': None, 
                            'regex': r'^(Y|N).*$', 
                            'func': lambda x: re.sub(r'^(Y|N).*$', r'\1', x.upper())=='Y', 
                            'type': bool, 
                            'value': None, 
                            }

multipart_num_var_dict = {  'name': 'multipart_num', 
                            'prompt': 'How many different mosaic observations are there for this galaxy? (input an integer)', 
                            'prompt_fields': None, 
                            'regex': r'^([0-9]+)$', 
                            'func': lambda x: int(x), 
                            'type': int, 
                            'value': None, 
                            }

phase_center_var_dict = {   'name': 'phase_center', 
                            'prompt': 'Please input R.A. and Dec. of the phase center of {} (J2000 epoch, FK5 frame)', 
                            'prompt_fields': ['this mosaic observation'], 
                            'regex': r'^([0-9hms:.+-]+) +([0-9dms:.+-]+)$', 
                            'func': lambda x: SkyCoord(re.sub(r'^([0-9hms:.+-]+) +([0-9dms:.+-]+)$', r'\1 \2', x), unit='deg'), 
                            'type': SkyCoord, 
                            'value': None, 
                            }

ms_data_path_var_dict = {   'name': 'ms_data_path', 
                            'prompt': 'Please input the ALMA calibrated Measurement Set paths for {} (no whitespace)', 
                            'prompt_fields': ['this mosaic observation'], 
                            'regex': r'^([^ ]+)$', 
                            'func': lambda x: os.path.abspath(x), 
                            'type': str, 
                            'value': None, 
                            }

array_tag_var_dict = {      'name': 'array_tag', 
                            'prompt': 'Which antenna array is used for {} (e.g., 7m or 12m, no whitespace)?', 
                            'prompt_fields': ['this mosaic observation'], 
                            'regex': r'^([^ ]+)$', 
                            'func': lambda x: x, 
                            'type': str, 
                            'value': None, 
                            }

project_tag_var_dict = {    'name': 'project_tag', 
                            'prompt': 'Which project tag for {} (simple project nick name, no whitespace)?', 
                            'prompt_fields': ['this mosaic observation'], 
                            'regex': r'^([^ ]+)$', 
                            'func': lambda x: x, 
                            'type': str, 
                            'value': None, 
                            }

more_ms_data_var_dict = {   'name': 'more_ms_data', 
                            'prompt': 'More ALMA calibrated Measurement Set to add? [Y/N]', 
                            'prompt_fields': None, 
                            'regex': r'^(Y|N).*$', 
                            'func': lambda x: re.sub(r'^(Y|N).*$', r'\1', x.upper())=='Y', 
                            'type': bool, 
                            'value': None, 
                            }

galaxy_vsys_var_dict = {    'name': 'galaxy_vsys', 
                            'prompt': 'Please input the system velocity of {} (in units of km/s)', 
                            'prompt_fields': ['this galaxy'], 
                            'regex': r'^([0-9eE.]+)$', 
                            'func': lambda x: float(x), 
                            'type': float, 
                            'value': None, 
                            }

galaxy_vwidth_var_dict = {  'name': 'galaxy_vwidth', 
                            'prompt': 'Please input the expected line velocity width of {} (in units of km/s)', 
                            'prompt_fields': ['this galaxy'], 
                            'regex': r'^([0-9eE.]+)$', 
                            'func': lambda x: float(x), 
                            'type': float, 
                            'value': None, 
                            }

line_products_var_dict = {  'name': 'line_products', 
                            'prompt': 'What spectral lines to process (select from co21, 13co, c18o, CI, HI, input the selected ones in one line with whitespace separated)?', 
                            'prompt_fields': None, 
                            'regex': r'^((co21|13co|c18o|CI|HI) *)+$', 
                            'func': lambda x: re.findall(r'\b(co21|13co|c18o|CI|HI)\b$', x, re.IGNORECASE), 
                            'type': list, 
                            'value': None, 
                            } #<TODO># line_list.py

channel_kms_var_dict = {    'name': 'channel_kms', 
                            'prompt': 'Please input the targetting line channel width in velocity (in units of km/s)', 
                            'prompt_fields': None, 
                            'regex': r'^([0-9eE.]+)$', 
                            'func': lambda x: float(x), 
                            'type': float, 
                            'value': None, 
                            }

do_continuum_var_dict = {   'name': 'do_continuum', 
                            'prompt': 'Should we process continuum for {}? [Y/N]', 
                            'prompt_fields': ['this mosaic observation'], 
                            'regex': r'^(Y|N).*$', 
                            'func': lambda x: re.sub(r'^(Y|N).*$', r'\1', x.upper())=='Y', 
                            'type': bool, 
                            'value': None, 
                            }

key_dir_ok_var_dict = {     'name': 'key_dir_ok', 
                            'prompt': 'We will take current directory as the directory to create our pipeline keys under a "phangsalma_keys" sub-directory, is that OK? [Y/N]', 
                            'prompt_fields': None, 
                            'regex': r'^(Y|N).*$', 
                            'func': lambda x: re.sub(r'^(Y|N).*$', r'\1', x.upper())=='Y', 
                            'type': bool, 
                            'value': None, 
                            }

key_dir_path_var_dict = {   'name': 'key_dir', 
                            'prompt': 'Please input your preferred configuration key directory (no whitespace)', 
                            'prompt_fields': None, 
                            'regex': r'^([^ ]+)$', 
                            'func': lambda x: os.path.abspath(x)+os.sep, 
                            'type': str, 
                            'value': None, 
                            }

work_dir_ok_var_dict = {    'name': 'work_dir_ok', 
                            'prompt': 'We will take current directory as the working directory and create our "reduction" sub-directories to work inside it, is that OK? [Y/N]', 
                            'prompt_fields': None, 
                            'regex': r'^(Y|N).*$', 
                            'func': lambda x: re.sub(r'^(Y|N).*$', r'\1', x.upper())=='Y', 
                            'type': bool, 
                            'value': None, 
                            }

work_dir_path_var_dict = {  'name': 'work_dir', 
                            'prompt': 'Please input your preferred working directory (no whitespace)', 
                            'prompt_fields': None, 
                            'regex': r'^([^ ]+)$', 
                            'func': lambda x: os.path.abspath(x)+os.sep, 
                            'type': str, 
                            'value': None, 
                            }



# 
# https://stackoverflow.com/questions/5403138/how-to-set-a-default-editable-string-for-raw-input
# 
def rlinput(prompt, prefill=''):
   readline.set_startup_hook(lambda: readline.insert_text(prefill))
   try:
      return raw_input(prompt)
   finally:
      readline.set_startup_hook()



# Define the function to read user input
def read_user_input(var_dict, repeat = True, check_file_existence = False, check_dir_existence = False):
    # 
    global confparfile
    global confpar
    global confparsection
    # 
    var_value = None
    last_var_str = None
    # 
    if confparsection is not None:
        if var_dict['name'] in confpar[confparsection]:
            last_var_str = confpar[confparsection][var_dict['name']]
    # 
    while True:
        # 
        if var_dict['prompt_fields'] is None:
            var_prompt = var_dict['prompt'] + ' : '
        else:
            var_prompt = var_dict['prompt'].format(*var_dict['prompt_fields']) + ' : '
        # 
        if last_var_str is not None:
            var_str = rlinput(var_prompt, last_var_str)
        else:
            var_str = raw_input(var_prompt)
        # 
        var_ok = False
        if re.match(var_dict['regex'], var_str, re.IGNORECASE):
            var_value = var_dict['func'](var_str)
            print('%s = %r'%(var_dict['name'], var_value))
            var_ok = True
            # 
            if check_dir_existence:
                if not os.path.isdir(var_value):
                    print('Error! Directory "%s" does not exist!'%(var_value))
                    var_ok = False
            # 
            if check_file_existence:
                if not os.path.isfile(var_value):
                    print('Error! File "%s" does not exist!'%(var_value))
                    var_ok = False
        # 
        if var_ok:
            break
        else:
            if repeat:
                print('Sorry, the input "%s" is invalid, please try again.'%(var_str))
            else:
                print('Sorry, the input "%s" is invalid, and we will not repeat our asking. We will just return None.'%(var_str))
                break
    # 
    # update confpar[confparsection]
    if confparsection is not None:
        print('Updating confpar[%r][%r] = %r'%(confparsection, var_dict['name'], var_str))
        confpar[confparsection][var_dict['name']] = var_str
        with open(confparfile, 'w') as confparfp:
            confpar.write(confparfp)
    # 
    return var_value



# Read galaxy name and properties from user command line input
galaxy_name = read_user_input(galaxy_name_var_dict)

# Check whether this galaxy is in our confpar record, so that we can display last inputs
confparsection = galaxy_name
if not (galaxy_name in confpar.sections()):
    confpar[confparsection] = {}

# Read galaxy properties from user command line input
has_multipart = read_user_input(has_multipart_var_dict)
galaxy_multiparts = []
if has_multipart:
    multipart_num = read_user_input(multipart_num_var_dict)
    galaxy_multiparts = ['%s_%d'%(galaxy_name, k) for k in range(multipart_num)]
else:
    multipart_num = 1

# Read galaxy RA Dec and ms data file paths, write to 'target_definitions_list' and 'ms_file_key_list'
target_definitions_list = []
ms_file_key_list = []
for i in range(multipart_num):
    
    # read in galaxy RA Dec
    target_definitions_dict = {}
    
    target_name = galaxy_name if (len(galaxy_multiparts) == 0) else galaxy_multiparts[i]
    phase_center_var_dict_copy = copy.deepcopy(phase_center_var_dict)
    phase_center_var_dict_copy['prompt_fields'][0] = 'the mosaic observation "%s"'%(target_name)
    phase_center = read_user_input(phase_center_var_dict_copy)
    target_definitions_dict['target'] = target_name
    target_definitions_dict['ra'] = phase_center.ra.to_string(sep='hms', precision=3)
    target_definitions_dict['dec'] = phase_center.dec.to_string(sep='dms', precision=3)
    target_definitions_list.append(target_definitions_dict)
    
    # read in ms data file paths
    more_ms_data = True
    ms_data_paths = []
    array_tags = []
    project_tags = []
    while more_ms_data:
        ms_data_path_var_dict_copy = copy.deepcopy(ms_data_path_var_dict)
        ms_data_path_var_dict_copy['prompt_fields'][0] = 'the mosaic observation "%s"'%(target_name)
        ms_data_path = read_user_input(ms_data_path_var_dict_copy, check_dir_existence=True)
        ms_data_paths.append(ms_data_path)
        # array tag
        array_tag = read_user_input(array_tag_var_dict)
        array_tags.append(array_tag)
        # project tag
        project_tag = read_user_input(project_tag_var_dict)
        project_tags.append(project_tag)
        # more to add?
        more_ms_data = read_user_input(more_ms_data_var_dict)
        if not more_ms_data:
            break
    # re-numbering array tags
    array_tag_check_list = list(set(array_tags)) # keep unique arrag_tags
    for array_tag in array_tag_check_list:
        if array_tags.count(array_tag) > 1:
            found_indices = [t for t in range(len(array_tags)) if array_tags[t]==array_tag]
            for j in range(len(found_indices)):
                array_tags[found_indices[t]] = array_tags[found_indices[j]]+'_%d'%(j+1) # update array_tag with number suffix
    
    # add to dict
    for j in range(len(ms_data_paths)):
        ms_file_key_dict = {}
        ms_file_key_dict['target'] = target_name
        ms_file_key_dict['project_tag'] = project_tags[j]
        ms_file_key_dict['array_tag'] = array_tags[j]
        ms_file_key_dict['measurement_set'] = ms_data_paths[j]
        ms_file_key_list.append(ms_file_key_dict)




# Read galaxy vsys and vwidth and add to 'target_definitions_list'
galaxy_vsys = read_user_input(galaxy_vsys_var_dict)
galaxy_vwidth = read_user_input(galaxy_vwidth_var_dict)

for i in range(len(target_definitions_list)):
    target_definitions_list[i]['vsys'] = galaxy_vsys
    target_definitions_list[i]['vwidth'] = galaxy_vwidth



# Read line product
line_products = read_user_input(line_products_var_dict)
channel_kms = read_user_input(channel_kms_var_dict)
do_continuum = read_user_input(do_continuum_var_dict)




# Read key_dir
key_dir_ok = read_user_input(key_dir_ok_var_dict)
if key_dir_ok:
    key_dir_path = os.getcwd()+os.sep+"phangsalma_keys"+os.sep
else:
    key_dir_path = read_user_input(key_dir_path_var_dict)



# Read work_dir
work_dir_ok = read_user_input(work_dir_ok_var_dict)
if work_dir_ok:
    work_dir_path = os.getcwd()+os.sep+"reduction"+os.sep
else:
    work_dir_path = read_user_input(work_dir_path_var_dict)




# Find the common beginning in a list of strings and store as 'ms_root'?
all_ms_files = [t['measurement_set'] for t in ms_file_key_list]
if len(all_ms_files) > 1:
    k = 1
    while np.all([t.startswith(all_ms_files[0][0:k]) for t in all_ms_files]):
        k+=1
else:
    k = len(all_ms_files[0])
ms_root = os.path.dirname(all_ms_files[0][0:k]) # here we do not do k-1 so that there is a trailing letter which is not in common, and we use dirname() to get rid of the tail.
ms_root_strlen = len(ms_root)
for i in range(len(ms_file_key_list)):
    ms_file_key_list[i]['measurement_set'] = ms_file_key_list[i]['measurement_set'][ms_root_strlen:]
    if ms_file_key_list[i]['measurement_set'].startswith(os.sep):
        ms_file_key_list[i]['measurement_set'] = ms_file_key_list[i]['measurement_set'][1:]





# We will create key files in a folder named "phangsalma_keys" under current directory

location_keys = {}
location_keys['key_dir'] = key_dir_path
location_keys['ms_root'] = ms_root
location_keys['cleanmask_root'] = work_dir_path+'cleanmasks'+os.sep
location_keys['singledish_root'] = work_dir_path+'singledish'+os.sep
location_keys['imaging_root'] = work_dir_path+'imaging'+os.sep
location_keys['postprocess_root'] = work_dir_path+'postprocess'+os.sep
location_keys['product_root'] = work_dir_path+'product'+os.sep
location_keys['release_root'] = work_dir_path+'release'+os.sep

data_file_keys = {}
data_file_keys['ms_key'] = 'ms_file_key.txt'
data_file_keys['singledish_key'] = 'singledish_key.txt'
data_file_keys['cleanmask_key'] = 'cleanmask_key.txt'

target_and_config_keys = {}
target_and_config_keys['config_key'] = 'config_definitions.txt'
target_and_config_keys['target_key'] = 'target_definitions.txt'
target_and_config_keys['imaging_key'] = 'imaging_recipes.txt'
target_and_config_keys['linmos_key'] = 'linearmosaic_definitions.txt'
target_and_config_keys['override_key'] = 'overrides.txt'
target_and_config_keys['dir_key'] = 'dir_key.txt'



# Create output directory
if not os.path.isdir(location_keys['key_dir']):
    os.makedirs(location_keys['key_dir'])


# Read/Write 'master_key.txt'
# 
#if not os.path.isfile(location_keys['key_dir']+'master_key.txt'):
#    shutil.copy(template_key_dir+'master_key.txt', location_keys['key_dir']+'master_key.txt')
#    print('Initialized "%s"'%(location_keys['key_dir']+'master_key.txt'))
# 
#master_keys = location_keys + data_file_keys + target_and_config_keys
# 
master_key_file_lines = []
with open(template_key_dir+'master_key.txt', 'r') as fp:
    master_key_file_lines = fp.readlines()
# 
master_key_file_lines_to_write_1 = []
master_key_file_lines_to_write_2 = []
master_key_file_lines_to_write_3 = []
check_point_1 = True
check_point_2 = True
check_point_3 = True
for master_key_file_line in master_key_file_lines:
    if master_key_file_line.startswith('key_dir'):
        check_point_1 = False
    elif master_key_file_line.startswith('ms_key'):
        check_point_2 = False
    elif master_key_file_line.startswith('config_key'):
        check_point_3 = False
    elif master_key_file_line.startswith('#') or master_key_file_line.strip() == '':
        if check_point_1:
            master_key_file_lines_to_write_1.append(master_key_file_line)
        elif check_point_2:
            master_key_file_lines_to_write_2.append(master_key_file_line)
        elif check_point_3:
            master_key_file_lines_to_write_3.append(master_key_file_line)
#keys_to_write = []
#for master_key_file_line in master_key_file_lines:
#    item_to_write = []
#    key_to_write = None
#    item_to_overwrite = []
#    key_to_overwrite = None
#    for key in master_keys:
#        if re.match(r'^(%s)\s+([^\s]+)\s*$'%(key), master_key_file_line):
#            item_read = re.sub(r'^(%s)\s+([^\s]+)\s*$'%(key), r'\2', master_key_file_line).strip()
#            items = master_keys[key]
#            if np.isscalar(items):
#                item = items
#                if not item.endswith(os.sep): 
#                    item += os.sep
#                if item != item_read:
#                    key_to_overwrite = key
#                    item_to_overwrite = [item_read, item] # overwrite this line with the same key but different item
#            else:
#                for item in items:
#                    if not item.endswith(os.sep): 
#                        item += os.sep
#                    if item == item_read:
#                        continue # skip this line with the same key but different item
#                    key_to_write = key
#                    item_to_write.append(item)
#            break # find the key for current line, overwrite it or skip it or add new lines
#    if key_to_write is None and key_to_overwrite is None:
#        master_key_file_lines_to_write.append(master_key_file_line) # keep original line content
#    elif key_to_overwrite is not None: 
#        item_read = item_to_overwrite[0]
#        item = item_to_overwrite[1]
#        key = key_to_overwrite
#        keys_to_write.append(key)
#        print('Overwritting key %s "%s", original "%s"'%(key, item, item_read))
#        master_key_file_lines_to_write.append('%-20s %s\n'%(key, item)) # insert line content
#    elif key_to_write is not None:
#        key = key_to_write
#        keys_to_write.append(key)
#        for item in item_to_write:
#            print('Adding key %s "%s"'%(key, item))
#            master_key_file_lines_to_write.append('%-20s %s\n'%(key, item)) # insert line content
    
if True:
    with open(location_keys['key_dir']+'master_key.txt', 'w') as fp:
        # 
        for master_key_file_line in master_key_file_lines_to_write_1:
            fp.write(master_key_file_line)
        # 
        for key in location_keys:
            items = location_keys[key]
            if np.isscalar(items):
                items = [items]
            for item in items:
                if not item.endswith(os.sep): 
                    item += os.sep
                item_str = '%-20s %s\n'%(key, item)
                fp.write(item_str)
        # 
        for master_key_file_line in master_key_file_lines_to_write_2:
            fp.write(master_key_file_line)
        # 
        #fp.write('\n')
        for key in data_file_keys:
            item = data_file_keys[key]
            item_str = '%-20s %s\n'%(key, item)
            fp.write(item_str)
        # 
        for master_key_file_line in master_key_file_lines_to_write_3:
            fp.write(master_key_file_line)
        # 
        #fp.write('\n')
        for key in target_and_config_keys:
            item = target_and_config_keys[key]
            item_str = '%-20s %s\n'%(key, item)
            fp.write(item_str)
     
    print('Written to "%s"'%(location_keys['key_dir']+'master_key.txt'))



# function to write fixed width table lines into fp
def write_list_of_dict_to_table_file(input_list, fp):
    if len(input_list) == 0:
        return
    colwidths = {} # column width to write
    colnames = list(input_list[0].keys())
    for key in colnames:
        colwidths[key] = max([len(str(t[key])) for t in input_list])
    for i in range(len(input_list)):
        for key in colnames:
            colformat = '{:%d}'%(colwidths[key]+1)
            if key != colnames[0]:
                colformat = '  '+colformat
            fp.write(colformat.format(input_list[i][key])) # write each column of each row
        fp.write('\n')




# Write 'ms_key' and other config files
# 
key_filename = data_file_keys['ms_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    shutil.copy(template_key_dir+os.sep+key_filename, key_filepath)
    print('Initialized "%s"'%(key_filepath))
# 
if len(ms_file_key_list) > 0:
    with open(key_filepath, 'a') as fp:
        fp.write('\n')
        write_list_of_dict_to_table_file(ms_file_key_list, fp)
    # 
    print('Written to "%s"'%(key_filepath))



# 
key_filename = data_file_keys['singledish_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    #with open(key_filepath, 'w') as fp:
    #    fp.write('# column 1: galaxy name\n')
    #    fp.write('# column 2: singledish image fits file (relative to singledish_root)\n')
    #    fp.write('# column 3: product\n')
    #    fp.write('# \n')
    shutil.copy(template_key_dir+os.sep+key_filename, key_filepath)
    print('Initialized "%s"'%(key_filepath))



# 
key_filename = data_file_keys['cleanmask_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    #with open(key_filepath, 'w') as fp:
    #    fp.write('# column 1: galaxy name\n')
    #    fp.write('# column 2: product\n')
    #    fp.write('# column 3: clean mask image fits file (relative to cleanmask_root)\n')
    #    fp.write('# \n')
    #    #<TODO># 
    shutil.copy(template_key_dir+os.sep+key_filename, key_filepath)
    print('Initialized "%s"'%(key_filepath))



# 
key_filename = target_and_config_keys['config_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    #with open(key_filepath, 'w') as fp:
    #    fp.write('# Column 1: type of field\n')
    #    fp.write('# Column 2: name (this can be anything, it appears in file names)\n')
    #    fp.write('# Column 3: parameters (write as literal dictionaries with no spaces)\n')
    #    fp.write('# \n')
    shutil.copy(template_key_dir+os.sep+key_filename, key_filepath)
    print('Initialized "%s"'%(key_filepath))
# 
if True:
    with open(key_filepath, 'a') as fp:
        colwidth1 = 15
        colwidth2 = 15
        colformat = '{:%d} {:%d} {}'%(colwidth1, colwidth2)
        if len(line_products) > 0:
            for line_product in line_products:
                fp.write('\n')
                fp.write(colformat.format("line_product", line_product, "{'line_tag':'%s','channel_kms':%s}"%(line_product, channel_kms)))
                fp.write('\n')
        if do_continuum:
            if len(line_products) > 0:
                lines_to_flag_str = str(line_products).replace(' ','') # remove whitespace as the third column does not allow whitespace
                fp.write('\n')
                fp.write(colformat.format("cont_product", "cont", "{'lines_to_flag':%s}"%(lines_to_flag_str)))
                fp.write('\n')
            else:
                fp.write('\n')
                fp.write(colformat.format("cont_product", "cont", "{'lines_to_flag':[]}"))
                fp.write('\n')
        # 
        array_tags = [t['array_tag'] for t in ms_file_key_list]
        array_tags = list(set(array_tags))
        if '12m' in array_tags:
            fp.write("\n")
            fp.write("interf_config   12m             {'array_tags':['12m']}\n")
            fp.write("interf_config   12m             {'res_min_arcsec':0.5,'res_max_arcsec':7.5,'res_step_factor':1.1}\n")
            fp.write("interf_config   12m             {'clean_scales_arcsec':[0,1,2.5,5.0]}\n")
        if '12m' in array_tags and '7m' in array_tags:
            fp.write("\n")
            fp.write("interf_config   12m+7m          {'array_tags':['12m','7m']}\n")
            fp.write("interf_config   12m+7m          {'res_min_arcsec':0.5,'res_max_arcsec':7.5,'res_step_factor':1.1}\n")
            fp.write("interf_config   12m+7m          {'clean_scales_arcsec':[0,1,2.5,5.0,10.0]}\n")
        if '7m' in array_tags:
            fp.write("\n")
            fp.write("interf_config   7m              {'array_tags':['7m']}\n")
            fp.write("interf_config   7m              {'res_min_arcsec':5.0,'res_max_arcsec':15.0,'res_step_factor':1.1}\n")
            fp.write("interf_config   7m              {'clean_scales_arcsec':[0,5.0,10.0]}\n")
        if '12m' in array_tags and '7m' in array_tags and 'tp' in array_tags:
            fp.write("\n")
            fp.write("feather_config  12m+7m+tp       {'interf_config':'12m+7m'}\n")
            fp.write("feather_config  12m+7m+tp       {'res_min_arcsec':0.5,'res_max_arcsec':7.5,'res_step_factor':1.1}\n")
        if '7m' in array_tags and 'tp' in array_tags:
            fp.write("\n")
            fp.write("feather_config  7m+tp           {'interf_config':'7m'}\n")
            fp.write("feather_config  7m+tp           {'res_min_arcsec':5.0,'res_max_arcsec':15.0,'res_step_factor':1.1}\n")
        # 
    # 
    print('Written to "%s"'%(key_filepath))



# 
# 
key_filename = target_and_config_keys['target_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    #with open(key_filepath, 'w') as fp:
    #    fp.write('# Column 1: target name\n')
    #    fp.write('# Column 2: phase center r.a. string\n')
    #    fp.write('# Column 3: phase center dec string\n')
    #    fp.write('# Column 4: source velocity [km/s]\n')
    #    fp.write('# Column 5: velocity width [km/s]\n')
    shutil.copy(template_key_dir+os.sep+key_filename, key_filepath)
    print('Initialized "%s"'%(key_filepath))
# 
if len(target_definitions_list) > 0:
    with open(key_filepath, 'a') as fp:
        fp.write('\n')
        write_list_of_dict_to_table_file(target_definitions_list, fp)
    # 
    print('Written to "%s"'%(key_filepath))



# 
key_filename = target_and_config_keys['imaging_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    shutil.copy(template_key_dir+os.sep+key_filename, key_filepath)
    print('Initialized "%s"'%(key_filepath))



# 
key_filename = target_and_config_keys['linmos_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    shutil.copy(template_key_dir+os.sep+key_filename, key_filepath)
    print('Initialized "%s"'%(key_filepath))
# 
if has_multipart:
    with open(key_filepath, 'a') as fp:
        fp.write('\n')
        colwidth1 = len(galaxy_name)
        colwidth2 = max([len(t) for t in galaxy_multiparts])
        colwidth1 = colwidth1 if colwidth1>15 else 15
        colwidth2 = colwidth2 if colwidth2>15 else 15
        colformat = '  {:%d}  {:%d}'%(colwidth1+1, colwidth2+1)
        for k in range(len(galaxy_multiparts)):
            fp.write(colformat.format(galaxy_name, galaxy_multiparts[k]))
            fp.write('\n')
    print('Written to "%s"'%(key_filepath))



# 
key_filename = target_and_config_keys['dir_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    shutil.copy(template_key_dir+os.sep+key_filename, key_filepath)
    print('Initialized "%s"'%(key_filepath))
# 
if has_multipart:
    with open(key_filepath, 'a') as fp:
        fp.write('\n')
        colwidth1 = max([len(t) for t in galaxy_multiparts])
        colwidth2 = len(galaxy_name)
        colwidth1 = colwidth1 if colwidth1>15 else 15
        colwidth2 = colwidth2 if colwidth2>15 else 15
        colformat = '  {:%d}  {:%d}'%(colwidth1+1, colwidth2+1)
        for k in range(len(galaxy_multiparts)):
            fp.write(colformat.format(galaxy_multiparts[k], galaxy_name))
            fp.write('\n')
    print('Written to "%s"'%(key_filepath))



# 
key_filename = target_and_config_keys['override_key']
key_filepath = location_keys['key_dir']+key_filename
# 
if not os.path.isfile(key_filepath):
    with open(key_filepath, 'w') as fp:
        fp.write('# column 1: keyword - image name or other file name\n')
        fp.write('# column 2: parameter to override\n')
        fp.write('# column 3: new value\n')
        fp.write('# \n')
        #<TODO># 
    print('Written to "%s"'%(key_filepath))



