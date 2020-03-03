"""
Utilities for reading our pipeline-specific keys.
"""

import os, ast

# There's room to further reduce code and redundancy here. For now,
# this may not be worth the time.

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

##############################################################
# Testing
##############################################################

def test_readers(root=''):
    """
    Test all of the readers by having them read the template
    keys. Hardcoded to the pipeline directory setup.
    """

    test = read_ms_key(root+'key_templates/ms_file_key.txt')
    test = read_singledish_key(root+'key_templates/singledish_key.txt')
    test = read_cleanmask_key(root+'key_templates/cleanmask_key.txt')

    test = read_target_key(root+'key_templates/target_definitions.txt')
    test = read_dir_key(root+'key_templates/dir_key.txt')
    test = read_linmos_key(root+'key_templates/linearmosaic_definitions.txt')

    test = read_config_key(root+'key_templates/config_definitions.txt')

    test = read_override_key(root+'key_templates/overrides.txt')


##############################################################
# Helper Functions
##############################################################

def batch_read(key_list=[], reader_function=None, key_dir='', existing_dict=None):
    """
    Read one set of keys.
    """

    if existing_dict is None:
        output_dict = {}
    else:
        output_dict = existing_dict

    for this_key in key_list:
        this_fname = key_dir + this_key
        if os.path.isfile(this_fname) is False:
            logger.error("I tried to read key "+this_fname+" but it does not exist.")
            continue
        output_dict = reader_function(fname=this_fname, existing_dict=output_dict)

    return(output_dict)

def skip_line(line='', comment='#', delim=None, expected_words=None, expected_format=None):
    if line == '\n':
        return(True)
    if len(line.strip()) == 0:
        return(True)
    if line[0] == comment:
        return(True)
    if expected_words is not None:
        if delim is None:
            words = line.split()
        else:
            words = line.split(delim)
        if len(words) != expected_words:
            logger.warning("Skipping line because it does not match expected format.")
            logger.warning("Expected "+str(expected_words)+" entries. Got "+str(len(words)))
            if expected_format is not None:
                logger.warning("Expected format is: "+expected_format)
            logger.warning("Line is: ")
            logger.warning(line)
            return(True)
    return(False)

def parse_one_line(line='', delim=None):
    """
    Simple helper routine to clean up our readers.
    """
    if delim is None:
        return(line.split())
    else:
        return(line.split(delim))

##############################################################
# UV Data Key
##############################################################
    
# Very specific format

def read_ms_key(fname='', existing_dict=None, delim=None):
    """
    Read a measurement set key.
    """
    
    # Initialize the dictionary

    if existing_dict is None:
        out_dict = {}
    else:
        out_dict = existing_dict

    # Check file existence

    if os.path.isfile(fname) is False:
        logger.error("I tried to read key "+fname+" but it does not exist.")
        return(out_dict)

    logger.info("Reading: "+fname)

    # Expected Format

    expected_words = 6
    expected_format = "target project science_target array_tag obsnum filename"

    # Open File

    infile = open(fname, 'r')
    
    # Loop over the lines    
    lines_read = 0
    while True:
        line  = infile.readline()
        if len(line) == 0:
            break
        if skip_line(line, expected_words=expected_words, delim=delim, expected_format=expected_format):
            continue

        this_target, this_proj, this_field, this_array, this_obsnum, this_file = \
            parse_one_line(line, delim=delim)
        
        # Check if the target is new

        if (this_target in out_dict.keys()) == False:
            out_dict[this_target] = {}

        # Check if the project tag exists already

        if (this_proj not in out_dict[this_target].keys()):
            out_dict[this_target][this_proj] = {}

        # Check if the array tag exists already

        if this_array not in out_dict[this_target][this_proj].keys():
            out_dict[this_target][this_proj][this_array] = {}

        # Check if the obsnum already exists - this may indicate a conflict

        if this_obsnum in out_dict[this_target][this_proj][this_array].keys():
            logger.warning("Possible double entry for. Current line is: ")
            logger.warning(line)
            logger.warning("I will overwrite with the current entry.")

        # Add the science target and filename to the dictionary
        out_dict[this_target][this_proj][this_array][this_obsnum]={}

        out_dict[this_target][this_proj][this_array][this_obsnum]['file'] = this_file
        out_dict[this_target][this_proj][this_array][this_obsnum]['field'] = this_field

        lines_read += 1
            
    infile.close() 

    logger.info("Read "+str(lines_read)+" lines into an ms dictionary.")

    return(out_dict)

##############################################################
# Singledish Data Key and Clean Mask Key
##############################################################

# Right now these are the same. They could split in the future.

def read_targetproductfile_key(fname='', existing_dict=None, delim=None):
    """
    Read a key that maps target and product to file.
    """

    # Initialize the dictionary

    if existing_dict is None:
        out_dict = {}
    else:
        out_dict = existing_dict

    # Check file existence

    if os.path.isfile(fname) is False:
        logger.error("I tried to read key "+fname+" but it does not exist.")
        return(out_dict)
    
    logger.info("Reading: "+fname)
                
    # Expected Format

    expected_words = 3
    expected_format = "target product filename"

    # Open file

    infile = open(fname, 'r')

    # Loop over the lines    
    lines_read = 0
    while True:
        line  = infile.readline()
        if len(line) == 0:
            break
        if skip_line(line, expected_words=expected_words, delim=delim, expected_format=expected_format):
            continue

        # Parse input

        this_target, this_product, this_file = parse_one_line(line, delim=delim)

        if this_target not in out_dict.keys():
            out_dict[this_target] = {}

        if this_product in out_dict[this_target].keys():
            logger.warning("Possible double entry for: "+str(this_target)+' '+str(this_product))
            logger.warning("I will overwrite with the current entry.")

        out_dict[this_target][this_product] = this_file
        lines_read += 1

    # Close and return

    infile.close()

    logger.info("Read "+str(lines_read)+" lines into a target/product/file dictionary.")
    
    return(out_dict)

def read_singledish_key(fname='', existing_dict=None, delim=None):
    """
    Read a single dish key.
    """

    out_dict = read_targetproductfile_key(fname=fname, existing_dict=existing_dict, delim=delim)
    return(out_dict)

def read_cleanmask_key(fname='', existing_dict=None, delim=None):
    """
    Read a clean mask dish key.
    """

    out_dict = read_targetproductfile_key(fname=fname, existing_dict=existing_dict, delim=delim)
    return(out_dict)

##############################################################
# Directory Mapping Key or Linear Mosaic Key
##############################################################

# Right now these are the same, they could split in the future.

def read_dir_key(fname='', existing_dict=None, delim=None):
    """
    Read a directory mapping key.
    """

    out_dict = read_nametoname_key(fname=fname, existing_dict=existing_dict, delim=delim)
    return(out_dict)

def read_linmos_key(fname='', existing_dict=None, delim=None):
    """
    Read a linear mosaic mapping key.
    """

    out_dict = read_nametoname_key(fname=fname, existing_dict=existing_dict, delim=delim,
                                   as_list=True)
    return(out_dict)

def read_nametoname_key(fname='', existing_dict=None, delim=None, as_list=False):
    """
    Read a key that maps one name to another.
    """

    # Initialize the dictionary

    if existing_dict is None:
        out_dict = {}
    else:
        out_dict = existing_dict

    # Check file existence

    if os.path.isfile(fname) is False:
        logger.error("I tried to read key "+fname+" but it does not exist.")
        return(out_dict)

    logger.info("Reading: "+fname)
    
    # Expected Format

    expected_words = 2
    expected_format = "input_name output_name"

    # Open File

    infile = open(fname, 'r')
    
    # Loop over the lines    
    lines_read = 0
    while True:
        line  = infile.readline()
        if len(line) == 0:
            break
        if skip_line(line, expected_words=expected_words, delim=delim, expected_format=expected_format):
            continue

        # Parse input

        this_input, this_output = parse_one_line(line, delim=delim)

        if as_list:
            if this_input not in out_dict.keys():
                out_dict[this_input] = []
            out_dict[this_input].append(this_output)
        else:
            if this_input in out_dict.keys():
                logger.warning("Possible double entry for: "+str(this_input)+' '+str(this_output))
                logger.warning("I will overwrite with the current entry.")                
            out_dict[this_input] = this_output

        lines_read += 1
        
    infile.close()

    logger.info("Read "+str(lines_read)+" lines into a name-to-name mapping dictionary.")

    return(out_dict)

##############################################################
# Target Definition Key
##############################################################

# Very specific format

def read_target_key(fname='', existing_dict=None, delim=None):
    """
    Read a target key.
    """
        
    # Check file existence

    if os.path.isfile(fname) is False:
        logger.error("I tried to read key "+fname+" but it does not exist.")
        return(existing_dict)
    
    logger.info("Reading: "+fname)
    
    # Expected Format

    expected_words = 5
    expected_format = "target raphasectr decphasectr vsys vwidth"

    # Initialize the dictionary

    if existing_dict is None:
        out_dict = {}
    else:
        out_dict = existing_dict

    infile = open(fname, 'r')

    # Loop over the lines    
    lines_read = 0
    while True:
        line  = infile.readline()
        if len(line) == 0:
            break
        if skip_line(line, expected_words=expected_words, delim=delim, expected_format=expected_format):
            continue

        this_target, this_ra, this_dec, this_vsys, this_vwidth = \
            parse_one_line(line, delim=delim)
        
        if (this_target in out_dict.keys()):
            logger.warning("Possible double entry detected. Using newest entry. Line is:")
            logger.warning(line)
 
        out_dict[this_target] = {}
        out_dict[this_target]['rastring'] = this_ra
        out_dict[this_target]['decstring'] = this_dec
        out_dict[this_target]['vsys'] = float(this_vsys)
        out_dict[this_target]['vwidth'] = float(this_vwidth)

        lines_read += 1

    infile.close()

    logger.info("Read "+str(lines_read)+" lines into a target definition dictionary.")
    
    return(out_dict)

##############################################################
# Config Definition Key
##############################################################

# Very specific format

def read_config_key(fname='', existing_dict=None, delim=None):
    """
    Read a configuration key.
    """

    # Check file existence

    if os.path.isfile(fname) is False:
        logger.error("I tried to read key "+fname+" but it does not exist.")
        return(existing_dict)

    logger.info("Reading: "+fname)

    # Expected Format

    expected_words = 3
    expected_format = "config_type config_name params_as_dict"

    # Open File

    infile = open(fname, 'r')

    # Initialize the dictionary

    if existing_dict is None:
        out_dict = {}
    else:
        out_dict = existing_dict

    # Loop over the lines    
    lines_read = 0
    while True:
        line  = infile.readline()
        if len(line) == 0:
            break
        if skip_line(line, expected_words=expected_words, delim=delim, expected_format=expected_format):
            continue

        this_type, this_value, this_params = parse_one_line(line, delim=delim)

        # Check if the type of entry is new
        if (this_type in out_dict.keys()) == False:
            out_dict[this_type] = {}

        # Initialize a configuration on the first entry - configs can have several lines
        if (this_value not in out_dict[this_type].keys()):
            out_dict[this_type][this_value] = {}

        # Parse the parameters as a literal
        try:
            this_params_dict = ast.literal_eval(this_params)
        except:
            logger.error("Could not parse parameters as a dictionary. Line is: ")
            logger.error(line)
            continue

        # Now read in parameters. To do this, define templates for
        # expected fields and data types for each type of
        # configuration. Check to match these.

        if this_type == "interf_config":
            expected_params = {
                'array_tags':[],
                'res_min_arcsec':0.0,
                'res_max_arcsec':0.0,
                'res_step_factor':1.0,
                'clean_scales_arcsec':[]}

        if this_type == "feather_config":
            expected_params = {
                'interf_config':'',
                'res_min_arcsec':0.0,
                'res_max_arcsec':0.0,
                'res_step_factor':1.0}
            
        if this_type == "line_product":
            expected_params = {
                'line_tag':'',
                'channel_kms':0.0}
                
        if this_type == "cont_product":
            expected_params = {
                'lines_to_flag':[]}                

        # Check configs for expected name and data type
                
        for this_key in this_params_dict.keys():
            if this_key not in expected_params.keys():
                logger.error('Got an unexpected parameter key. Line is:')
                logger.error(line)
                continue
            if type(this_params_dict[this_key]) != type(expected_params[this_key]):
                logger.error('Got an unexpected parameter type for parameter '+str(this_key)+'. Line is:')
                logger.error(line)
                continue
            if this_key in out_dict[this_type][this_value].keys():
                logger.debug("Got a repeat parameter definition for "+this_type+" "+this_value)
                logger.debug("Parameter "+this_key+" repeats. Using the latest value.")
                    
            out_dict[this_type][this_value][this_key] = this_params_dict[this_key]

        lines_read += 1

    infile.close()

    logger.info("Read "+str(lines_read)+" lines into a configuration definition dictionary.")

    return(out_dict)

##############################################################
# Override Key
##############################################################

# Format still nebulous

def read_override_key(fname='', existing_dict=None, delim=None):
    """
    Read a file of hand-set overrides.
    """
    
    # Check file existence

    if os.path.isfile(fname) is False:
        logger.error("I tried to read key "+fname+" but it does not exist.")
        return(existing_dict)

    logger.info("Reading: "+fname)

    # Expected Format

    expected_words = 3
    expected_format = "filename param new_value"

    # Open File

    infile = open(fname, 'r')
    
    # Initialize the dictionary

    if existing_dict is None:
        out_dict = {}
    else:
        out_dict = existing_dict

    # Loop over the lines    
    lines_read = 0
    while True:
        line  = infile.readline()
        if len(line) == 0:
            break
        if skip_line(line, expected_words=expected_words, delim=delim, expected_format=expected_format):
            continue

        this_filename, this_param, this_value = parse_one_line(line, delim=delim)

        if this_filename not in out_dict.keys():
            out_dict[this_filename] = {}
        
        out_dict[this_filename][this_param] = this_value
        lines_read += 1

    infile.close()

    logger.info("Read "+str(lines_read)+" lines into an override dictionary.")
    
    return(out_dict)


