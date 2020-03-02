"""
Parts of the PHANGS pipeline that handle the targets, data files,
etc. This is the program that navigates the galaxy list, directory
structure, etc. This part is pure python.
"""

import os, sys, re
import glob
import ast
import numpy as np
from math import floor

try:
    import line_list as ll
except ImportError:
    from phangsPipeline import line_list as ll

try:
    import utilsLists as list_utils
except ImportError:
    from phangsPipeline import utilsLists as list_utils

try:
    import utilsKeyReaders as key_readers
except ImportError:
    from phangsPipeline import utilsKeyReaders as key_readers

try:
    import utilsFilenames as fnames
except ImportError:
    from phangsPipeline import utilsFilenames as fnames

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

VALID_IMAGING_STAGES = ['dirty','multiscale','singlescale']

class KeyHandler:
    """
    Class to handle data files that indicate the names and data sets
    associated with reducing a large ALMA imaging project.
    """

    def __init__(self,
                 master_key = 'key_templates/master_key.txt',
                 quiet=False,
                 dochecks=True,
                 ):

        self._dochecks = dochecks

        self._master_key = None

        self._target_dict = None
        self._config_dict = None
        self._imaging_dict = None
        
        self._ms_dict = None
        self._sd_dict = None
        self._cleanmask_dict = None
        self._linmos_dict = None

        self._dir_for_target = None
        self._override_dict = None

        self.build_key_handler(master_key)

##############################################################
# FILE READING
##############################################################

#region Initialize the key handler and read files.

    def build_key_handler(self,master_key=None):
        """
        Construct the key handler object.
        """

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&&%&%&%&%&%&%&%&%&%&%")
        logger.info("Initializing the PHANGS-ALMA pipeline KeyHandler.")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&&%&%&%&%&%&%&%&%&%&%")
        logger.info("")

        if os.path.isfile(master_key) is False:
            logger.error("Master key "+master_key+" not found. Aborting.")
            raise Exception("Master key "+master_key+" not found. Aborting.")
            return(False)
        pwd = os.getcwd()
        self._master_key = pwd + '/' + master_key
        self._read_master_key()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Reading individual key files.")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")

        self._read_all_keys()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Running checks and cross-links.")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")

        if self._dochecks:
            self.check_ms_existence()

        if self._dochecks:
            self.check_sd_existence()

        self._target_list = []
        self._missing_targets = []
        self._build_target_list()
        if self._dochecks:
            self.print_missing_targets()

        self._expand_dir_key()
        if self._dochecks:
            self.check_dir_existence()

        self._build_whole_target_list()
        self._map_targets_to_mosaics()
        self._map_configs()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Printing configurations.")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")

        self.print_configs()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Printing products.")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")

        self.print_products()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("Master key reading and checks complete.")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&")
        logger.info("")

##############################################################
# READ THE MASTER KEY
##############################################################

# This steers the whole pipeline process and infrastructure and links
# from fields into attributes of the object.

    def _read_master_key(self):
        """
        Read the master key.
        """
        logger.info("Reading the master key.")
        
        logger.info("Master key file: "+self._master_key)
        fname = self._master_key
        infile = open(fname, 'r')
        
        # Initialize

        self._key_dir = None
        self._imaging_root = os.getcwd()+'/../imaging/'
        self._postprocess_root = os.getcwd()+'/../postprocess/'
        self._product_root = os.getcwd()+'/../product/'
        self._release_root = os.getcwd()+'/../release/'

        self._ms_roots = []
        self._sd_roots = []
        self._cleanmask_roots = []
        self._ms_keys = []
        self._sd_keys = []
        self._cleanmask_keys = []

        self._config_keys = []
        self._target_keys = []
        self._linmos_keys = []
        self._dir_keys = []
        self._imaging_keys = []
        self._override_keys = []

        first_key_dir = True
        first_imaging_root = True
        first_postprocess_root = True
        first_product_root = True
        first_release_root = True

        lines_read = 0
        while True:
            line  = infile.readline()
            if len(line) == 0:
                break
            if line[0] == '#' or line == '\n':
                continue

            words = line.split()
            # All key entries go key-value
            if len(words) != 2:
                continue
            
            this_key = words[0]
            this_value = words[1]

            if this_key == 'imaging_root':
                self._imaging_root = this_value
                if first_imaging_root:
                    first_imaging_root = False
                else:
                    logger.warning("Multiple imaging_root definitions. Using the last one.")
                lines_read += 1

            if this_key == 'postprocess_root':
                self._postprocess_root = this_value
                if first_postprocess_root:
                    first_postprocess_root = False
                else:
                    logger.warning("Multiple postprocess_root definitions. Using the last one.")
                lines_read += 1

            if this_key == 'product_root':
                self._product_root = this_value
                if first_product_root:
                    first_product_root = False
                else:
                    logger.warning("Multiple product_root definitions. Using the last one.")
                lines_read += 1

            if this_key == 'release_root':
                self._release_root = this_value
                if first_release_root:
                    first_release_root = False
                else:
                    logger.warning("Multiple release_root definitions. Using the last one.")
                lines_read += 1
            
            if this_key == 'key_dir':
                self._key_dir = this_value
                if first_key_dir:
                    first_key_dir = False
                else:
                    logger.warning("Multiple key directory definitions. Using the last one.")
                lines_read += 1
            
            if this_key == 'ms_root':
                self._ms_roots.append(this_value)
                lines_read += 1

            if this_key == 'singledish_root':
                self._sd_roots.append(this_value)
                lines_read += 1

            if this_key == 'cleanmask_root':
                self._cleanmask_roots.append(this_value)
                lines_read += 1
            
            if this_key == 'ms_key':
                self._ms_keys.append(this_value)
                lines_read += 1

            if this_key == 'config_key':
                self._config_keys.append(this_value)
                lines_read += 1

            if this_key == 'cleanmask_key':
                self._cleanmask_keys.append(this_value)
                lines_read += 1
            
            if this_key == 'dir_key':
                self._dir_keys.append(this_value)
                lines_read += 1
            
            if this_key == 'target_key':
                self._target_keys.append(this_value)
                lines_read += 1

            if this_key == 'imaging_key':
                self._imaging_keys.append(this_value)
                lines_read += 1

            if this_key == 'override_key':
                self._override_keys.append(this_value)
                lines_read += 1
            
            if this_key == 'linmos_key':
                self._linmos_keys.append(this_value)
                lines_read += 1

            if this_key == 'singledish_key':
                self._sd_keys.append(this_value)
                lines_read += 1

        logger.info("Successfully imported "+str(lines_read)+" key/value pairs.")
        
        infile.close()

        if self._dochecks:
            self.check_key_existence()

        return(True)

    def check_key_existence(self):
        """
        Check file existence for the input keys defined in the master file.
        """

        logger.info("------------------------------------")
        logger.info("Checking the existence of key files.")
        logger.info("------------------------------------")
        
        all_valid = True
        errors = 0

        self._key_dir_exists = os.path.isdir(self._key_dir)
        if not self._key_dir_exists:
            logger.error("Missing the key directory. Currently set to "+ self._key_dir)
            logger.error("I need the key directory to proceed. Set key_dir in your master_key file.")
            all_valid = False
            errors += 1
            return(all_valid)

        self._imaging_root_exists = os.path.isdir(self._imaging_root)
        if not self._imaging_root_exists:
            logger.error("The imaging root directory does not exist. Currently set to "+ self._imaging_root)
            logger.error("I need the imaging root to proceed. Set imaging_root in your master_key file.")
            all_valid = False
            errors += 1
            return(all_valid)

        all_key_lists = \
            [self._ms_keys, self._dir_keys, self._target_keys, self._override_keys, self._imaging_keys,
             self._linmos_keys, self._sd_keys, self._config_keys, self._cleanmask_keys]
        for this_list in all_key_lists:
            for this_key in this_list:
                this_key_exists = os.path.isfile(self._key_dir+this_key)
                if not this_key_exists:
                    all_valid = False
                    errors += 1
                    logger.error("key "+ this_key+ " is defined but does not exist in "+self._key_dir)

        if all_valid:
            logger.info("Checked file existence and all files found.")
        else:
            logger.error("Checked file existence. Found "+str(errors)+"errors.")
        
        return(all_valid)

##############################################################
# READ THE INDIVIDUAL KEYS INTO DICTIONARIES
##############################################################

    def _initialize_imaging_dict(self):
        """
        Initialize the imaging dictionary. Requires the config dict to
        exist already to work successfully.
        """

        logger.info("Initializing the imaging dictionary using known configs and stages.")

        full_dict = {}

        # Loop over configs
        for this_config in self.get_interf_configs():
            full_dict[this_config] = {}
        
            # Loop over line products
            for this_product in self.get_line_products():
                full_dict[this_config][this_product] = {}
            
                # Loop over stages
                for this_stage in VALID_IMAGING_STAGES:
                    full_dict[this_config][this_product][this_stage] = []

            # Loop over continuum products
            for this_product in self.get_continuum_products():
                full_dict[this_config][this_product] = {}
                # Loop over stages
                for this_stage in VALID_IMAGING_STAGES:
                    full_dict[this_config][this_product][this_stage] = []

        self._imaging_dict = full_dict

        return()

    def _read_imaging_key(self, fname=None, existing_dict=None, delim=None):
        """
        Read an imaging recipe key. Needs to be in the keyHandler
        object because the order of entries matters AND wild cards
        matter and these interact.
        """
    
        # Check file existence
        
        if os.path.isfile(fname) is False:
            logger.error("I tried to read key "+fname+" but it does not exist.")
            return(existing_dict)

        logger.info("Reading: "+fname)

        # Expected Format

        expected_words = 4
        expected_format = "config product stage recipe"

        # Open File

        infile = open(fname, 'r')
    
        # Initialize the dictionary

        if existing_dict is None:
            if self._imaging_dict is None:
                self._initialize_imaging_dict()
            out_dict = self._imaging_dict
        else:
            out_dict = existing_dict

        # Loop over the lines    
        lines_read = 0
        while True:
            line  = infile.readline()
            if len(line) == 0:
                break
            if key_readers.skip_line(line, expected_words=expected_words, 
                                     delim=delim, expected_format=expected_format):
                continue

            this_config, this_product, this_stage, this_recipe = \
                key_readers.parse_one_line(line, delim=delim)

            if this_config.lower() != 'all':
                if this_config not in self.get_interf_configs():
                    logger.warning("Config not recognized. Line is:")
                    logger.warning(line)
                    continue

            if this_product.lower() != 'all_line' and this_product.lower() != 'all_cont':
                if (this_product not in self.get_line_products()) and \
                        (this_product not in self.get_continuum_products()):
                    logger.warning("Product not recognized. Line is:")
                    logger.warning(line)
                    continue

            if this_stage.lower() != 'all' and this_stage.lower() not in VALID_IMAGING_STAGES:
                logger.warning("Imaging stage not recognized. Line is:")
                logger.warning(line)
                continue
                
            # Just to make the last bit clean, force everything into a list

            if this_config.lower() == 'all':
                config_list = self.get_interf_configs()
            else:
                config_list = [this_config]

            if this_product.lower() == 'all_line':
                product_list = self.get_line_products()
            elif this_product.lower() == 'all_cont':
                product_list = self.get_continuum_products()
            else:
                product_list = [this_product]

            if this_stage == 'all':
                stage_list = VALID_IMAGING_STAGES
            else:
                stage_list = [this_stage]

            for each_config in config_list:
                for each_product in product_list:
                    for each_stage in stage_list:
                        out_dict[each_config][each_product][each_stage] = this_recipe
                        #out_dict[each_config][each_product][each_stage] = this_recipe
                        #logger.debug('self._imaging_dict[%r][%r][%r] = %r'%(each_config, each_product, each_stage, this_recipe))

            lines_read += 1
        
        self._imaging_dict = out_dict

        infile.close()

        logger.info("Read "+str(lines_read)+" lines into the imaging recipe dictionary.")
    
        return out_dict


# Mostly these wrap around the programs in utilsKeyReaders , which
# parse individual key files into dictionaries.        

    def _read_all_keys(self):
        """
        Read the other keys into dictionary attributes.
        """

        self._ms_dict = key_readers.batch_read(
            key_list=self._ms_keys, reader_function=key_readers.read_ms_key,
            key_dir=self._key_dir)

        self._cleanmask_dict = key_readers.batch_read(
            key_list=self._cleanmask_keys, reader_function=key_readers.read_cleanmask_key,
            key_dir=self._key_dir)

        self._sd_dict = key_readers.batch_read(
            key_list=self._sd_keys, reader_function=key_readers.read_singledish_key,
            key_dir=self._key_dir)

        self._config_dict = key_readers.batch_read(
            key_list=self._config_keys, reader_function=key_readers.read_config_key,
            key_dir=self._key_dir)

        self._target_dict = key_readers.batch_read(
            key_list=self._target_keys, reader_function=key_readers.read_target_key,
            key_dir=self._key_dir)

        self._linmos_dict = key_readers.batch_read(
            key_list=self._linmos_keys, reader_function=key_readers.read_linmos_key,
            key_dir=self._key_dir)

        self._dir_for_target = key_readers.batch_read(
            key_list=self._dir_keys, reader_function=key_readers.read_dir_key,
            key_dir=self._key_dir)

        self._initialize_imaging_dict()
        self._imaging_dict = key_readers.batch_read(
            key_list=self._imaging_keys, reader_function=self._read_imaging_key,
            key_dir=self._key_dir, existing_dict=self._imaging_dict)

        self._override_dict = key_readers.batch_read(
            key_list=self._override_keys, reader_function=key_readers.read_override_key,
            key_dir=self._key_dir)

#endregion

##############################################################
# LINKING AND KEY BUILDING
##############################################################

#region Programs to cross-link and build internal structures

    def _build_target_list(self, check=True):
        """
        After the keys have been read in, create a complete target
        list and check for inconsistencies or missing entries.
        """

        logger.info("Building the target list.")
        
        if self._target_dict is None:
            logger.error("I don't have a target dictionary. Can't proceed.")
            return

        self._target_list = list(self._target_dict.keys())
        self._target_list.sort()
        
        self._missing_targets = []

        missing_targets = []

        if self._ms_dict is not None:
            ms_targets = self._ms_dict.keys()
            for target in ms_targets:
                if target not in self._target_list:
                    logger.error(target+ " is in the measurement set key but not the target list.")
                    if target not in missing_targets:
                        missing_targets.append(target)

        if self._dir_for_target is not None:
            dir_targets = self._dir_for_target.keys()
            for target in dir_targets:
                if target not in self._target_list:
                    logger.error(target+ " is in the directory key but not the target list.")
                    if target not in missing_targets:
                        missing_targets.append(target)

        if self._sd_dict is not None:
            sd_targets = self._sd_dict.keys()
            for target in sd_targets:
                if target not in self._target_list:
                    logger.error(target+ " is in the single dish key but not the target list.")
                    if target not in missing_targets:
                        missing_targets.append(target)

        if self._linmos_dict is not None:
            linmos_targets = self._linmos_dict.keys()
            for target in linmos_targets:
                if target not in self._target_list:
                    logger.error(target+ " is in the linear mosaic key but not the target list.")
                    if target not in missing_targets:
                        missing_targets.append(target)
            
        self._missing_targets = missing_targets

        logger.info("Total of "+str(len(self._target_list))+" targets.")
        n_missing = len(self._missing_targets)
        if n_missing == 0:
            logger.info("No cases found where I expect a target but lack a definition.")
        else:
            logger.error(str(n_missing)+" cases where I expected a target definition but didn't find one.")

        return()
        
    def _build_whole_target_list(self):
        """
        Create a list where the parts that will go into a linear
        mosaic are gone and only whole targets remain.
        """

        logger.info("Building the list of whole targets.")

        if self._target_dict is None:
            logger.error("I don't have a target dictionary. Can't proceed.")
            return()

        if self._target_list is None:
            self._build_target_list()

        self._whole_target_list = list(self._target_list)

        if self._linmos_dict is None:
            return()

        linmos_targets = self._linmos_dict.keys()
        for this_mosaic in linmos_targets:
            if this_mosaic not in self._whole_target_list:
                self._whole_target_list.append(this_mosaic)
            for this_part in self._linmos_dict[this_mosaic]:
                if this_part in self._whole_target_list:
                    self._whole_target_list.remove(this_part)

        self._whole_target_list.sort()

        logger.info("Total of "+str(len(self._whole_target_list))+" 'whole' targets.")

        return()

    def _map_targets_to_mosaics(self):
        """
        Create a dictionary that maps from individual targets to mosaics.
        """

        logger.info("Mapping targets to linear mosaics.")

        if self._target_dict is None:
            logger.error("I don't have a target dictionary. Can't proceed.")
            return()

        if self._target_list is None:
            self._build_target_list()

        self._mosaic_assign_dict = {}

        if self._linmos_dict is None:
            return()

        linmos_targets = self._linmos_dict.keys()
        for this_mosaic in linmos_targets:
            part_list = self._linmos_dict[this_mosaic]
            for this_part in part_list:
                self._mosaic_assign_dict[this_part] = this_mosaic

        logger.info("Total of "+str(len(self._mosaic_assign_dict.keys()))+ \
                        " targets assigned to mosaics.")

        return()

    def _map_configs(self):
        """
        Map interferferometric and feather configurations to one another.
        """

        logger.info("Cross-matching interferometric and feather configs")

        if 'interf_config' not in self._config_dict.keys():
            return()

        # Initialize
        for interf_config in self._config_dict['interf_config'].keys():
            
            self._config_dict['interf_config'][interf_config]['feather_config'] = None
            
            if 'feather_config' not in self._config_dict.keys():
                continue
            
            for feather_config in self._config_dict['feather_config'].keys():
                if self._config_dict['feather_config'][feather_config]['interf_config'] != interf_config:
                    continue
                self._config_dict['interf_config'][interf_config]['feather_config'] = feather_config

        return()

#endregion

#region Programs to run checks on the keyHandler and the data

    def print_missing_targets(self):
        """
        Print the targets missing a definition in the target list key.
        """

        if len(self._target_list) == 0:
            self._build_target_list()

        logger.info("-------------------------")
        logger.info("Printing missing targets.")
        logger.info("-------------------------")

        if len(self._missing_targets) == 0:
            logger.info("No missing targets.")
            return(None)
            
        logger.error("Missing a total of "+str(len(self._missing_targets))+" target definitions.")
        for target in self._missing_targets:
            logger.error("missing: "+target)
            
        return()
    
    def check_ms_existence(self):
        """
        Check that the measurement sets in the ms key all exist in the
        specified directories.
        """

        logger.info("-------------------------------------------")
        logger.info("Checking the existence of measurement sets.")
        logger.info("-------------------------------------------")

        if self._ms_dict is None:
            return()

        found_count = 0
        missing_count = 0
        for target in self._ms_dict.keys():
            for project_tag in self._ms_dict[target].keys():
                for array_tag in self._ms_dict[target][project_tag].keys():                    
                    for obs_tag in self._ms_dict[target][project_tag][array_tag].keys():
                        found = False
                        local_found_count = 0
                        for ms_root in self._ms_roots:
                            this_ms = ms_root + self._ms_dict[target][project_tag][array_tag][obs_tag]['file']
                            if os.path.isdir(this_ms):
                                found = True
                                found_count += 1
                                local_found_count += 1
                            if local_found_count > 1:                        
                                logger.error("Found multiple copies of ms for "+target+" "+project_tag+" "+array_tag)
                        if found:
                            continue
                        missing_count += 1
                        logger.error("Missing ms for "+target+" "+project_tag+" "+array_tag)
        
        logger.info("Verified the existence of "+str(found_count)+" measurement sets.")
        if missing_count == 0:
            logger.info("No measurement sets found to be missing.")
        else:
            logger.error("Missing "+str(missing_count)+" measurement set key entries.")

        return()

    def check_sd_existence(self):
        """
        Check that the FITS files in the singledish key all exist in
        the specified directories.
        """

        logger.info("-------------------------------------------")
        logger.info("Checking the existence of single dish data.")
        logger.info("-------------------------------------------")

        if self._sd_dict is None:
            return()

        found_count = 0
        missing_count = 0
        for target in self._sd_dict.keys():
            for product in self._sd_dict[target].keys():
                found = False
                local_found_count = 0
                for this_root in self._sd_roots:
                    this_fname = this_root + self._sd_dict[target][product]
                    if os.path.isfile(this_fname):
                        found = True
                        found_count += 1
                        local_found_count += 1
                if local_found_count > 1:                        
                    logger.error("Found multiple copies of singledish data for "+target+" "+product)
                if found:
                    continue
                missing_count += 1
                logger.error("Missing singledish data for "+target+" "+product)
        
        logger.info("Verified the existence of "+str(found_count)+" single dish data sets.")
        if missing_count == 0:
            logger.info("No single dish data found to be missing.")
        else:
            logger.error("Missing "+str(missing_count)+" single dish key entries.")

        return()

    def _expand_dir_key(self):
        """
        After the keys have been read in, work out the
        directory/target mapping in more detail.
        """
        
        if self._target_list is None:
            self._build_target_list()

        for target in self._target_list:
            if target in self._dir_for_target.keys():
                continue
            self._dir_for_target[target] = target
            
        self._targets_for_dir = {}
        for target in self._target_list:
            this_dir = self._dir_for_target[target]
            if this_dir in self._targets_for_dir.keys():
                current_targets = self._targets_for_dir[this_dir]
                if target in current_targets:
                    continue
                current_targets.append(target)
                current_targets.sort()
                self._targets_for_dir[this_dir] = current_targets
            else:
                self._targets_for_dir[this_dir] = [target]
        
        return()

    def check_dir_existence(self, imaging=True, postprocess=True, product=True):
        """
        Check the existence of the directories for imaging and post-processing.
        """
        
        logger.info("--------------------------------------")
        logger.info("Checking the existence of directories.")
        logger.info("--------------------------------------")

        dir_list = self._targets_for_dir.keys()

        found_dirs=0
        missing_dirs=[]
        for this_dir in dir_list:

            if imaging:
                if os.path.isdir(self._imaging_root+this_dir):
                    found_dirs += 1
                else:
                    logger.warning("Missing imaging directory :"+self._imaging_root+this_dir)
                    missing_dirs.append(self._imaging_root+this_dir)

            if postprocess:
                if os.path.isdir(self._postprocess_root+this_dir):
                    found_dirs += 1
                else:
                    logging.warning("Missing post-processing directory :"+self._postprocess_root+this_dir)
                    missing_dirs.append(self._postprocess_root+this_dir)

            if product:
                if os.path.isdir(self._product_root+this_dir):
                    found_dirs += 1
                else:
                    logging.warning("Missing product directory :"+self._product_root+this_dir)
                    missing_dirs.append(self._product_root+this_dir)
        
        logging.info("Found "+str(found_dirs)+" directories.")

        missing_count = (len(missing_dirs))
        if missing_count == 0:
            logger.info("No directories appear to be missing.")
        else:
            logger.warning("Missing "+str(missing_count)+" directories. Returning that list.")

        return(missing_dirs)
        
#endregion

#region Access data and lists
    
    def _get_dir_for_target(self, target=None, changeto=False, 
                            imaging=False, postprocess=False, product=False):
        """
        Return the imaging working directory given a target name. If
        changeto is true, then change directory to that location.
        """

        if target == None:
            logging.error("Please specify a target.")
            return(None)

        if target not in self._dir_for_target.keys():
            logging.error("Target "+target+" is not in the list of targets.")
            return(None)

        if (imaging and postprocess) or (imaging and product) or \
                (product and postprocess):
            logging.error("Multiple flags set, pick only one type of directory.")
            return(None)
            
        if imaging:
            this_dir = self._imaging_root + self._dir_for_target[target]+'/'
        elif postprocess:
            this_dir = self._postprocess_root + self._dir_for_target[target]+'/'
        elif product:
            this_dir = self._product_root + self._dir_for_target[target]+'/'
        else:
            logging.error("Pick a type of directory. None found.")
            return(None)

        if changeto:
            if not os.path.isdir(this_dir):
                logging.error("You requested to change to "+this_dir+" but it does not exist.")
                return(this_dir)
            os.chdir(this_dir)
            return(this_dir)
        return(this_dir)

    def get_imaging_dir_for_target(self, target=None, changeto=False):
        """
        Return the imaging working directory given a target name. If
        changeto is true, then change directory to that location.
        """
        return(self._get_dir_for_target(target=target, changeto=changeto, imaging=True))

    def get_postprocess_dir_for_target(self, target=None, changeto=False):
        """
        Return the postprocessing working directory given a target name. If
        changeto is true, then change directory to that location.
        """
        return(self._get_dir_for_target(target=target, changeto=changeto, postprocess=True))

    def get_product_dir_for_target(self, target=None, changeto=False):
        """
        Return the product working directory given a target name. If
        changeto is true, then change directory to that location.
        """
        return(self._get_dir_for_target(target=target, changeto=changeto, product=True))

    def get_targets(self, only=None, skip=None, first=None, last=None):
        """
        List the full set of targets.
        
        Modified by keywords only, skip, first, last. Will only return
        targets in only, skip targets in skip, and return targets
        alphabetically after first and before last.
        """
        this_target_list = \
            list_utils.select_from_list(self._target_list, first=first, last=last \
                                       , skip=skip, only=only, loose=True)
        return(this_target_list)

    def get_targets_in_ms_key(self, only=None, skip=None, first=None, last=None):
        """
        List all targets that have uv data.

        Modified by keywords only, skip, first, last. Will only return
        targets in only, skip targets in skip, and return targets
        alphabetically after first and before last.
        """
        targets = self._target_list
        targets_with_ms = []
        for target in targets:
            if target in self._ms_dict.keys():
                targets_with_ms.append(target)

        this_target_list = \
            list_utils.select_from_list(targets_with_ms, first=first, last=last \
                                       , skip=skip, only=only, loose=True)
        return(this_target_list)

    def get_linmos_targets(self, only=None, skip=None, first=None, last=None):
        """
        List all linear mosaics.

        Modified by keywords only, skip, first, last. Will only return
        targets in only, skip targets in skip, and return targets
        alphabetically after first and before last.
        """
        if self._linmos_dict is None:
            return(None)
        mosaic_targets = list(self._linmos_dict.keys())

        this_target_list = \
            list_utils.select_from_list(mosaic_targets, first=first, last=last \
                                       , skip=skip, only=only, loose=True)
        return(this_target_list)


    def get_whole_targets(self, only=None, skip=None, first=None, last=None):
        """
        List only full galaxy names (no parts, e.g., _1 or _2). Very
        similar to the directory list.

        Modified by keywords only, skip, first, last. Will only return
        targets in only, skip targets in skip, and return targets
        alphabetically after first and before last.
        """
        this_whole_target_list = \
            list_utils.select_from_list(self._whole_target_list, first=first, last=last \
                                       , skip=skip, only=only, loose=True)
        return(this_whole_target_list)

    def get_interf_configs(self, only=None, skip=None):
        """
        Get a list of interferometer array configurations

        Modified by keywords only, skip. Will only return targets in
        only, skip targets in skip.
        """
        if type(self._config_dict) != type({}):
            return(None)
        if 'interf_config' not in self._config_dict.keys():
            return(None)
        interf_configs = self._config_dict['interf_config'].keys()
        this_list = \
            list_utils.select_from_list(interf_configs, skip=skip, only=only, loose=True)
        return(this_list)

    def get_feather_configs(self, only=None, skip=None):
        """
        Get a list of feathered single dish + interferometer array
        configruations.

        Modified by keywords only, skip. Will only return targets in
        only, skip targets in skip.
        """
        if type(self._config_dict) != type({}):
            return(None)
        if 'feather_config' not in self._config_dict.keys():
            return(None)
        feather_configs = self._config_dict['feather_config'].keys()
        this_list = \
            list_utils.select_from_list(feather_configs, skip=skip, only=only, loose=True)
        return(this_list)

    def get_line_products(self, only=None, skip=None):
        """
        Get a list of line 'products,' i.e., line plus velocity
        resolution combinations.

        Modified by keywords only, skip. Will only return targets in
        only, skip targets in skip.
        """
        if type(self._config_dict) != type({}):
            return(None)
        if 'line_product' not in self._config_dict.keys():
            return(None)
        line_products = self._config_dict['line_product'].keys()
        this_list = \
            list_utils.select_from_list(line_products, skip=skip, only=only, loose=True)
        return(this_list)

    def get_continuum_products(self, only=None, skip=None):
        """
        Get a list of continuum 'products'.

        Modified by keywords only, skip. Will only return targets in
        only, skip targets in skip.
        """
        if type(self._config_dict) != type({}):
            return(None)
        if 'cont_product' not in self._config_dict.keys():
            return(None)
        cont_products = self._config_dict['cont_product'].keys()
        this_list = \
            list_utils.select_from_list(cont_products, skip=skip, only=only, loose=True)
        return(this_list)

    def is_target_linmos(self, target=None):
        """
        Return True or False based on whether the target is a linear
        mosaic. True means that this target is the OUTPUT of a linear
        mosaic operation.
        """
        if target == None:
            return(False)

        if self._linmos_dict is None:
            logging.error("No linear mosaic dictionary defined.")
            return(False)

        if target in self._linmos_dict.keys():
            return(True)

        return(False)

    def get_parts_for_linmos(self, target=None):
        """
        Return the parts for a linear mosaic.
        """
        if target == None:
            logging.error("Please specify a target.")
            return(None)

        if self._linmos_dict is None:
            logging.error("No linear mosaic dictionary defined.")
            return(None)

        if target not in self._linmos_dict.keys():
            logging.error("No linear mosaic defined for "+target)
            return(None)
        
        return(self._linmos_dict[target])

    def is_target_in_mosaic(self, target):
        """
        Return true or false depending on whether the target is in a linear mosaic.
        """
        return(target in self._mosaic_assign_dict.keys())

    def get_imaging_recipes(self, config=None, product=None, stage=None):
        """
        Return the imaging recipe for the input config and product.
        """
        if config is None or product is None:
            return None
        if config in self._imaging_dict:
            if product in self._imaging_dict[config]:
                imaging_recipe_dir = self._key_dir #<TODO><DL># imaging_recipe_dir
                if not imaging_recipe_dir.endswith(os.sep):
                    imaging_recipe_dir += os.sep
                if stage is not None:
                    if stage in self._imaging_dict[config][product]:
                        return imaging_recipe_dir+self._imaging_dict[config][product][stage] #<TODO><DL># multiple recipes for one stage?
                else:
                    return [imaging_recipe_dir+self._imaging_dict[config][product][t] for t in VALID_IMAGING_STAGES]
        return None

    def set_dochecks(self, dochecks=True):
        """
        Set the feedback level.
        """
        self._dochecks = dochecks
        return

    def get_system_velocity_and_velocity_width_for_target(
        self, 
        target=None,
        ):
        """
        """
        if target is None:
            logging.error("Please specify a target.")
            raise Exception("Please specify a target.")
        
        target_vsys = None
        target_vwidth = None
        if target in self._target_dict:
            if 'vsys' in self._target_dict[target] and 'vwidth' in self._target_dict[target]:
                target_vsys = self._target_dict[target]['vsys']
                target_vwidth = self._target_dict[target]['vwidth']
        
        if target_vsys is None or target_vwidth is None:
            logging.error('No vsys or vwidth values set for the target '+target+'. Please check your "target_definitions.txt".')
            raise Exception('No vsys or vwidth values set for the target '+target)
        
        return(target_vsys, target_vwidth)
    
    def get_phasecenter_for_target(
        self,
        target=None,
        ):
        """
        Return the strings (ra, dec) that define the phase center for
        a target in the target dictionary.
        """
        if target is None:
            logging.error("Please specify a target.")
            raise Exception("Please specify a target.")
        
        rastring = None
        decstring = None
        if target in self._target_dict:
            if 'rastring' in self._target_dict[target]:
                rastring = self._target_dict[target]['rastring']
            if 'decstring' in self._target_dict[target]:
                decstring = self._target_dict[target]['decstring']
        
        if rastring is None or decstring is None:
            logging.error('Missing phase center for target ', target)
            raise Exception('Missing phase center for target ', target)
        
        return(rastring, decstring)

    def get_channel_width_for_line_product(
        self, 
        product=None,
        ):
        """
        """
        if product is None:
            logging.error("Please specify a product.")
            raise Exception("Please specify a product.")
            return None
        
        channel_kms = None
        if 'line_product' in self._config_dict:
            if product in self._config_dict['line_product']:
                if 'channel_kms' in self._config_dict['line_product'][product]:
                    channel_kms = self._config_dict['line_product'][product]['channel_kms']
                #<TODO># maybe we can allow user to set 'channel_width' in units of channels in the "config_definitions.txt" at some point. 
        
        if channel_kms is None:
            logging.error('No channel_kms value set for the input line product '+product+'. Please check your "config_definitions.txt".')
            raise Exception('No channel_kms value set for the input line product '+product)
        
        return channel_kms

    def get_line_tag_for_line_product(
        self, 
        product=None,
        ):
        """
        """
        if product is None:
            logging.error("Please specify a product.")
            raise Exception("Please specify a product.")
            return None
        
        line_tag = None
        if 'line_product' in self._config_dict:
            if product in self._config_dict['line_product']:
                if 'line_tag' in self._config_dict['line_product'][product]:
                    line_tag = self._config_dict['line_product'][product]['line_tag']
        
        if line_tag is None:
            logging.error('No line_tag value set for the input line product '+product+'. Please check your "config_definitions.txt".')
            raise Exception('No line_tag value set for the input line product '+product)
        
        return line_tag

    def get_lines_to_flag_for_continuum_product(
        self, 
        product=None,
        ):
        """
        """
        if product is None:
            logging.error("Please specify a product.")
            raise Exception("Please specify a product.")
            return None
        
        lines_to_flag = []
        if 'cont_product' in self._config_dict:
            if product in self._config_dict['cont_product']:
                if 'lines_to_flag' in self._config_dict['cont_product'][product]:
                    lines_to_flag = self._config_dict['cont_product'][product]['lines_to_flag']
        
        if len(lines_to_flag) == 0:
            logging.error('No lines to flag for the input continuum product '+product)
            #raise Exception('No lines to flag for the input continuum product '+product)
            
        return lines_to_flag
    
    def get_array_tags_for_config(
        self, 
        config=None, 
        ):
        """
        """
        if config is None:
            logging.error("Please specify a config.")
            return None
        
        if 'interf_config' in self._config_dict:
            if config in self._config_dict['interf_config']:
                if 'array_tags' in self._config_dict['interf_config'][config]:
                    return self._config_dict['interf_config'][config]['array_tags']
        
        return None
        
    def loop_over_input_ms(
        self, 
        target=None,
        config=None,
        ):
        """
        Loop over the the target name, project tag, array tag, and
        obsnum for each input visibility file. If a target is supplied
        then restrict to that target, trying to match to a linear
        mosaic if the target is not represented in the dictionary. If
        a configuration is supplied, then only include array_tags that
        contribute to that configuration.
        """
    
        # if user has input a target, match to it
        if target is not None:
            if target in self._ms_dict.keys():
                targets = [target]
            else:
                targets = self.get_parts_for_linmos(target=target)
        
        # if user has input a config, match to it
        if config is not None:
            matching_arraytags = self.get_array_tags_for_config(config)
        else:
            matching_arraytags = []
        
        for this_target in self._ms_dict.keys():

            # if user has input targets, match it
            if target is not None:
                if this_target not in targets:
                    continue

            for this_project in self._ms_dict[this_target].keys():
                for this_arraytag in self._ms_dict[this_target][this_project].keys():

                    # if user has input a config, match it
                    if config is not None:
                        if not (this_arraytag in matching_arraytags):
                            continue
                    for this_obsnum in self._ms_dict[this_target][this_project][this_arraytag].keys():
                        yield this_target, this_project, this_arraytag, this_obsnum

    def get_file_for_input_ms(
        self, 
        target=None, 
        project=None, 
        array_tag=None, 
        obsnum=None,
        ):
        """
        Return the full file path to of an input measurement set given
        a target, project, array_tag, obsnum combination.
        """

        # Verify input

        if target is None:
            logging.error("Please specify a target.")
            return None
        if project is None:
            logging.error("Please specify a project.")
            return None
        if array_tag is None:
            logging.error("Please specify an array_tag.")
            return None
        if obsnum is None:
            logging.error("Please specify an obsnum.")
            return None
        
        # Check that the provided input is in the ms_dict 

        check_valid = False
        if target in self._ms_dict.keys():
            if project in self._ms_dict[target].keys():
                if array_tag in self._ms_dict[target][project].keys():
                    if obsnum in self._ms_dict[target][project][array_tag].keys():
                        check_valid = True

        if not check_valid:
            logger.error('No target '+target+' project '+project+' array_tag '+ \
                             array_tag+' obsnum '+str(obsnum)+' defined in ms dictionary.')
            raise Exception('No target '+target+' project '+project+' array_tag '+ \
                                array_tag+' obsnum '+str(obsnum)+' defined in ms dictionary.')

        # Else locate the full path to the file

        ms_file_path = ''
        file_paths = []
        for ms_root in self._ms_roots:
            file_path = ms_root + self._ms_dict[target][project][array_tag][obsnum]['file']
            if os.path.isdir(file_path):
                file_paths.append(file_path)

        # Error check the results (redundant with read-in checks but good to have)
                
        if len(file_paths) == 0:
            logging.warning('Could not find measurement set for target '+ \
                                target+' project '+project+' array_tag '+ \
                                array_tag+' obsnum '+str(obsnum)+'!')
        else:
            if len(file_paths) > 1:
                logging.warning('Found multiple measurement sets for target '+ \
                                    target+' project '+project+' array_tag '+ \
                                    array_tag+' obsnum '+str(obsnum)+ \
                                    ':\n'+'\n'.join(file_paths)+'\n'+'Returning the first one.')
                
        ms_file_path = file_paths[0]
        return(ms_file_path)
        
    def get_field_for_input_ms(
        self, 
        target=None, 
        project=None, 
        array_tag=None, 
        obsnum=None,
        ):
        """
        Return the science field given a target, project, array_tag,
        obsnum combination.
        """

        # Verify input

        if target is None:
            logging.error("Please specify a target.")
            return None
        if project is None:
            logging.error("Please specify a project.")
            return None
        if array_tag is None:
            logging.error("Please specify an array_tag.")
            return None
        if obsnum is None:
            logging.error("Please specify an obsnum.")
            return None
        
        # Check that the provided input is in the ms_dict 

        check_valid = False
        if target in self._ms_dict.keys():
            if project in self._ms_dict[target].keys():
                if array_tag in self._ms_dict[target][project].keys():
                    if obsnum in self._ms_dict[target][project][array_tag].keys():
                        check_valid = True

        if not check_valid:
            logger.error('No target '+target+' project '+project+' array_tag '+ \
                             array_tag+' obsnum '+str(obsnum)+' defined in ms dictionary.')
            raise Exception('No target '+target+' project '+project+' array_tag '+ \
                                array_tag+' obsnum '+str(obsnum)+' defined in ms dictionary.')

        return(self._ms_dict[target][project][array_tag][obsnum]['field'])

    def has_singledish(
        self,
        target=None, 
        product=None,
        ):
        """
        Return true or false indicating if this target and product
        combination has associated single dish data.
        """

        if self._sd_dict is None:
            logging.error("No single dish key defined.")
            return(None)

        if target == None:
            logging.error("Please specify a target.")
            return(None)

        if product == None:
            logging.error("Please specify a product.")
            return(None)

        if target not in self._sd_dict.keys():
            return(False)
        
        this_dict = self._sd_dict[target]
        
        if product not in this_dict.keys():
            return(False)

        return(True)
        
    def get_sd_filename(
        self, 
        target=None, 
        product=None,
        ):
        """
        Return the single dish filename for a target and product combination.
        """
        
        if self._sd_dict is None:
            logging.error("No single dish key defined.")
            return(None)

        if target == None:
            logging.error("Please specify a target.")
            return(None)

        if product == None:
            logging.error("Please specify a product.")
            return(None)

        if target not in self._sd_dict.keys():
            logging.warning("Not in single dish keys: "+target)
            return(None)

        this_dict = self._sd_dict[target]
        
        if product not in this_dict.keys():
            logging.warning("Product not found for "+target+" : "+product)
            return(None)            
        
        found = False
        found_count = 0
        last_found_file = None
        for this_root in self._sd_roots:
            this_fname = this_root + this_dict[product]
            if os.path.isfile(this_fname) or os.path.isdir(this_fname):
                found = True
                found_count += 1
                last_found_file = this_fname
        
        if found_count > 1:
            logging.error("Found multiple copies of single dish data for "+target+" "+product)
            logging.error("Returning last one, but this is likely an error.")
            return(last_found_file)
        
        if found_count == 0:
            logging.error("Did not find single dish data for "+target+" "+product)
            return(None)
        
        return(last_found_file)

    def get_cleanmask_filename(
        self,
        target=None, 
        product=None,
        ):
        """
        Get the file name of the clean mask associated with a target and product.
        """
        
        if self._cleanmask_dict is None:
            logging.error("No cleanmask dictionary defined.")
            return(None)

        if target == None:
            logging.error("Please specify a target.")
            return(None)

        if target not in self._cleanmask_dict.keys():
            logging.warning("Not in cleanmask keys: "+target)
            return(None)

        this_dict = self._cleanmask_dict[target]

        if 'all' in this_dict.keys():
            this_product = 'all'
        else:
            this_product = product
        
        if this_product not in this_dict.keys():
            logging.warning("Cleanmask not found for "+target+" and product "+str(this_product))
            return(None)            
        
        found = False
        found_count = 0
        last_found_file = None
        for this_root in self._cleanmask_roots:
            this_fname = this_root + this_dict[this_product]
            if os.path.isfile(this_fname) or os.path.isdir(this_fname):
                found = True
                found_count += 1
                last_found_file = this_fname
        
        if found_count > 1:
            logging.error("Found multiple copies of cleanmask for "+target+" "+str(this_product))
            logging.error("Returning last one, but this is likely an error.")
            return(last_found_file)
        
        if found_count == 0:
            logging.error("Did not find a cleanmask for "+target+" "+str(this_product))
            return(None)
        
        logger.debug('Using clean mask file "'+os.path.basename(last_found_file)+'" for target "'+target+'" and product "'+product+'"')
        
        return(last_found_file)

        return ()

    def get_feather_config_for_interf_config(
        self,
        interf_config=None
        ):
        """
        Get the interferometric configuration to go with a feather configuration.
        """
        
        if interf_config is None:
            return(None)

        if 'interf_config' not in self._config_dict.keys():
            return(None)

        interf_config_dict = self._config_dict['interf_config']

        if interf_config not in interf_config_dict.keys():
            return(None)

        if 'feather_config' not in interf_config_dict[interf_config].keys():
            return(None)

        return(interf_config_dict[interf_config]['feather_config'])

    def get_interf_config_for_feather_config(
        self,
        feather_config=None
        ):
        """
        Get the feather configuration to go with an interferometric configuration.
        """
        
        if feather_config is None:
            return(None)

        if 'interf_config' not in self._config_dict.keys():
            return(None)

        feather_config_dict = self._config_dict['feather_config']

        if feather_config not in feather_config_dict.keys():
            return(None)

        if 'interf_config' not in feather_config_dict[feather_config].keys():
            return(None)

        return(feather_config_dict[feather_config]['interf_config'])

    def get_clean_scales_for_config(
        self,
        config=None,
        ):
        """
        Return the angular scales used for multiscale clean for an
        interferometric configuration.
        """

        if config is None:
            return(None)

        if config in self._config_dict['interf_config'].keys():
            this_dict = self._config_dict['interf_config'][config]
        else:
            return(None)
                
        return(this_dict['clean_scales_arcsec'])


    def get_res_for_config(
        self,
        config=None,
        ):
        """
        Return the resolutions used for convolution and product
        creation for an interferometric or feather configuration.
        """
        
        if config is None:
            return(None)

        if config in self._config_dict['interf_config'].keys():
            this_dict = self._config_dict['interf_config'][config]
        elif config in self._config_dict['feather_config'].keys():
            this_dict = self._config_dict['feather_config'][config]
        else:
            return(None)
        
        min_res = this_dict['res_min_arcsec']
        max_res = this_dict['res_max_arcsec']
        step = this_dict['res_step_factor']

        max_steps = np.log10(max_res/min_res)/np.log10(step)+1
        res_array = min_res*step**(np.arange(0.,max_steps,1))
        res_array = res_array[res_array <= max_res]

        return(res_array)

    def get_tag_for_res(
        self,
        res=None,
        ):
        """
        Return a string in the format we use for filenames given a float resolution in arcseconds.
        """

        sigfigs = 2

        res_str = str(int(floor(res)))+'p'+str(int(round((res-floor(res))*10**sigfigs))).zfill(sigfigs)
        
        return(res_str)

    def print_configs(
        self
        ):
        """
        Print out the configurations for inspection.
        """

        if self._config_dict is None:
            return()

        logger.info("Interferometric Configurations")
        for this_config in self._config_dict['interf_config'].keys():
            logger.info("... "+this_config)
            this_arrays = self._config_dict['interf_config'][this_config]['array_tags']
            this_other_config = self._config_dict['interf_config'][this_config]['feather_config']
            this_min_res = self._config_dict['interf_config'][this_config]['res_min_arcsec']
            this_max_res = self._config_dict['interf_config'][this_config]['res_max_arcsec']
            res_step_factor = self._config_dict['interf_config'][this_config]['res_step_factor']
            scales_for_clean = self._config_dict['interf_config'][this_config]['clean_scales_arcsec']
            logger.info("... ... includes arrays "+str(this_arrays))
            logger.info("... ... maps to feather config "+str(this_other_config))
            logger.info("... ... minimum, maximum resolution for products "+str(this_min_res)+' '+str(this_max_res))
            logger.info("... ... step resolution by this factor for products "+str(res_step_factor))
            logger.info("... ... clean these scales in arcsec "+str(scales_for_clean))

        if 'feather_config' in self._config_dict:
            logger.info("Feather Configurations")
            for this_config in self._config_dict['feather_config'].keys():
                logger.info("... "+this_config)
                this_other_config = self._config_dict['feather_config'][this_config]['interf_config']
                this_min_res = self._config_dict['feather_config'][this_config]['res_min_arcsec']
                this_max_res = self._config_dict['feather_config'][this_config]['res_max_arcsec']
                res_step_factor = self._config_dict['feather_config'][this_config]['res_step_factor']
                logger.info("... ... maps to interferometer config "+str(this_other_config))
                logger.info("... ... minimum, maximum resolution for products "+str(this_min_res)+' '+str(this_max_res))
                logger.info("... ... step resolution by this factor for products "+str(res_step_factor))

        return()

    def print_products(
        self
        ):
        """
        Print out the products for inspection.
        """

        if self._config_dict is None:
            return()

        logger.info("Continuum data products")
        for this_product in self._config_dict['cont_product'].keys():
            logger.info("... "+this_product)
            lines_to_flag = self._config_dict['cont_product'][this_product]['lines_to_flag']
            logger.info("... ... lines to flag "+str(lines_to_flag))

        logger.info("Line data products")
        for this_product in self._config_dict['line_product'].keys():
            logger.info("... "+this_product)
            channel_width = self._config_dict['line_product'][this_product]['channel_kms']
            line_name = self._config_dict['line_product'][this_product]['line_tag']
            logger.info("... ... channel width [km/s] "+str(channel_width))
            logger.info("... ... line name code "+str(line_name))

        return()        

    def has_overrides_for_key(self, key=None):
        """Check whether the override dictionary contains the input key or not.
        
        The key is usually a file name or a galaxy name.
        """
        
        if self._override_dict is None:
            return False

        # check key
        if key is None:
            return False
        
        # check key
        if key in self._override_dict:
            return True
        else:
            return False
    
    def get_overrides(self, key=None, param=None, default=None):
        """
        Check the override dictionary for entries given some key,
        parameter pair.
        """
        
        if self._override_dict is None:
            return default

        # check key and param
        if key is None and param is None:
            logger.error('Please input key and param to get overrides!')
            raise Exception('Please input key and param to get overrides!')
            return default
        
        # check key in override dict or not
        if not (key in self._override_dict):
            #logger.debug('Warning! The input key "'+str(key)+'" is not in the overrides file!')
            return default
        
        # if user has input key only, return the dict
        if param is None:
            #logger.debug('Warning! Only has key input when getting overrides! Will return the override dictionary for this key "'+str(key)+'".')
            return self._override_dict[key] #<TODO># dzliu: what if user has just input a key?
        
        if not (param in self._override_dict[key]):
            #logger.debug('Warning! The input key "'+str(key)+'" does not have a param "'+str(param)+'" in the overrides file! Return default value '+str(default))
            return default
        
        logger.debug('Overriding key "'+str(key)+'" param "'+str(param)+'" value '+str(self._override_dict[key][param])+', default ' +str(default) )
        return self._override_dict[key][param]

#endregion
    
#region Manipulate files and file structure

    def make_missing_directories(self, imaging=False, postprocess=False, product=False):
        """
        Make any missing imaging or postprocessing directories.
        """
        
        if not imaging and not postprocess and not product:
            logging.error("Set either imaging or postprocess or product to True. Returning.")
            return(False)

        if imaging:
            if not os.path.isdir(self._imaging_root):
                logging.info("Missing imaging root directory.")
                logging.info("Create: "+self._imaging_root)
                #return(False)
                os.makedirs(self._imaging_root)

        if postprocess:
            if not os.path.isdir(self._postprocess_root):
                logging.info("Missing postprocess root directory.")
                logging.info("Create: "+self._postprocess_root)
                #return(False)
                os.makedirs(self._postprocess_root)

        if product:
            if not os.path.isdir(self._product_root):
                logging.info("Missing product root directory.")
                logging.info("Create: "+self._product_root)
                #return(False)
                os.makedirs(self._product_root)

        missing_dirs = self.check_dir_existence(imaging=imaging, postprocess=postprocess, product=product)
        made_directories = 0
        for this_missing_dir in missing_dirs:
            made_directories += 1
            os.makedirs(this_missing_dir)
        
        missing_dirs = self.check_dir_existence(imaging=imaging, postprocess=postprocess, product=product)

        logging.info("Made "+str(made_directories)+" directories. Now "+str(len(missing_dirs))+" missing.")

        if len(missing_dirs) == 0:
            return(True)

        return(False)

#endregion
