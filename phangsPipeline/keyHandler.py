"""
Parts of the PHANGS pipeline that handle the targets, data files,
etc. This is the program that navigates the galaxy list, directory
structure, etc. This part is pure python.
"""

import os
import glob
try:
    import line_list as ll
except ImportError:
    from phangsPipeline import line_list as ll

try:
    import utils
except ImportError:
    from phangsPipeline import utils as utils

class KeyHandler:
    """
    Class to handle data files that indicate the names and data sets
    associated with reducing a large ALMA imaging project.
    """

    def __init__(self,
                 master_key = 'phangsalma_keys/master_key.txt',
                 quiet=False,
                 dochecks=True,
                 ):

        self._quiet = quiet
        self._dochecks = dochecks

        self._master_key = None

        self._target_dict = None
        self._config_dict = None
        
        self._ms_dict = None
        self._singledish_dict = None
        self._cleanmask_dict = None
        self._linearmosaic_dict = None

        self._dir_for_target = None
        self._override_dict = None

        self.build_key_handler(master_key)

#region Initialize the key handler and read files.

    def build_key_handler(self,master_key=None):
        """
        Construct the key handler object.
        """

        if not self._quiet:
            print("-------------------------------------------------")
            print("Initializing the PHANGS-ALMA pipeline KeyHandler.")
            print("-------------------------------------------------")

        if not self._quiet:
            print("")
            print("READING THE MASTER KEY.")
            print("")

        if os.path.isfile(master_key) is False:
            print("Master key "+master_key+" not found. Aborting.")
            return False
        pwd = os.getcwd()
        self._master_key = pwd + '/' + master_key
        self._read_master_key()

        if not self._quiet:
            print("")
            print("READING INDIVIDUAL FILES.")
            print("")

        self._read_config_keys()
        self._read_ms_keys()
        self._read_dir_keys()
        self._read_target_keys()    
        self._read_singledish_keys()
        self._read_linearmosaic_keys() 
        self._read_override_keys()
        
        if not self._quiet:
            print("")
            print("RUNNING CHECKS AND CROSS-LINKS.")
            print("")

        if self._dochecks and not self._quiet:
            self.check_ms_existence()

        if self._dochecks and not self._quiet:
            self.check_singledish_existence()

        self._target_list = []
        self._missing_targets = []
        self._build_target_list()
        if self._dochecks and not self._quiet:
            self.check_target_existence()

        self._expand_dir_key()
        if self._dochecks and not self._quiet:
            self.check_dir_existence()

        self._build_whole_target_list()
        self._map_targets_to_mosaics()

        if not self._quiet:
            print("")
            print("COMPLETE.")
            print("")

    def _read_master_key(self):
        """
        Read the master key.
        """
        
        fname = self._master_key
        infile = open(fname, 'r')
        
        # Initialize

        self._key_dir = None
        self._imaging_root = os.getcwd()+'/../'
        self._postprocess_root = None
        self._release_root = None

        self._ms_roots = []
        self._singledish_roots = []
        self._ms_keys = []
        self._singledish_keys = []
        self._cleanmask_keys = []

        self._config_keys = []
        self._target_keys = []
        self._linearmosaic_keys = []
        self._dir_keys = []
        self._override_keys = []

        first_imaging_root = True
        first_key_dir = True
        first_postprocess_root = True
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
                    print("Warning, multiple imaging_root definitions. Using the last one.")
                lines_read += 1

            if this_key == 'postprocess_root':
                self._postprocess_root = this_value
                if first_postprocess_root:
                    first_postprocess_root = False
                else:
                    print("Warning, multiple postprocess_root definitions. Using the last one.")
                lines_read += 1

            if this_key == 'release_root':
                self._release_root = this_value
                if first_release_root:
                    first_release_root = False
                else:
                    print("Warning, multiple release_root definitions. Using the last one.")
                lines_read += 1
            
            if this_key == 'key_dir':
                self._key_dir = this_value
                if first_key_dir:
                    first_key_dir = False
                else:
                    print("Warning, multiple key directory definitions. Using the last one.")
                lines_read += 1
            
            if this_key == 'ms_root':
                self._ms_roots.append(this_value)
                lines_read += 1

            if this_key == 'singledish_root':
                self._singledish_roots.append(this_value)
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

            if this_key == 'override_key':
                self._override_keys.append(this_value)
                lines_read += 1
            
            if this_key == 'linearmosaic_key':
                self._linearmosaic_keys.append(this_value)
                lines_read += 1

            if this_key == 'singledish_key':
                self._singledish_keys.append(this_value)
                lines_read += 1

        if not self._quiet:
            print("Read in master key: "+self._master_key)
            print("... successfully imported "+str(lines_read)+" key/value pairs.")
        
        infile.close()

        if self._dochecks:
            self.check_key_existence()

    def check_key_existence(self):
        """
        Check file existence for the input keys defined in the master file.
        """

        all_valid = True
        errors = 0

        self._key_dir_exists = os.path.isdir(self._key_dir)
        if not self._key_dir_exists:
            if not self._quiet:
                print("Missing the key directory. Currently set to ", self._key_dir)
                print("I need this to proceed. Set key_dir in your master_key file.")
            all_valid = False
            errors += 1
            return(all_valid)

        self._imaging_root_exists = os.path.isdir(self._imaging_root)
        if not self._imaging_root_exists:
            if not self._quiet:
                print("The imaging root directory does not exist. Currently set to ", self._imaging_root)
                print("I need this to proceed. Set imaging_root in your master_key file.")
            all_valid = False
            errors += 1
            return(all_valid)

        all_key_lists = \
            [self._ms_keys, self._dir_keys, self._target_keys, self._override_keys,
             self._linearmosaic_keys, self._singledish_keys, self._config_keys, self._cleanmask_keys]
        for this_list in all_key_lists:
            for this_key in this_list:
                this_key_exists = os.path.isfile(self._key_dir+this_key)
                if not this_key_exists:
                    all_valid = False
                    errors += 1
                    if not self._quiet:
                        print("key ", this_key, " is defined but does not exist in "+self._key_dir)

        if not self._quiet:
            if all_valid:
                print("Checked file existence and all files found.")
            else:
                print("Checked file existence. Found "+str(errors)+"errors.")
        
        return(all_valid)

    def _read_ms_keys(self):
        """
        Read all of the measurement set keys.
        """

        self._ms_dict = {}
        for this_key in self._ms_keys:
            this_fname = self._key_dir + this_key
            if os.path.isfile(this_fname) is False:
                if not self._quiet:
                    print("I tried to read key "+this_fname+" but it does not exist.")
                return()
            self._read_one_ms_key(fname=this_fname)

        return()

    def _read_one_ms_key(self, fname=None):
        """
        Read one measurement set key. Called by _read_ms_keys. 
        """

        # should not get to this, but check just in case
        if os.path.isfile(fname) is False:
            if not self._quiet:
                print("I tried to read key "+fname+" but it does not exist.")
            return()

        infile = open(fname, 'r')
        
        lines_read = 0
        while True:
            line  = infile.readline()
            if len(line) == 0:
                break
            if line[0] == '#' or line == '\n':
                continue

            words = line.split()

            if len(words) < 4:
                if not self._quiet:
                    print("Skipping line because it does not match ms_key format.")
                    print("Format is space delimited: target project_tag array_tag filename .")
                    print("Line is: ")
                    print(line)
                continue
            if len(words) > 4:
                if not self._quiet:
                    print("Skipping line because it does not match ms_key format.")
                    print("Format is space delimited: target project_tag array_tag filename .")
                    print("Line is: ")
                    print(line)
                continue

            this_target = words[0]
            this_proj_tag = words[1]
            this_array_tag = words[2]
            this_file = words[3]

            # Initialize the dictionary the first time we get a result.
            if self._ms_dict is None:
                self._ms_dict = {}

            # Check if the target is new
            if (this_target in self._ms_dict.keys()) == False:
                self._ms_dict[this_target] = {}

            # Check if the project tag exists already
            if (this_proj_tag not in self._ms_dict[this_target].keys()):
                self._ms_dict[this_target][this_proj_tag] = {}

            # Check if the array tag exists already
            if this_array_tag in self._ms_dict[this_target][this_proj_tag].keys():
                print("Possible double entry/mistake for: ", this_target, this_proj_tag, this_array_tag)

            # Add the filename to the dictionary
            self._ms_dict[this_target][this_proj_tag][this_array_tag] = this_file
            lines_read += 1
            
        infile.close() 

        if not self._quiet:
            print("Read "+str(lines_read)+" lines into the ms dictionary.")

        return()

    def _read_dir_keys(self):
        """
        Read all of the directory keys.
        """

        self._dir_for_target = {}
        for this_key in self._dir_keys:
            this_fname = self._key_dir + this_key
            if os.path.isfile(this_fname) is False:
                if not self._quiet:
                    print("I tried to read key "+this_fname+" but it does not exist.")
                return()
            self._read_one_dir_key(fname=this_fname)

        return()

    def _read_one_dir_key(self, fname=None):
        """
        Read one directory key. This is called by read_dir_keys.
        """

        # should not get to this, but check just in case
        if os.path.isfile(fname) is False:
            if not self._quiet:
                print("I tried to read key "+fname+" but it does not exist.")
            return()

        infile = open(fname, 'r')
    
        lines_read = 0

        while True:

            line  = infile.readline()    
            if len(line) == 0:
                break
            if line[0] == '#' or line == '\n':
                continue

            words = line.split()
                
            if len(words) < 2:
                if not self._quiet:
                    print("Skipping line because it does not match dir_key format.")
                    print("Format is space delimited: target directory .")
                    print("Line is: ")
                    print(line)
                continue
            if len(words) > 2:
                if not self._quiet:
                    print("Skipping line because it does not match dir_key format.")
                    print("Format is space delimited: target directory .")
                    print("Line is: ")
                    print(line)
                continue

            this_target = words[0]
            this_dir = words[1]    

            # Initialize the dictionary the first time we get a result.
            if self._dir_for_target is None:
                self._dir_for_target = {}

            self._dir_for_target[this_target] = this_dir
            lines_read += 1
        
        infile.close()

        if not self._quiet:
            print("Read "+str(lines_read)+" lines into the directory mapping dictionary.")

        return()

    def _read_target_keys(self):
        """
        Read all of the target keys.
        """
        
        self._target_dict = {}

        for this_key in self._target_keys:
            this_fname = self._key_dir + this_key
            if os.path.isfile(this_fname) is False:
                if not self._quiet:
                    print("I tried to read key "+fname+" but it does not exist.")
                return()
            self._read_one_target_key(fname=this_fname)

        return()

    def _read_one_target_key(self, fname=None):
        """
        Read one target key. Called by read_target_keys.
        """
        
        # should not get to this, but check just in case
        if os.path.isfile(fname) is False:
            self.update_files()
            return()

        infile = open(fname, 'r')

        lines_read = 0
        while True:
            line  = infile.readline()    
            if len(line) == 0:
                break
            if line[0] == '#' or line == '\n':
                continue
            words = line.split()

            if len(words) < 5:
                if not self._quiet:
                    print("Skipping line because it does not match target_key format.")
                    print("Format is space delimited: target RA Dec vsys vwidth .")
                    print("Line is: ")
                    print(line)
                continue
            if len(words) > 5:
                if not self._quiet:
                    print("Skipping line because it does not match target_key format.")
                    print("Format is space delimited: target RA Dec vsys vwidth .")
                    print("Line is: ")
                    print(line)
                continue

            # Initialize the dictionary the first time we get a result.
            if self._target_dict is None:
                self._target_dict = {}

            this_target = words[0]
            this_ra = words[1]
            this_dec = words[2]
            this_vsys = words[3]
            this_vwidth = words[4]
        
            self._target_dict[this_target] = {}
            self._target_dict[this_target]['rastring'] = this_ra
            self._target_dict[this_target]['decstring'] = this_dec
            self._target_dict[this_target]['vsys'] = float(this_vsys)
            self._target_dict[this_target]['vwidth'] = float(this_vwidth)

            lines_read += 1

        infile.close()

        if not self._quiet:
            print("Read "+str(lines_read)+" lines into the target definition dictionary.")
    
        return()

    def _read_singledish_keys(self):
        """
        Read all of the single dish keys.
        """

        self._singledish_dict = {}
        for this_key in self._singledish_keys:
            this_fname = self._key_dir + this_key
            if os.path.isfile(this_fname) is False:
                if not self._quiet:
                    print("I tried to read key "+this_fname+" but it does not exist.")
                return()
            self._read_one_singledish_key(fname=this_fname)

        return()

    def _read_one_singledish_key(self,fname=None):
        """
        Read one single dish key. Called by read_single_dish_keys()
        """
                
        # should not get to this, but check just in case
        if os.path.isfile(fname) is False:
            self.update_files()
            return()

        infile = open(fname, 'r')

        lines_read = 0
        while True:

            line = infile.readline()    

            if len(line) == 0:
                break

            if line[0] == '#' or line == '\n':
                continue

            words = line.split()

            if len(words) < 3:
                if not self._quiet:
                    print("Skipping line because it does not match singledish_key format.")
                    print("Format is space delimited: target filename product/line .")
                    print("Line is: ")
                    print(line)
                continue
            if len(words) > 3:
                if not self._quiet:
                    print("Skipping line because it does not match singledish_key format.")
                    print("Format is space delimited: target filename product/line .")
                    print("Line is: ")
                    print(line)
                continue

            this_target = words[0]
            this_file = words[1]
            this_product = words[2]

            # Initialize the dictionary the first time we get a result.
            if self._singledish_dict is None:
                self._singledish_dict = {}
        
            if this_target not in self._singledish_dict.keys():
                self._singledish_dict[this_target] = {}

            self._singledish_dict[this_target][this_product] = this_file
            lines_read += 1

        infile.close()

        if not self._quiet:
            print("Read "+str(lines_read)+" lines into the single dish dictionary.")
    
        return()

    def _read_config_keys(self):
        """
        Read all of the configuration definitions.
        """

        self._config_dict = {}
        for this_key in self._config_keys:
            this_fname = self._key_dir + this_key
            if os.path.isfile(this_fname) is False:
                if not self._quiet:
                    print("I tried to read key "+this_fname+" but it does not exist.")
                return()
            self._read_one_config_key(fname=this_fname)

        return()

    def _read_one_config_key(self, fname=None):
        """
        Read one configuration key. Called by _read_config_keys().
        """
        
        # should not get to this, but check just in case
        if os.path.isfile(fname) is False:
            self.update_files()
            return()

        infile = open(fname, 'r')

        lines_read = 0
        while True:
            line = infile.readline()    
            if len(line) == 0:
                break
            if line[0] == '#' or line == '\n':
                continue

            words = line.split()
            if len(words) < 3:
                if not self._quiet:
                    print("Skipping line because it does not match configuration definition format.")
                    print("Format is space delimited: vis_file param new_value .")
                    print("Line is: ")
                    print(line)
                continue
            if len(words) > 3:
                if not self._quiet:
                    print("Skipping line because it does not match configuration definition format.")
                    print("Format is space delimited: vis_file param new_value .")
                    print("Line is: ")
                    print(line)
                continue

            this_type = words[0]
            this_value = words[1]
            this_params = words[2]

            # Initialize the dictionary the first time we get a result.
            if self._config_dict is None:
                self._config_dict = {}

            # Check if the type of entry is new
            if (this_type in self._config_dict.keys()) == False:
                self._config_dict[this_type] = {}

            # Check if the particular entry already exists
            if (this_value not in self._config_dict[this_type].keys()):
                self._config_dict[this_type][this_value] = {}
            else:
                print("Got a repeat configuration definition of "+this_type+" "+this_value)
                print("... using the latest one, but this may be a problem.")

            # Now do logic on particular types of setups.

            if this_type == "interf_config":
                arrays_to_include = this_params.split(',')
                self._config_dict[this_type][this_value]['arrays'] = \
                    arrays_to_include

            if this_type == "feather_config":
                arrays_to_include = this_params
                self._config_dict[this_type][this_value]['interf_configs'] = \
                    arrays_to_include

            if this_type == "line_product":
                defn = this_params.split(',')
                self._config_dict[this_type][this_value]['line'] = \
                    defn[0]
                channel = float((defn[1].split('kms'))[0])
                self._config_dict[this_type][this_value]['channel'] = \
                    channel
                
            if this_type == "cont_product":
                lines_to_flag = this_params.split(',')
                self._config_dict[this_type][this_value]['lines_to_flag'] = \
                    lines_to_flag

            lines_read += 1

        infile.close()

        if not self._quiet:
            print("Read "+str(lines_read)+" lines into the configuration definition dictionary.")

        return()

    def _read_linearmosaic_keys(self):
        """
        Read all of the linear mosaic keys.
        """

        self._linearmosaic_dict = {}
        for this_key in self._linearmosaic_keys:
            this_fname = self._key_dir + this_key
            if os.path.isfile(this_fname) is False:
                if not self._quiet:
                    print("I tried to read key "+this_fname+" but it does not exist.")
                return()
            self._read_one_linearmosaic_key(fname=this_fname)

        return()

    def _read_one_linearmosaic_key(self, fname=None):
        """
        Read one file containing the assignments of targets to linear
        mosaics. Called by read_linearmosaic_keys().
        """
                
        # should not get to this, but check just in case
        if os.path.isfile(fname) is False:
            self.update_files()
            return()

        infile = open(fname, 'r')
            
        lines_read = 0
        while True:
            line  = infile.readline()    
            if len(line) == 0:
                break
            if line[0] == '#' or line == '\n':
                continue
            words = line.split()

            if len(words) < 2:
                if not self._quiet:
                    print("Skipping line because it does not match linearmosaic key format.")
                    print("Format is space delimited: mosaic_target assigned_target .")
                    print("Line is: ")
                    print(line)
                continue
            if len(words) > 2:
                if not self._quiet:
                    print("Skipping line because it does not match linearmosaic key format.")
                    print("Format is space delimited: mosaic_target assigned_target .")
                    print("Line is: ")
                    print(line)
                continue

            this_mosaic_name = words[0]
            this_target_name = words[1]
    
            # Initialize the dictionary the first time we get a result.
            if self._linearmosaic_dict is None:
                self._linearmosaic_dict = {}

            if this_mosaic_name not in self._linearmosaic_dict.keys():
                self._linearmosaic_dict[this_mosaic_name] = [this_target_name]
            else:
                current_targets = self._linearmosaic_dict[this_mosaic_name]
                if this_target_name not in current_targets:
                    current_targets.append(this_target_name)
                self._linearmosaic_dict[this_mosaic_name] = current_targets

            lines_read += 1

        infile.close()

        if not self._quiet:
            print("Read "+str(lines_read)+" lines into the linear mosaic definition dictionary.")

        return()

    def _read_override_keys(self):
        """
        Read all of the override keys.
        """

        self._override_dict = {}
        for this_key in self._override_keys:
            this_fname = self._key_dir + this_key
            if os.path.isfile(this_fname) is False:
                if not self._quiet:
                    print("I tried to read key "+this_fname+" but it does not exist.")
                return()
            self._read_one_override_key(fname=this_fname)

        return()

    def _read_one_override_key(self, fname=None):
        """
        Read one file of hand-set overrides. Called by read_override_keys().
        """
        
        # should not get to this, but check just in case
        if os.path.isfile(fname) is False:
            self.update_files()
            return()

        infile = open(fname, 'r')

        lines_read = 0
        while True:
            line = infile.readline()    
            if len(line) == 0:
                break
            if line[0] == '#' or line == '\n':
                continue

            words = line.split()
            if len(words) < 3:
                if not self._quiet:
                    print("Skipping line because it does not match override format.")
                    print("Format is space delimited: vis_file param new_value .")
                    print("Line is: ")
                    print(line)
                continue
            if len(words) > 3:
                if not self._quiet:
                    print("Skipping line because it does not match override format.")
                    print("Format is space delimited: vis_file param new_value .")
                    print("Line is: ")
                    print(line)
                continue

            this_target = words[0]
            this_param = words[1]
            this_value = words[2]

            # Initialize the dictionary the first time we get a result
            if self._override_dict is None:
                self._override_dict = {}

            if this_target not in self._override_dict.keys():
                self._override_dict[this_target] = {}
            self._override_dict[this_target][this_param] = this_value
            lines_read += 1

        infile.close()

        if not self._quiet:
            print("Read "+str(lines_read)+" lines into the override dictionary.")
    
        return()

#endregion

#region Programs to cross-link and build internal structures

    def _build_target_list(self, check=True):
        """
        After the keys have been read in, create a complete target
        list and check for inconsistencies or missing entries.
        """
        
        if self._target_dict is None:
            print("I don't have a target dictionary. Can't proceed.")
            return

        self._target_list = list(self._target_dict.keys())
        self._target_list.sort()
        
        self._missing_targets = []

        missing_targets = []

        if self._ms_dict is not None:
            ms_targets = self._ms_dict.keys()
            for target in ms_targets:
                if target not in self._target_list:
                    if not self._quiet:
                        print(target, " is in the measurement set key but not the target list.")
                    if target not in missing_targets:
                        missing_targets.append(target)

        if self._dir_for_target is not None:
            dir_targets = self._dir_for_target.keys()
            for target in dir_targets:
                if target not in self._target_list:
                    if not self._quiet:
                        print(target, " is in the directory key but not the target list.")
                    if target not in missing_targets:
                        missing_targets.append(target)

        if self._singledish_dict is not None:
            singledish_targets = self._singledish_dict.keys()
            for target in singledish_targets:
                if target not in self._target_list:
                    if not self._quiet:
                        print(target, " is in the single dish key but not the target list.")
                    if target not in missing_targets:
                        missing_targets.append(target)

        if self._linearmosaic_dict is not None:
            linearmosaic_targets = self._linearmosaic_dict.keys()
            for target in linearmosaic_targets:
                if target not in self._target_list:
                    if not self._quiet:
                        print(target, " is in the linear mosaic key but not the target list.")
                    if target not in missing_targets:
                        missing_targets.append(target)
            
        self._missing_targets = missing_targets

        return()
        
    def _build_whole_target_list(self):
        """
        Create a list where the parts that will go into a linear
        mosaic are gone and only whole targets remain.
        """

        if self._target_dict is None:
            print("I don't have a target dictionary. Can't proceed.")
            return()

        if self._target_list is None:
            self._build_target_list()

        self._whole_target_list = list(self._target_list)

        if self._linearmosaic_dict is None:
            return()

        linearmosaic_targets = self._linearmosaic_dict.keys()
        for this_mosaic in linearmosaic_targets:
            if this_mosaic not in self._whole_target_list:
                self._whole_target_list.append(this_mosaic)
            for this_part in self._linearmosaic_dict[this_mosaic]:
                if this_part in self._whole_target_list:
                    self._whole_target_list.remove(this_part)

        self._whole_target_list.sort()

        if not self._quiet:
            print("Total of "+str(len(self._whole_target_list))+" 'whole' targets.")

        return()

    def _map_targets_to_mosaics(self):
        """
        Create a dictionary that maps from individual targets to mosaics.
        """

        if self._target_dict is None:
            print("I don't have a target dictionary. Can't proceed.")
            return()

        if self._target_list is None:
            self._build_target_list()

        self._mosaic_assign_dict = {}

        if self._linearmosaic_dict is None:
            return()

        linearmosaic_targets = self._linearmosaic_dict.keys()
        for this_mosaic in linearmosaic_targets:
            part_list = self._linearmosaic_dict[this_mosaic]
            for this_part in part_list:
                self._mosaic_assign_dict[this_part] = this_mosaic

        if not self._quiet:
            print("Total of "+str(len(self._mosaic_assign_dict.keys()))+ \
                      " targets assigned to mosaics.")

        return()

#endregion

#region Programs to run checks on the keyHandler and the data

    def check_target_existence(self):
        """
        Print the targets missing a definition in the target list key.
        """

        if self._missing_targets is None:
            print("Missing target list not constructed.")
            
        print("Missing a total of "+str(len(self._missing_targets))+" target definitions.")
        for target in self._missing_targets:
            print("missing: "+target)
            
        return()
    
    def check_ms_existence(self):
        """
        Check that the measurement sets in the ms key all exist in the
        specified directories.
        """

        if self._ms_dict is None:
            return()

        found_count = 0
        missing_count = 0
        for target in self._ms_dict.keys():
            for project_tag in self._ms_dict[target].keys():
                for array_tag in self._ms_dict[target][project_tag].keys():                    
                    found = False
                    local_found_count = 0
                    for ms_root in self._ms_roots:
                        this_ms = ms_root + self._ms_dict[target][project_tag][array_tag]
                        if os.path.isdir(this_ms):
                            found = True
                            found_count += 1
                            local_found_count += 1
                    if local_found_count > 1 and not self._quiet:                        
                        print("Found multiple copies of ms for "+target+" "+project_tag+" "+array_tag)
                    if found:
                        continue
                    missing_count += 1
                    if not self._quiet:
                        print("Missing ms for "+target+" "+project_tag+" "+array_tag)
        
        if not self._quiet:
            print("Verified the existence of "+str(found_count)+" measurement sets. Missing "+str(missing_count)+" measurement set key entries.")
        return()

    def check_singledish_existence(self):
        """
        Check that the FITS files in the singledish key all exist in
        the specified directories.
        """

        if self._singledish_dict is None:
            return()

        found_count = 0
        missing_count = 0
        for target in self._singledish_dict.keys():
            for product in self._singledish_dict[target].keys():
                found = False
                local_found_count = 0
                for this_root in self._singledish_roots:
                    this_fname = this_root + self._singledish_dict[target][product]
                    if os.path.isfile(this_fname):
                        found = True
                        found_count += 1
                        local_found_count += 1
                if local_found_count > 1 and not self._quiet:                        
                    print("Found multiple copies of singledish data for "+target+" "+product)
                if found:
                    continue
                missing_count += 1
                if not self._quiet:
                    print("Missing singledish data for "+target+" "+product)
        
        if not self._quiet:
            print("Verified the existence of "+str(found_count)+" single dish data sets. Missing "+str(missing_count)+" single dish key entries.")
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

    def check_dir_existence(self, imaging=True, postprocess=True):
        """
        Check the existence of the directories for imaging and post-processing.
        """
        
        dir_list = self._targets_for_dir.keys()

        found_dirs=0
        missing_dirs=[]
        for this_dir in dir_list:

            if imaging:
                if os.path.isdir(self._imaging_root+this_dir):
                    found_dirs += 1
                else:
                    if not self._quiet:
                        print("Missing imaging directory :"+self._imaging_root+this_dir)
                    missing_dirs.append(self._imaging_root+this_dir)

            if postprocess:
                if os.path.isdir(self._postprocess_root+this_dir):
                    found_dirs += 1
                else:
                    if not self._quiet:
                        print("Missing post-processing directory :"+self._postprocess_root+this_dir)
                    missing_dirs.append(self._postprocess_root+this_dir)
        
        print("Found "+str(found_dirs)+" directories. Missing "+str(len(missing_dirs))+" directories.")

        return(missing_dirs)
        
#endregion

#region Access data and lists
    
    def _get_dir_for_target(self, target=None, changeto=False, imaging=True):
        """
        Return the imaging working directory given a target name. If
        changeto is true, then change directory to that location.
        """

        if target == None:
            if not self._quiet:
                print("Please specify a target.")
            return(None)

        if target not in self._dir_for_target.keys():
            if not self._quiet:
                print("Target "+target+" is not in the list of targets.")
            return(None)
            
        if imaging:
            this_dir = self._imaging_root + self._dir_for_target[target]+'/'
        else:
            this_dir = self._postprocess_root + self._dir_for_target[target]+'/'

        if changeto:
            if not os.path.isdir(this_dir):
                if not self._quiet:
                    print("You requested to change to "+this_dir+" but it does not exist.")
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
        return(self._get_dir_for_target(target=target, changeto=changeto, imaging=False))

    def get_targets(self, only=None, skip=None, first=None, last=None):
        """
        List the full set of targets.
        """
        this_target_list = \
            utils.select_from_list(self._target_list, first=first, last=last \
                                       , skip=skip, only=only, loose=True)
        return(this_target_list)

    def get_targets_in_ms_key(self, only=None, skip=None, first=None, last=None):
        """
        List all targets that have uv data
        """
        targets = self._target_list
        targets_with_ms = []
        for target in targets:
            if target in self._ms_dict.keys():
                targets_with_ms.append(target)

        this_target_list = \
            utils.select_from_list(targets_with_ms, first=first, last=last \
                                       , skip=skip, only=only, loose=True)
        return(this_target_list)

    def get_linear_mosaic_targets(self, only=None, skip=None, first=None, last=None):
        """
        List all linear mosaics.
        """
        if self._linearmosaic_dict is None:
            return(None)
        mosaic_targets = list(self._linearmosaic_dict.keys())

        this_target_list = \
            utils.select_from_list(mosaic_targets, first=first, last=last \
                                       , skip=skip, only=only, loose=True)
        return(this_target_list)


    def get_whole_targets(self, only=None, skip=None, first=None, last=None):
        """
        List only full galaxy names (no parts, e.g., _1 or _2). Very
        similar to the directory list.
        """
        this_whole_target_list = \
            utils.select_from_list(self._whole_target_list, first=first, last=last \
                                       , skip=skip, only=only, loose=True)
        return(this_whole_target_list)

    def get_interf_configs(self, only=None, skip=None):
        """
        Get a list of interferometer array configurations
        """
        if type(self._config_dict) != type({}):
            return(None)
        if 'interf_config' not in self._config_dict.keys():
            return(None)
        interf_configs = self._config_dict['interf_config'].keys()
        this_list = \
            utils.select_from_list(interf_configs, skip=skip, only=only, loose=True)
        return(this_list)

    def get_feather_configs(self, only=None, skip=None):
        """
        Get a list of feathered single dish + interferometer array configruations.
        """
        if type(self._config_dict) != type({}):
            return(None)
        if 'feather_config' not in self._config_dict.keys():
            return(None)
        feather_configs = self._config_dict['feather_config'].keys()
        this_list = \
            utils.select_from_list(feather_configs, skip=skip, only=only, loose=True)
        return(this_list)

    def get_line_products(self, only=None, skip=None):
        """
        Get a list of line 'products,' i.e., line plus velocity
        resolution combinations.
        """
        if type(self._config_dict) != type({}):
            return(None)
        if 'line_product' not in self._config_dict.keys():
            return(None)
        line_products = self._config_dict['line_product'].keys()
        this_list = \
            utils.select_from_list(line_products, skip=skip, only=only, loose=True)
        return(this_list)

    def get_continuum_products(self, only=None, skip=None):
        """
        Get a list of continuum 'products'.
        """
        if type(self._config_dict) != type({}):
            return(None)
        if 'cont_product' not in self._config_dict.keys():
            return(None)
        cont_products = self._config_dict['cont_product'].keys()
        this_list = \
            utils.select_from_list(cont_products, skip=skip, only=only, loose=True)
        return(this_list)

    def get_parts_for_linear_mosaic(self, target=None):
        """
        Return the parts for a linear mosaic.
        """
        if target == None:
            if not self._quiet:
                print("Please specify a target.")
            return(None)

        if self._linearmosaic_dict is None:
            if not self._quiet:
                print("No linear mosaic dictionary defined.")
            return(None)

        if target not in self._linearmosaic_dict.keys():
            if not self._quiet:
                print("No linear mosaic defined for "+target)
            return(None)
        
        return(self._linearmosaic_dict[target])

    def is_target_in_mosaic(self, target):
        """
        Return true or false depending on whether the target is in a linear mosaic.
        """
        return(target in self._mosaic_assign_dict.keys())

    def set_quiet(self, quiet=True):
        """
        Set the feedback level.
        """
        self._quiet = quiet
        return

    def set_dochecks(self, dochecks=True):
        """
        Set the feedback level.
        """
        self._dochecks = dochecks
        return

#endregion
    
#region Manipulate files and file structure

    def make_missing_directories(self, imaging=False, postprocess=False):
        """
        Make any missing imaging or postprocessing directories.
        """
        
        if not imaging and not postprocess:
            if not self._quiet:
                print("Set either imaging or postprocess to True. Returning.")
            return(False)

        if imaging:
            if not os.path.isdir(self._imaging_root):
                if not self._quiet:
                    print("Missing imaging root directory. Returning.")
                    print("Create: "+self._imaging_root)
                return(False)

        if postprocess:
            if not os.path.isdir(self._postprocess_root):
                if not self._quiet:
                    print("Missing postprocess root directory. Returning.")
                    print("Create: "+self._postprocess_root)
                return(False)

        missing_dirs = self.check_dir_existence(imaging=imaging, postprocess=postprocess)
        made_directories = 0
        for this_missing_dir in missing_dirs:
            made_directories += 1
            os.mkdir(this_missing_dir)
        
        missing_dirs = self.check_dir_existence(imaging=imaging, postprocess=postprocess)

        if not self._quiet:
            print("Made "+str(made_directories)+" directories. Now "+str(len(missing_dirs))+" missing.")

        if len(missing_dirs) == 0:
            return(True)

        return(False)

#endregion
