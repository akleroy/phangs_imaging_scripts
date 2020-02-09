"""uvdataHandler

This module copies uv data from ms_dir to imaging dir and prepare things for imaging. 

This code needs to be run inside CASA. 

Example:
    $ casa
    from phangsPipeline import keyHandler as kh
    from phangsPipeline import uvdataHandler as uvh
    this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
    this_uvh = uvh.uvDataHandler(key_handler=this_kh)
    this_uvh.copy_data(target = 'ngc3627')
    
"""

import os, sys, re, shutil
import glob
import numpy as np
import casaCubeRoutines as ccr
import casaMosaicRoutines as cmr
import casaFeatherRoutines as cfr
import casaVisRoutines as cvr
reload(cvr) #<TODO><DEBUG># 

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = True

class uvDataHandler:
    """
    Class to manage ALMA calibrated uv data, aka. Measurement Sets. 
    """
    
    def __init__(
        self, 
        key_handler, 
        ):
        
        # check key_handler input
        try:
            key_handler._key_dir
            key_handler._imaging_root
            key_handler._ms_keys
            key_handler._dir_keys
            key_handler._target_keys
        except:
            logger.error('The input key_handler is invalid! Please check your master_key.txt and ms_file_key.txt!')
            raise Exception('Please input a valid key handler!')
        
        # prepare imaging root directory
        if not os.path.isdir(key_handler._imaging_root):
            logger.debug('Imaging root directory does not exist, creating '+key_handler._imaging_root)
            os.makedirs(key_handler._imaging_root)
            if not os.path.isdir(key_handler._imaging_root):
                logger.error('Failed to create the imaging root directory '+key_handler._imaging_root+'! Please check your file system writing permission!')
                raise Exception('Failed to create the imaging root directory! Please check your file system writing permission!')
        
        # check key_handler validity
        if not key_handler.check_key_existence():
            raise Exception('Incomplete key handler! Please check the error message above!')
        
        # store key_handler
        self.key_handler = key_handler
        
        # stage uv data for imaging
        #self.stage_imaging()
    
    # 
    def stage_imaging(
        self, 
        targets = None, 
        do_copy = True, 
        do_custom_scripts = True,
        do_extract_lines = True, 
        do_concat_lines = True, 
        do_extract_cont = True, 
        do_concat_cont = True, 
        do_cleanup = True, 
        ): 
        """
        This function stages the data before imaging process. 
        This includes: copying uv data, splitting science uv data, splitting line uv data.
        """
        
        # debug checks
        logger.debug('self.key_handler.get_targets() = '+str(self.key_handler.get_targets()))
        logger.debug('self.key_handler.get_targets_in_ms_key() = '+str(self.key_handler.get_targets_in_ms_key()))
        
        # 
        if targets is None:
            this_targets = self.key_handler.get_targets()
            #this_targets = self.key_handler.get_targets_in_ms_key()
        else:
            if np.isscalar(targets):
                targets = [targets]
            this_targets = self.key_handler.get_targets(only = targets)
            #this_targets = self.key_handler.get_targets_in_ms_key(only = targets)
        
        # check targets
        if len(this_targets) == 0:
            logger.error('No target found! Input targets = '+str(targets))
            raise Exception('No target found! Input targets = '+str(targets))
        
        # loop targets
        for this_target in this_targets:
            # 
            # renaming
            gal = this_target
            # 
            # copy uv data for each galaxy
            if do_copy:
                self.copy_data(gal = gal)
            # 
            # Optionally, run custom scripts at this stage. This could, for
            # example, flag data or carry out uv continuum subtraction. The
            # line and continuum extensions defined here point the subsequent
            # programs at the processed data.
            if do_custom_scripts:
                scripts_for_this_gal = glob.glob(os.path.join(self.key_handler._key_dir, 'custom_staging_scripts', gal+'_staging_script.py'))
                for this_script in scripts_for_this_gal:
                    execfile(this_script) #<TODO># 
            # 
            # Extract lines, includes regridding and rebinning to the velocity
            # grid specified in the text file keys. Runs statwt afterwards,
            # the result is a bunch of line-only data sets but still
            # execution-by-execution.
        
            if do_extract_lines:
                self.extract_phangs_lines(
                    gal=gal,
                    just_array=just_array,
                    quiet=False,
                    append_ext=line_ext,
                    lines=lines)
        
            # Concatenate the extracted lines into the measurement sets that
            # we will use for imaging. This step also makes a "channel 0"
            # measurement for each line.
        
            if do_concat_lines:
                self.concat_phangs_lines(
                    gal=gal,
                    just_array=just_array,
                    quiet=False,
                    lines=lines)
        
            # Extract the continuum, avoiding lines and averaging all
            # frequencies in each SPW together. This step also uses statwt to
            # empirically weight the data.
        
            if do_extract_cont:
                self.extract_phangs_continuum(
                    gal=gal,
                    just_array=just_array,
                    quiet=False,
                    do_statwt=True,
                    append_ext=cont_ext)
        
            if do_concat_cont:
                self.concat_phangs_continuum(
                    gal=gal,
                    just_array=just_array,
                    quiet=False)
            
            # Remove intermediate files. The big space-savers here are the
            # initial copies of the data. The data after frequency averaging
            # are smaller by a large factor (~10). For reference, re-copying
            # all of the PHANGS-ALMA LP takes less than a day on the OSU
            # system. Full line and continuum exraction takes longer. 
            
            if do_cleanup:
                self.cleanup_phangs_staging(
                    gal=gal,
                    just_array=just_array)
                
        # end of stage_imaging()
    
    
    # 
    def copy_data(self, 
                  gal = None,
                  just_proj = None,
                  just_array = None,
                  just_ms = None,
                  do_split = True,
                  do_statwt = False,
                  use_symlink = True, 
                  overwrite = False, 
                  quiet = False):
        """
        Copies ALMA measurement set (ms) data from its original location, which is specified in a
        text file ms_file_key.txt, to the imaging directory. 
        
        Then splits out only the science target.
        """
        
        # 
        # This code is adapted from phangsPipeline.py / imagingPipeline.py.
        # 
        # just_proj is like "886" indicating which ALMA project
        # just_array is like "12m" without suffix
        # just_ms is like "12m_1" with a suffix
        
        # 
        # check input galaxy
        if gal is None:
            logger.error("Please specify a galaxy to copy the uv data.")
            raise Exception('Please specify a galaxy to copy the uv data.')
            return
        
        # 
        # make sure just_proj is a list
        if just_proj is not None:
            if np.isscalar(just_proj):
                just_proj = [just_proj]
        
        # 
        # make sure just_array is a list
        if just_array is not None:
            if np.isscalar(just_array):
                just_array = [just_array]
        
        # 
        # make sure just_ms is a list
        if just_ms is not None:
            if np.isscalar(just_ms):
                just_ms = [just_ms]
        
        # 
        # check multipart names for the input galaxy
        #logger.debug('self.key_handler.get_parts_for_linmos(gal) = '+str(self.key_handler.get_parts_for_linmos(gal)))
        this_target_multipart_names = self.key_handler.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each target, if this target has no multipart, it is just its target name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # check ms dict
            #logger.debug('self.key_handler._ms_dict = '+str(self.key_handler._ms_dict))
            this_ms_dict = self.key_handler._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split science target ms data into the imaging dir
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    this_ms_data = this_ms_dict[this_proj_tag][this_array_tag]
                    this_imaging_dir = self.key_handler.get_imaging_dir_for_target(this_target_ms_name)
                    logger.debug('Target '+this_target_ms_name+', proj '+this_proj_tag+', array '+this_array_tag+', ms data '+this_ms_data+', imaging dir '+this_imaging_dir)
                    # 
                    # do some variable renaming
                    this_proj = this_proj_tag
                    this_array = re.sub(r'^(.*?)_([0-9]+)$', r'\1', this_array_tag)
                    this_ms = this_array_tag # to be compatible with original phangsPipeline.py
                    # 
                    # check user input just_proj, just_array and just_ms
                    if just_proj is not None:
                        if not (this_proj_tag in just_proj):
                            continue
                    if just_array is not None:
                        if not (this_array in just_array):
                            continue
                    if just_ms is not None:
                        if not (this_ms in just_ms):
                            continue
                    # 
                    # find ms data absolute path
                    this_ms_data_abspath = ''
                    for ms_root in self.key_handler._ms_roots:
                        if os.path.isdir(os.path.join(ms_root, this_ms_data)):
                            this_ms_data_abspath = os.path.abspath(os.path.join(ms_root, this_ms_data))
                    if this_ms_data_abspath == '':
                        logger.error('Could not find the measurement set "'+this_ms_data+'" under keyHandler._ms_roots "'+str(self.key_handler._ms_roots)+'"!')
                        raise Exception('Could not find the measurement set! Please check your ms_root in master_key.txt and the ms_file_key.txt!')
                    # 
                    # check imaging directory
                    if not os.path.isdir(this_imaging_dir):
                        logger.debug('Creating imaging directory '+this_imaging_dir)
                        os.makedirs(this_imaging_dir)
                        if not os.path.isdir(this_imaging_dir):
                            logger.error('Failed to create the imaging directory '+this_imaging_dir+'! Please check your file system writing permission!')
                            raise Exception('Failed to create the imaging directory! Please check your file system writing permission!')
                    # 
                    # change directory
                    current_dir = os.getcwd()
                    os.chdir(this_imaging_dir)
                    # 
                    # start copying data
                    if not quiet:
                        logger.info("--------------------------------------------------------")
                        logger.info("START: Copying the original data.")
                        logger.info("--------------------------------------------------------")
                        #logger.info("Galaxy: " + this_target_ms_name)
                        #logger.info("Project: " + this_proj_tag)
                        #logger.info("Array: " + this_array_tag)
                        #logger.info("Measurements Set: " + this_ms) this_ms = this_array_tag + suffix
                        logger.info("Outputting to " + this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'.ms')
                    
                    # 
                    # copy data and split science targets
                    cvr.split_science_targets(in_ms = this_ms_data_abspath, 
                                              out_ms = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'.ms', 
                                              do_split = do_split, 
                                              do_statwt = do_statwt, 
                                              use_symlink = use_symlink, 
                                              overwrite = overwrite, 
                                              quiet = quiet )
                    
                    # 
                    # change dir back
                    os.chdir(current_dir)
                    
                    # 
                    # print ending message
                    if not quiet:
                        logger.info("--------------------------------------------------------")
                        logger.info("END: Copying data from original location.")
                        logger.info("--------------------------------------------------------")
                # 
                # end for array
            # 
            # end for proj
        # 
        # end for this_target_ms_name
    # 
    # end of copy_data
    
    #
    #<TODO><20200209># def extract_phangs_lines
    #<TODO><20200209># def concat_phangs_lines
    #<TODO><20200209># def extract_phangs_continuum
    #<TODO><20200209># def concat_phangs_continuum
    #<TODO><20200209># def cleanup_phangs_staging




