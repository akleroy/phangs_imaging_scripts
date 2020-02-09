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
        ): 
        
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
            # copy uv data for each galaxy
            self.copy_data(gal = this_target_ms_name, 
                           just_proj = this_proj_tag, 
                           just_array = this_array_tag, 
                           just_ms = this_ms_data, 
                           )
                
    # 
    def copy_data(gal = None,
                  just_proj = None,
                  just_array = None,
                  just_ms = None,
                  do_split = True,
                  do_statwt = False,
                  data_dirs = None,
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
        # initialize data_dirs 
        if data_dirs is None:
            data_dirs = ['']
        
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
                    this_array = re.sub(r'^(.*?)_([0-9]+)$', r'\1', this_array_tag)
                    this_ms = this_array_tag # to be compatible with original phangsPipeline.py
                    logger.debug('Target '+this_target_ms_name+', proj '+this_proj_tag+', array '+this_array_tag+', ms data '+this_ms_data+', imaging dir '+this_imaging_dir)
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
                        logger.info("Galaxy: ", gal)
                        logger.info("Project: ", this_proj_tag)
                        logger.info("Array: ", this_array)
                        logger.info("Measurements Set: ", this_ms)
                    # 
                    # prepare the copied data folder name (copied_file is a data folder).
                    # If we are going to do some additional processing, make this an intermediate file ("_copied"). 
                    if do_split:
                        copied_file = gal+'_'+this_proj+'_'+this_ms+'_copied.ms'
                    else:
                        copied_file = gal+'_'+this_proj+'_'+this_ms+'.ms'
                    # 
                    # check existing copied data in the imaging directory
                    if os.path.isdir(copied_file):
                        if not overwrite:
                            logger.warning('Found existing copied data '+copied_file+', will not re-copy it.')
                            continue
                        else:
                            shutil.rmtree(copied_file)
                            if os.path.isdir(copied_file+'.flagversions'):
                                shutil.rmtree(copied_file+'.flagversions')
                    # 
                    # Copy. We could place a symbolic link here using ln -s
                    # instead, but instead I think the right move is to make
                    # the intermediate files and then clean them up. This
                    # avoids "touching" the original data at all.
                    command = 'cp -Lr '+in_file+' '+copied_file
                    print command
                    var = os.system(command)    
                    print var
                    
                    command = 'cp -Lr '+in_file+'.flagversions'+' '+copied_file+'.flagversions'
                    print command
                    var = os.system(command) 
                    print var
                    
                    # Call split and statwt if desired.

                    if do_split:

                        if quiet == False:
                            print "Splitting out science target data."

                        out_file = gal+'_'+this_proj+'_'+this_ms+'.ms'

                        os.system('rm -rf '+out_file)
                        os.system('rm -rf '+out_file+'.flagversions')
                        
                        # If present, we use the corrected column. If not,
                        # then we use the data column.

                        mytb = au.createCasaTool(tbtool)
                        mytb.open(copied_file)
                        colnames = mytb.colnames()
                        if colnames.count('CORRECTED_DATA') == 1:
                            print "Data has a CORRECTED column. Will use that."
                            use_column = 'CORRECTED'
                        else:
                            print "Data lacks a CORRECTED column. Will use DATA column."
                            use_column = 'DATA'
                        mytb.close()

                        split(vis=copied_file
                              , intent ='OBSERVE_TARGET#ON_SOURCE'
                              , datacolumn=use_column
                              , outputvis=out_file)        

                        os.system('rm -rf '+copied_file)
                        os.system('rm -rf '+copied_file+'.flagversions')

                    if do_statwt:

                        if quiet == False:
                            print "Using statwt to re-weight the data."

                        statwt(vis=out_file,
                               datacolumn='DATA')
                    
                    # 
                    # change dir back
                    os.chdir(current_dir)
                    
                    # 
                    # print ending message
                    if not quiet:
                        print "--------------------------------------------------------"
                        print "END: Copying data from original location."
                        print "--------------------------------------------------------"
                    

        
        
        
    # 
    # 

