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
    
    ############
    # __init__ #
    ############
    
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
    
    
    #################
    # stage_imaging #
    #################
    
    def stage_imaging(
        self, 
        targets = None, 
        do_copy = True, 
        do_split = True, 
        do_statwt = False, 
        use_symlink = True, 
        do_custom_scripts = True,
        do_extract_lines = True, 
        do_concat_lines = True, 
        do_extract_cont = True, 
        do_concat_cont = True, 
        do_cleanup = True, 
        just_proj = None, 
        just_array = None, 
        just_ms = None, 
        just_line = None, 
        overwrite = False, 
        quiet = False, 
        ): 
        """
        This function stages the data before imaging process. 
        This includes: copying uv data, splitting science uv data, splitting line uv data.
        """
        
        # debug checks
        logger.debug('self.key_handler.get_targets() = '+str(self.key_handler.get_targets()))
        logger.debug('self.key_handler.get_targets_in_ms_key() = '+str(self.key_handler.get_targets_in_ms_key()))
        
        # check input targets, if none then process all targets, otherwise make sure it is a list
        if targets is not None:
            if np.isscalar(targets):
                targets = [targets]
        gals = self.key_handler.get_targets(only = targets) # if input targets is None, it will return all targets.
        
        # check targets
        if len(gals) == 0:
            logger.error('No galaxy found with the input targets '+str(targets)+'!')
            raise Exception('No galaxy found with the input targets '+str(targets)+'!')
        
        # 
        # make sure just_proj is a list or None
        if just_proj is not None:
            if np.isscalar(just_proj):
                just_proj = [just_proj]
        
        # 
        # make sure just_array is a list or None
        if just_array is not None:
            if np.isscalar(just_array):
                just_array = [just_array]
        
        # 
        # make sure just_ms is a list or None
        if just_ms is not None:
            if np.isscalar(just_ms):
                just_ms = [just_ms]
        # 
        # make sure just_line is a list or None
        if just_line is not None:
            if np.isscalar(just_line):
                just_line = [just_line]
        
        # 
        # loop targets
        for gal in gals:
            # 
            # copy uv data for each galaxy
            if do_copy:
                self.copy_data(
                    gal = gal, 
                    just_proj = just_proj, 
                    just_array = just_array, 
                    just_ms = just_ms, 
                    do_split = do_split, 
                    do_statwt = do_statwt, 
                    use_symlink = use_symlink, 
                    overwrite = overwrite, 
                    quiet = quiet, 
                    )
            # 
            # 
            if do_custom_scripts:
                # Optionally, run custom scripts at this stage. This could, for
                # example, flag data or carry out uv continuum subtraction. The
                # line and continuum extensions defined here point the subsequent
                # programs at the processed data.
                self.custom_scripts(
                    gal = gal, 
                    )
            # 
            # 
            if do_extract_lines:
                # Extract lines, includes regridding and rebinning to the velocity
                # grid specified in the text file keys. Runs statwt afterwards,
                # the result is a bunch of line-only data sets but still
                # execution-by-execution.
                self.extract_lines(
                    gal = gal, 
                    just_array = just_array, 
                    just_line = just_line, 
                    quiet = quiet, 
                    overwrite = overwrite, 
                    )
            # 
            # 
            if do_concat_lines:
                # Concatenate the extracted lines into the measurement sets that
                # we will use for imaging. This step also makes a "channel 0"
                # measurement for each line.
                raise NotImplementedError()
                self.concat_phangs_lines(
                    gal=gal,
                    just_array=just_array,
                    quiet=quiet,
                    lines=just_line)
            # 
            # 
            if do_extract_cont:
                # Extract the continuum, avoiding lines and averaging all
                # frequencies in each SPW together. This step also uses statwt to
                # empirically weight the data.
                raise NotImplementedError()
                self.extract_phangs_continuum(
                    gal=gal,
                    just_array=just_array,
                    quiet=quiet,
                    do_statwt=True,
                    append_ext=cont_ext)
            # 
            # 
            if do_concat_cont:
                raise NotImplementedError()
                self.concat_phangs_continuum(
                    gal=gal,
                    just_array=just_array,
                    quiet=quiet)
            # 
            # 
            if do_cleanup:
                # Remove intermediate files. The big space-savers here are the
                # initial copies of the data. The data after frequency averaging
                # are smaller by a large factor (~10). For reference, re-copying
                # all of the PHANGS-ALMA LP takes less than a day on the OSU
                # system. Full line and continuum exraction takes longer. 
                self.cleanup_phangs_staging(
                    gal=gal,
                    just_array=just_array)
                
        # end of stage_imaging()
    
    
    #############
    # copy_data #
    #############
    
    def copy_data(
        self, 
        gal,
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
        # This code is adapted from the function with the same name in phangsPipeline.py / imagingPipeline.py.
        # 
        # just_proj is like "886" indicating which ALMA project
        # just_array is like "12m" without suffix
        # just_ms is like "12m_1" with a suffix
        # 
        # TODO 20200209 dzliu: gal is a scalar or list?
        
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
        this_target_multipart_names = self.key_handler.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self.key_handler._ms_dict = '+str(self.key_handler._ms_dict))
            this_ms_dict = self.key_handler._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split science target ms data into the imaging dir
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get ms data file path
                    this_ms_data = this_ms_dict[this_proj_tag][this_array_tag]
                    # 
                    # get imaging dir
                    this_imaging_dir = self.key_handler.get_imaging_dir_for_target(this_target_ms_name)
                    # 
                    # print debug info
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
    # end of copy_data()
    
    
    ##################
    # custom_scripts #
    ##################
    
    def custom_scripts(
        self, 
        gal, 
        quiet = False, 
        ): 
        """
        Optionally, run custom scripts at this stage. This could, for
        example, flag data or carry out uv continuum subtraction. The
        line and continuum extensions defined here point the subsequent
        programs at the processed data.
        """
        
        # 
        # check multipart names for the input galaxy
        this_target_multipart_names = self.key_handler.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get imaging dir
            this_imaging_dir = self.key_handler.get_imaging_dir_for_target(this_target_ms_name)
            # 
            # find custom staging scripts under config key directory
            scripts_for_this_gal = glob.glob(os.path.join(self.key_handler._key_dir, 'custom_staging_scripts', this_target_ms_name+'_staging_script.py'))
            # 
            # print starting message
            if not quiet:
                logger.info("--------------------------------------------------------")
                logger.info("START: Running custom staging scripts.")
                logger.info("--------------------------------------------------------")
                logger.info("Custom staging scripts: "+str(scripts_for_this_gal))
                logger.info("Imaging directory: "+str(this_imaging_dir))
            # 
            # change directory
            current_dir = os.getcwd()
            os.chdir(this_imaging_dir)
            # 
            # run custom staging scripts
            for this_script in scripts_for_this_gal:
                execfile(this_script) #<TODO># 
            # 
            # change dir back
            os.chdir(current_dir)
            # 
            # print ending message
            if not quiet:
                logger.info("--------------------------------------------------------")
                logger.info("END: Running custom staging scripts.")
                logger.info("--------------------------------------------------------")
    
    
    #################
    # extract_lines #
    #################
    
    def extract_lines(
        self, 
        gal, 
        just_array = None, 
        just_line = None, 
        ext = '', 
        append_ext = '', 
        quiet = False, 
        overwrite = False, 
        ):
        """
        Split line uv data from the input uv data measurement set. 
        """
        
        # 
        # This code is updated from the function extract_phangs_lines() in phangsPipeline.py / imagingPipeline.py.
        # 
        # dzliu: I renamed the argument 'lines' to 'just_line'
        # 
        
        # 
        # get lines
        lines = self.key_handler.get_line_products(only = just_line) # if just_line is None, then it will return all lines.
        if len(lines) == 0:
            if just_line is not None:
                logger.error('Error! Could not find the input line "'+str(just_line)+'" in the line names as defined in "config_definitions.txt"!')
                raise Exception('Error! Could not find the input line "'+str(just_line)+'" in the line names as defined in "config_definitions.txt"!')
            else:
                logger.error('Error! Could not find lines! Please check "config_definitions.txt"!')
                raise Exception('Error! Could not find lines! Please check "config_definitions.txt"!')
        # 
        # print starting message
        if not quiet:
            logger.info("--------------------------------------------------------")
            logger.info("START: Extracting spectral lines from data set.")
            logger.info("--------------------------------------------------------")
        # 
        # Loop and extract lines for each data set
        for line in lines:    
            # 
            # set target_width as defined by users in "config_definitions.txt"
            target_width = {}
            target_width['co21'] = 2.5
            target_width['13co21'] = 2.5
            target_width['c18o21'] = 6.0
            target_width[line] = self.key_handler._config_dict['line_product'][line]['channel']
            # 
            # set edge_for_statwt <TODO>
            edge_for_statwt = {}
            edge_for_statwt[line] = 25 # default
            edge_for_statwt['co21'] = 25
            edge_for_statwt['13co21'] = 25
            edge_for_statwt['c18o21'] = 20
            # 
            # calculate phangs chanwidth
            interp_to, rebin_fac = \
            self.calculate_phangs_chanwidth(
                gal = gal,
                line = line,
                just_array = just_array,
                ext = ext,
                target_width = target_width[line],
                quiet = quiet,
                )
            # 
            if interp_to is None or rebin_fac is None:
                logger.warning('Warning! Failed to extract the line "'+line+'" for galaxy "'+gal+'"!')
                raise Exception('Error! Failed to extract the line "'+line+'" for galaxy "'+gal+'"!')
                #continue
            # 
            # extract line data
            self.extract_line_for_galaxy(
                gal = gal, 
                line = line, 
                just_array = just_array, 
                ext = ext, 
                append_ext = append_ext, 
                chan_fine = interp_to, 
                rebin_factor = rebin_fac, 
                edge_for_statwt = edge_for_statwt[line],
                quiet = quiet, 
                overwrite = overwrite, 
                )
        # 
        # print ending message
        if not quiet:
            logger.info("--------------------------------------------------------")
            logger.info("END: Extracting spectral lines from data set.")
            logger.info("--------------------------------------------------------")
    
    
    ################# ##############################
    # extract_lines # # calculate_phangs_chanwidth #
    ################# ##############################
    
    def calculate_phangs_chanwidth(
        self, 
        gal, 
        line, 
        just_proj = None, 
        just_ms = None, 
        just_array = None, 
        ext = '', 
        append_ext = '', 
        target_width = 2.5, 
        quiet = False, 
        ):
        """Determine the channel width to use when splitting line data from a measurment set.
        """
        
        # 
        # This code is updated from the function with the same name in phangsPipeline.py / imagingPipeline.py.
        # 
        
        # 
        # set constants
        one_plus_eps = 1.0+1e-3
        
        # 
        # Initialize an empty list
        chanwidth_list = []
        vis_list = []
        
        # 
        # check multipart names for the input galaxy
        this_target_multipart_names = self.key_handler.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self.key_handler._ms_dict = '+str(self.key_handler._ms_dict))
            this_ms_dict = self.key_handler._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split line data
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get imaging dir
                    this_imaging_dir = self.key_handler.get_imaging_dir_for_target(this_target_ms_name)
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
                    # get copied ms data
                    this_vis = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+ext+'.ms'+append_ext
                    # 
                    # change directory
                    current_dir = os.getcwd()
                    os.chdir(this_imaging_dir)
                    # 
                    # get chanwidth for each ms data
                    this_chanwidth = cvr.chanwidth_for_line(in_ms = this_vis,
                                                            line = line,
                                                            gal = gal, 
                                                            key_handler = self.key_handler, 
                                                            quiet = quiet)
                    # 
                    if this_chanwidth is None:
                        continue
                    # 
                    # record all chanwidths
                    for chanwidth in this_chanwidth:
                        chanwidth_list.append(chanwidth)
                    vis_list.append(this_vis)
                    # 
                    # change dir back
                    os.chdir(current_dir)
        # 
        # No line found in any ms data? <TODO>
        if len(chanwidth_list) == 0:
            return None, None
        # 
        # Calculate the least common channel
        chanwidths = np.array(chanwidth_list)
        max_cw = np.max(chanwidths)
        min_cw = np.min(chanwidths)
        interpolate_cw = max_cw*one_plus_eps
        # 
        # Get the mosaic parameters for comparison
        #mosaic_parms = read_mosaic_key() 
        #if mosaic_parms.has_key(gal):
        #    vsys = mosaic_parms[gal]['vsys']
        #    vwidth = mosaic_parms[gal]['vwidth']
        # 
        # Get galaxy vsys and vwidth from self.key_handler
        gal_vsys = self.key_handler._target_dict[gal]['vsys']
        gal_vwidth = self.key_handler._target_dict[gal]['vwidth']
        # 
        # Rebinning factor
        rat = target_width / interpolate_cw
        rebin_fac = int(round(rat))
        if rebin_fac < 1:
            rebin_fac = 1
        # 
        if not quiet:
            logger.info("")
            logger.info("For galaxy: "+gal+" and line "+line)
            logger.info("... channel widths:")
            for ii in range(len(vis_list)):
                logger.info(str(chanwidth_list[ii]) + ' ... ' + str(vis_list[ii]))
            logger.info("... max is " + str(max_cw))
            logger.info("... min is " + str(min_cw))
            logger.info("... interpolate_to " + str(interpolate_cw))
            logger.info("... then rebin by " + str(rebin_fac))
            logger.info("... to final " + str(rebin_fac*interpolate_cw))
        # 
        # Return
        return interpolate_cw, rebin_fac
    
    
    ################# ###########################
    # extract_lines # # extract_line_for_galaxy #
    ################# ###########################
    
    def extract_line_for_galaxy(
        self, 
        gal, 
        line, 
        vsys = None, 
        vwidth = None, 
        just_proj = None,
        just_ms = None,
        just_array = None, 
        ext = None, 
        append_ext = None, 
        chan_fine = 0.5, 
        rebin_factor = 5, 
        edge_for_statwt = -1,
        do_statwt = True, 
        quiet = False, 
        overwrite = overwrite, 
        ):
        """Extract the uv data for a given line for all data sets for a galaxy. 
        """
        
        # 
        # This code is updated from the function with the same name in phangsPipeline.py / imagingPipeline.py.
        # 
        
        # 
        # get vsys and vwidth from key_handler as defined in "target_definitions.txt"
        gal_vsys = self.key_handler._target_dict[gal]['vsys']
        gal_vwidth = self.key_handler._target_dict[gal]['vwidth']
        # 
        # if user has not input a vsys, then we use what is defined in the key_handler
        if vsys is None:
            vsys = gal_vsys
        elif not np.isclose(vsys, gal_vsys):
            # if user has input a vsys, use it instead of the one in the key_handler, but report warning if the values are different
            logger.warning('Warning! User has input a vsys of '+str(vsys)+' km/s which is different from the vsys of '+str(gal_vsys)+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
        # 
        # if user has not input a vwidth, then we use what is defined in the key_handler
        if vwidth is None:
            vwidth = gal_vwidth
        elif not np.isclose(vwidth, gal_vwidth):
            # if user has input a vwidth, use it instead of the one in the key_handler, but report warning if the values are different
            logger.warning('Warning! User has input a vwidth of '+str(vwidth)+' km/s which is different from the vwidth of '+str(gal_vwidth)+' km/s for the galaxy "'+gal+'" in the key_handler as defined in "target_definitions.txt"!')
        # 
        # check multipart names for the input galaxy
        this_target_multipart_names = self.key_handler.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self.key_handler._ms_dict = '+str(self.key_handler._ms_dict))
            this_ms_dict = self.key_handler._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split line data
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get imaging dir
                    this_imaging_dir = self.key_handler.get_imaging_dir_for_target(this_target_ms_name)
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
                    # get copied ms data as in_file
                    in_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+ext+'.ms'+append_ext
                    out_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'_'+line+'.ms'
                    # 
                    # change directory
                    current_dir = os.getcwd()
                    os.chdir(this_imaging_dir)
                    # 
                    # get chanwidth for each ms data
                    lines_in_ms = cvr.list_lines_in_ms(in_ms = this_vis,
                                                       gal = gal, 
                                                       key_handler = self.key_handler, 
                                                       quiet = quiet)
                    # 
                    if lines_in_ms is None:
                        logger.warning('No lines found in the measurement set "'+in_file+'".')
                        continue
                    # 
                    if not (line in lines_in_ms):
                        logger.warning('Line "'+line+'" not found in the measurement set "'+in_file+'".')
                        continue
                    # 
                    cvr.extract_line(in_file = in_file, 
                                     out_file = out_file, 
                                     line = line, 
                                     vsys = vsys, 
                                     vwidth = vwidth, 
                                     gal = gal, 
                                     key_handler = self.key_handler, 
                                     chan_fine = chan_fine, 
                                     rebin_factor = rebin_factor, 
                                     do_statwt = do_statwt, 
                                     edge_for_statwt = edge_for_statwt, 
                                     quiet = quiet, 
                                     overwrite = overwrite, 
                                    )
                    # 
                    # change dir back
                    os.chdir(current_dir)
    
    
    
    
    
    # 
    # 
    # 
    def concat_phangs_lines(   
        gal=None,
        just_array='',
        ext='',
        quiet=False,
        lines=['co21', 'c18o21'],
        ):
        """
        Concatenate the extracted lines into a few aggregated measurement
        sets.
        """

        if quiet == False:
            print "--------------------------------------------------------"
            print "START: Concatenating spectral line measurements."
            print "--------------------------------------------------------"
            print ""
            print "Galaxy: "+gal

        for line in lines:    

            # Unless we just do the 12m, we build a 7m dataset
            if just_array != '12m':
                concat_line_for_gal(
                    gal=gal,
                    just_array = '7m',
                    tag='7m',
                    line=line,
                    do_chan0=True)

            # Unless we just do the 7m, we build a 12m dataset
            if just_array != '7m':
                concat_line_for_gal(
                    gal=gal,
                    just_array = '12m',
                    tag='12m',
                    line=line,
                    do_chan0=True)

            # This can probably be improved, but works for now. Check if
            # we lack either 12m or 7m data, in which case there is no
            # combined data set to make.

            has_7m = len(glob.glob(gal+'*7m*'+line+'*')) > 0
            has_12m = len(glob.glob(gal+'*12m*'+line+'*')) > 0
            if has_12m == False or has_7m == False:
                print "Missing 12m or 7m ... no combined set made."
                continue

            if just_array == '' or just_array == None:
                concat_line_for_gal(
                    gal=gal,
                    just_array = None,
                    tag='12m+7m',
                    line=line,
                    do_chan0=True)

        if quiet == False:
            print "--------------------------------------------------------"
            print "END: Concatenating spectral line measurements."
            print "--------------------------------------------------------"
    
    #
    #<TODO><20200209># def extract_phangs_lines
    #<TODO><20200209># def concat_phangs_lines
    #<TODO><20200209># def extract_phangs_continuum
    #<TODO><20200209># def concat_phangs_continuum
    #<TODO><20200209># def cleanup_phangs_staging




