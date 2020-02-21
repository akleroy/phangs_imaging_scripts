"""
UVDataHandler
    
The PHANGS pipeline to handle staging and pre-processing of uv data
before imaging. Works through a single big class (the
UVDataHandler). This needs to be attached to a keyHandler to access
the target, product, and configuration keys and locations of the uv
data.

There should not be any direct calls to CASA in this
code.. Eventually, this should be able to run without CASA enabled
(though it won't be able to call any of the CASA-specific
routines). Right now, just avoid direct calls to CASA from this class.

To run the individual routines, this code needs to be run inside
CASA. See an example application inside stage_7m_co21.py .

Example:

    $ casa
    from phangsPipeline import keyHandler as kh
    from phangsPipeline import uvdataHandler as uvh
    this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
    this_uvh = uvh.UVDataHandler(key_handler=this_kh)
    this_uvh.copy_data(target = 'ngc3627')

"""

import os, sys, re, shutil
import glob
import numpy as np
import handlerTemplate

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Replace this with a check on load in the future
casa_enabled = True

if casa_enabled:
    import casaVisRoutines as cvr
    #<TODO><DEBUG># 
    reload(cvr)

class UVDataHandler(handlerTemplate.HandlerTemplate):
    """
    Class to manipulate calibrated ALMA visibility data (measurement
    sets), extracting lines, combining multiple data sets, and
    carrying out other steps in prepration for imaging.
    """
    
    ############
    # __init__ #
    ############
    
    def __init__(
            self, 
            key_handler = None,
            dry_run = False,):
        """
        """
        # Can't use super and keep python2/3 agnostic
        handlerTemplate.HandlerTemplate.__init__(self,key_handler = key_handler, dry_run = dry_run)

        
#region 

    ###########################################
    # Define file names for various products. #
    ###########################################

    def _fname_dict(
            self):
        """
        """
        pass
    
    ##########################################
    # Tasks - individual operations on data. #
    ##########################################
    
    def task_copy_data(
            self,
            target = None,
            product = None,
            config = None):
        """
        """
        pass

    def task_run_custom_scripts(
            self,
            target = None,
            product = None,
            config = None):
        """
        """
        pass

    def task_extract_line(
            self,
            target = None,
            product = None,
            config = None):
        """
        """
        pass

    def task_extract_continuum(
            self,
            target = None,
            product = None,
            config = None):
        """
        """
        pass

    def task_concat_uvdata(
            self,
            target = None,
            product = None,
            config = None):
        """
        """
        pass    
    
    ###################################
    # Recipes - combinations of tasks #
    ###################################

    ######################################
    # Loop through all steps and targets #
    ######################################
    
    def loop_postprocess(
        self,
        do_copy=False,
        do_custom=False,
        do_extract_line=False,
        do_extract_cont=False,
        do_concat_line=False,
        do_concat_cont=False,
        make_directories=True,
        ):
        """
        Loops over the full set of targets, products, and configurations
        to run the uv data processing. Toggle the parts of the loop
        using the do_XXX booleans. Other choices affect the algorithms
        used.
        """

        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
 
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)

        if make_directories:
            self._kh.make_missing_directories(imaging=True)        
        
        if do_copy:
                
            for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_interf=True):
                
                pass


        if do_custom:

            for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_interf=True):
                
                pass

        if do_extract_line:

            for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_line=True,just_interf=True):
                
                pass

        if do_extract_cont:

            for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_cont=True,just_interf=True):
                
                pass
        
        if do_concat_line:

            for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_line=True,just_interf=True):
                
                pass

        if do_concat_cont:

            for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_cont=True,just_interf=True):
                
                pass

        return()
        
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
        this_target_multipart_names = self._kh.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self._kh._ms_dict = '+str(self._kh._ms_dict))
            this_ms_dict = self._kh._ms_dict[this_target_ms_name]
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
                    this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target_ms_name)
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
                    for ms_root in self._kh._ms_roots:
                        if os.path.isdir(os.path.join(ms_root, this_ms_data)):
                            this_ms_data_abspath = os.path.abspath(os.path.join(ms_root, this_ms_data))
                    if this_ms_data_abspath == '':
                        logger.error('Could not find the measurement set "'+this_ms_data+'" under keyHandler._ms_roots "'+str(self._kh._ms_roots)+'"!')
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
                        logger.info("")
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
                    cvr.split_science_targets(in_file = this_ms_data_abspath, 
                                              out_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'.ms', 
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
        this_target_multipart_names = self._kh.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get imaging dir
            this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target_ms_name)
            # 
            # find custom staging scripts under config key directory
            scripts_for_this_gal = glob.glob(os.path.join(self._kh._key_dir, 'custom_staging_scripts', this_target_ms_name+'_staging_script.py'))
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
        lines = self._kh.get_line_products(only = just_line) # if just_line is None, then it will return all lines.
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
            logger.info("")
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
            target_width[line] = self._kh._config_dict['line_product'][line]['channel']
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
            self.calculate_chanwidth(
                gal = gal,
                line = line,
                just_array = just_array,
                ext = ext,
                target_width = target_width[line],
                quiet = quiet, 
                overwrite = overwrite, 
                )
            # 
            if interp_to is None or rebin_fac is None:
                logger.warning('Warning! Failed to calculate chanwidth for the line "'+line+'" for galaxy "'+gal+'"!')
                #raise Exception('Error! Failed to calculate chanwidth the line "'+line+'" for galaxy "'+gal+'"!')
            else:
                # 
                #<TODO><20200210># we can move the call of calculate_chanwidth() inside extract_line_for_galaxy(), so that it does not have to be re-run everytime when not overwriting. 
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
        overwrite = False, 
        ):
        """Extract the uv data for a given line for all data sets for a galaxy. 
        """
        
        # 
        # This code is updated from the function with the same name in phangsPipeline.py / imagingPipeline.py.
        # 
        
        # 
        # get vsys and vwidth from key_handler as defined in "target_definitions.txt"
        gal_vsys = self._kh._target_dict[gal]['vsys']
        gal_vwidth = self._kh._target_dict[gal]['vwidth']
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
        this_target_multipart_names = self._kh.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self._kh._ms_dict = '+str(self._kh._ms_dict))
            this_ms_dict = self._kh._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split line data
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get imaging dir
                    this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target_ms_name)
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
                    line_list_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+ext+'_list_lines_in_ms.txt'+append_ext
                    # 
                    # change directory
                    current_dir = os.getcwd()
                    os.chdir(this_imaging_dir)
                    # 
                    # get lines_in_ms for each ms data. Use line_list_file as a cache.
                    if os.path.isfile(line_list_file) and overwrite:
                        os.remove(line_list_file)
                    if not os.path.isfile(line_list_file):
                        lines_in_ms = cvr.list_lines_in_ms(in_file = in_file,
                                                           gal = gal, 
                                                           key_handler = self._kh, 
                                                           quiet = quiet)
                        if lines_in_ms is None:
                            lines_in_ms = []
                        np.savetxt(line_list_file, lines_in_ms, fmt='%s', delimiter=',')
                    else:
                        lines_in_ms = np.genfromtxt(line_list_file, dtype=np.str,  delimiter=',')
                        if lines_in_ms is not None:
                            lines_in_ms = lines_in_ms.tolist()
                    # 
                    if len(lines_in_ms) == 0:
                        logger.warning('No lines found in the measurement set "'+in_file+'".')
                        lines_in_ms = []
                    # 
                    if not (line in lines_in_ms):
                        logger.warning('Line "'+line+'" not found in the measurement set "'+in_file+'".')
                    else:
                        cvr.extract_line(in_file = in_file, 
                                         out_file = out_file, 
                                         line = line, 
                                         vsys = vsys, 
                                         vwidth = vwidth, 
                                         gal = gal, 
                                         key_handler = self._kh, 
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
    
    
    ################
    # concat_lines #
    ################
    
    def concat_lines(   
        self, 
        gal, 
        just_array = None, 
        just_line = None, 
        ext = '', 
        quiet = False, 
        overwrite = False, 
        ):
        """Concatenate the extracted lines into a few aggregated measurement sets.
        """
        
        # 
        # This code is updated from the function concat_phangs_lines() in phangsPipeline.py / imagingPipeline.py.
        # 
        # dzliu: I renamed the argument 'lines' to 'just_line'
        # dzliu: I updated the just_array behavior. It can be a str, a list, or None. It can be '7m', ['12m','7m'], or '12m+7m'. 
        # dzliu: I renamed concat_line_for_gal() to concat_line_for_galaxy()
        # dzliu: I added overwrite option, when calling concat_line_for_galaxy().
        # 
        
        # 
        # make sure just_array is a list. If it is None, we will process both '7m' and '12m'.
        if just_array is not None:
            if np.isscalar(just_array):
                just_array = [just_array]
        
        # 
        # make sure just_line is a list
        if just_line is not None:
            if np.isscalar(just_line):
                just_line = [just_line]
        
        # 
        # get lines
        lines = self._kh.get_line_products(only = just_line) # if just_line is None, then it will return all lines.
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
            logger.info("")
            logger.info("--------------------------------------------------------")
            logger.info("START: Concatenating spectral line measurements.")
            logger.info("--------------------------------------------------------")
            logger.info("Galaxy: "+gal)
        # 
        # Loop and extract lines for each data set
        for line in lines: 
            
            ## Unless we just do the 12m, we build a 7m dataset
            #if just_array != '12m':
            # Build 7m dataset
            if (just_array is None) or (np.any([re.match(r'(\b|_)7m(\b|_)', t, re.IGNORECASE) for t in just_array])):
                self.concat_line_for_galaxy(
                    gal = gal,
                    just_array = '7m',
                    tag = '7m',
                    line = line,
                    do_chan0 = True, 
                    overwrite = overwrite, 
                    )
            
            ## Unless we just do the 7m, we build a 12m dataset
            #if just_array != '7m':
            # Build 12m dataset
            if (just_array is None) or (np.any([re.match(r'(\b|_)12m(\b|_)', t, re.IGNORECASE) for t in just_array])):
                self.concat_line_for_galaxy(
                    gal = gal,
                    just_array = '12m',
                    tag = '12m',
                    line = line,
                    do_chan0 = True, 
                    overwrite = overwrite, 
                    )
            
            # This can probably be improved, but works for now. Check if
            # we lack either 12m or 7m data, in which case there is no
            # combined data set to make.
            
            has_7m = (len(glob.glob(gal+'*7m*'+line+'*')) > 0)
            has_12m = (len(glob.glob(gal+'*12m*'+line+'*')) > 0)
            if (not has_12m) or (not has_7m):
                logger.warning('Missing 12m or 7m for gal "'+gal+'" and line "'+line+'" ... no combined set made.')
                continue
            
            # Build combined 12m+7m data
            if (just_array is None) or (np.any([re.match(r'(\b|_)7m(\b|_)', t, re.IGNORECASE) for t in just_array]) and \
                                        np.any([re.match(r'(\b|_)12m(\b|_)', t, re.IGNORECASE) for t in just_array])):
                self.concat_line_for_galaxy(
                    gal = gal,
                    just_array = None,
                    tag = '12m+7m',
                    line = line,
                    do_chan0 = True, 
                    overwrite = overwrite, 
                    )
        
        if not quiet:
            logger.info("--------------------------------------------------------")
            logger.info("END: Concatenating spectral line measurements.")
            logger.info("--------------------------------------------------------")
    
    
    
    ################ ##########################
    # concat_lines # # concat_line_for_galaxy #
    ################ ##########################
    
    def concat_line_for_galaxy(
        self, 
        gal, 
        line, 
        just_proj = None,
        just_ms = None,
        just_array = None, 
        tag = '',
        do_chan0 = True,
        quiet = False, 
        overwrite = False, 
        ):
        """Concatenate the uv data for a given line for all data sets for a galaxy. 
        """
        
        # 
        # This code is updated from the function concat_line_for_gal() in phangsPipeline.py / imagingPipeline.py.
        # 
        
        
        # 
        # get imaging dir
        this_gal_imaging_dir = self._kh.get_imaging_dir_for_target(gal)
        logger.info('Imaging dir: '+this_gal_imaging_dir)
        # 
        # prepare the output file name
        if tag != '':
            out_file =  os.path.join(this_gal_imaging_dir, gal+'_'+tag+'_'+line+'.ms')
            chan0_vis = os.path.join(this_gal_imaging_dir, gal+'_'+tag+'_'+line+'_chan0.ms')
        else:
            out_file =  os.path.join(this_gal_imaging_dir, gal+'_'+line+'.ms')
            chan0_vis = os.path.join(this_gal_imaging_dir, gal+'_'+line+'_chan0.ms')
        # 
        # check existing output data
        if do_chan0:
            if os.path.isdir(out_file) and os.path.isdir(chan0_vis) and not overwrite:
                logger.info('Found existing output data "'+os.path.basename(out_file)+'" and "'+os.path.basename(chan0_vis)+'", will not overwrite them.')
                return
        else:
            if os.path.isdir(out_file) and not overwrite:
                logger.info('Found existing output data "'+os.path.basename(out_file)+'", will not overwrite it.')
                return
        # 
        # prepare the concatenating file list
        files_to_concat = []
        # 
        # check multipart names for the input galaxy
        this_target_multipart_names = self._kh.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self._kh._ms_dict = '+str(self._kh._ms_dict))
            this_ms_dict = self._kh._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split line data
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get imaging dir
                    this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target_ms_name)
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
                    # get extracted line data as this_in_file
                    this_in_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'_'+line+'.ms'
                    this_in_file = os.path.join(this_imaging_dir, this_in_file)
                    # 
                    # check this_in_file
                    if not os.path.isdir(this_in_file):
                        logger.warning('Warning! The line uv data measurement set "'+this_in_file+'" was not found! It should be extracted in previous steps.')
                    else:
                        # record this_in_file
                        files_to_concat.append(this_in_file)
        # 
        # end of for each this_target_multipart_names
        if len(files_to_concat) == 0:
            logger.info('No files to concatenate for gal "'+gal+'" and line "'+line+'". Returning.')
            return
        else:
            logger.info('Concatenating '+str(len(files_to_concat))+' measurement sets for gal "'+gal+'" and line "'+line+'"')
            logger.debug('files_to_concat: '+str(files_to_concat))
        # 
        # change directory to imaging dir
        current_dir = os.getcwd()
        os.chdir(this_gal_imaging_dir)
        # 
        # Concatenate all of the relevant files
        cvr.concat_ms(in_file_list = files_to_concat, 
                      out_file = out_file, 
                      do_chan0 = do_chan0, 
                      quiet = quiet, 
                      overwrite = overwrite, 
                      )
        # 
        # change dir back
        os.chdir(current_dir)
        # 

        
    
    
    ################# #######################
    # extract_lines # # calculate_chanwidth #
    ################# #######################
    
    def calculate_chanwidth(
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
        overwrite = False, 
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
        this_target_multipart_names = self._kh.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self._kh._ms_dict = '+str(self._kh._ms_dict))
            this_ms_dict = self._kh._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split line data
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get imaging dir
                    this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target_ms_name)
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
                    this_chanwidth_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+ext+'_chanwidth_for_line_'+line+'.txt'
                    # 
                    # change directory
                    current_dir = os.getcwd()
                    os.chdir(this_imaging_dir)
                    # 
                    # get chanwidth for each ms data. Use this_chanwidth_file as a cache.
                    #overwrite = True #<TODO><DEBUG># 
                    if os.path.isfile(this_chanwidth_file) and overwrite:
                        os.remove(this_chanwidth_file)
                    if not os.path.isfile(this_chanwidth_file):
                        this_chanwidth = cvr.chanwidth_for_line(in_file = this_vis,
                                                                line = line,
                                                                gal = gal, 
                                                                key_handler = self._kh, 
                                                                quiet = quiet)
                        if this_chanwidth is None:
                            this_chanwidth = []
                        np.savetxt(this_chanwidth_file, this_chanwidth, delimiter=',')
                    else:
                        this_chanwidth = np.loadtxt(this_chanwidth_file, delimiter=',', ndmin=1)
                        if this_chanwidth is not None:
                            this_chanwidth = this_chanwidth.tolist()
                    # 
                    # change dir back
                    os.chdir(current_dir)
                    # 
                    if this_chanwidth is None:
                        continue
                    # 
                    #logger.debug('type(this_chanwidth) = '+str(type(this_chanwidth))+', this_chanwidth = '+str(this_chanwidth)) #<DEBUG>#
                    if len(this_chanwidth) == 0:
                        continue
                    # 
                    # record all chanwidths
                    for chanwidth in this_chanwidth:
                        chanwidth_list.append(chanwidth)
                    vis_list.append(this_vis)
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
        # Get galaxy vsys and vwidth from self._kh
        gal_vsys = self._kh._target_dict[gal]['vsys']
        gal_vwidth = self._kh._target_dict[gal]['vwidth']
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
                logger.info('... ' + str(chanwidth_list[ii]) + ' ... ' + str(vis_list[ii]))
            logger.info("... max is " + str(max_cw))
            logger.info("... min is " + str(min_cw))
            logger.info("... interpolate_to " + str(interpolate_cw))
            logger.info("... then rebin by " + str(rebin_fac))
            logger.info("... to final " + str(rebin_fac*interpolate_cw))
        # 
        # Return
        return interpolate_cw, rebin_fac
    
    
    
    #####################
    # extract_continuum #
    #####################
    
    def extract_continuum(
        self, 
        gal, 
        vsys = None, 
        vwidth = None, 
        just_array = None, 
        ext = '', 
        append_ext = '', 
        do_statwt = True, 
        quiet = False, 
        overwrite = False, 
        ):
        """
        Split continuum uv data from the input uv data measurement set. 
        
        Line channels will be identified automatically.
        """
        
        # 
        # This code is updated from the function extract_phangs_continuum() in phangsPipeline.py / imagingPipeline.py.
        # 
        # 
        
        # Best practice here regarding statwt isn't obvious - it's the
        # continuum, so there are no signal free channels. I think we just
        # have to hope that the signal does not swamp the noise during the
        # statwt or consider turning off the statwt in high S/N continuum
        # cases.
        
        # 
        # get lines to flag
        #lines_to_flag = line_list.lines_co+line_list.lines_13co+line_list.lines_c18o
        lines_to_flag = ['co', '13co','c18o'] # dzliu: the line_list module is not imported here, but inside casaVisRoutines.py. So here we provide a str list here, then deal with it inside casaVisRoutines.py.
        # 
        # print starting message
        if not quiet:
            logger.info("")
            logger.info("--------------------------------------------------------")
            logger.info("START: Extracting continuum from data set.")
            logger.info("--------------------------------------------------------")
        # 
        # extract continuum data
        self.extract_continuum_for_galaxy(
            gal = gal, 
            lines_to_flag = lines_to_flag, 
            vsys = vsys, 
            vwidth = vwidth, 
            just_array = just_array, 
            ext = ext, 
            append_ext = append_ext, 
            do_statwt = do_statwt, 
            do_collapse = True, 
            quiet = quiet, 
            overwrite = overwrite, 
            )
        # 
        # print ending message
        if not quiet:
            logger.info("--------------------------------------------------------")
            logger.info("END: Extracting continuum from data set.")
            logger.info("--------------------------------------------------------")
    
    
    
    ##################### ################################
    # extract_continuum # # extract_continuum_for_galaxy #
    ##################### ################################
    
    def extract_continuum_for_galaxy(
        self, 
        gal, 
        lines_to_flag, 
        vsys = None, 
        vwidth = None, 
        just_proj = None,
        just_ms = None,
        just_array = None, 
        ext = None, 
        append_ext = None, 
        do_statwt = True, 
        do_collapse = True, 
        quiet = False, 
        overwrite = False, 
        ):
        """Extract the continuum uv data for all data sets for a galaxy. 
        """
        
        # 
        # This code is updated from the function with the same name in phangsPipeline.py / imagingPipeline.py.
        # 
        
        # 
        # get vsys and vwidth from key_handler as defined in "target_definitions.txt"
        gal_vsys = self._kh._target_dict[gal]['vsys']
        gal_vwidth = self._kh._target_dict[gal]['vwidth']
        # 
        # if user has not input a vsys, then we use what is defined in the key_handler
        if vsys is None:
            vsys = gal_vsys
        # 
        # if user has not input a vwidth, then we use what is defined in the key_handler
        if vwidth is None:
            vwidth = gal_vwidth
        # 
        # check multipart names for the input galaxy
        this_target_multipart_names = self._kh.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self._kh._ms_dict = '+str(self._kh._ms_dict))
            this_ms_dict = self._kh._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split line data
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get imaging dir
                    this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target_ms_name)
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
                    out_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'_cont.ms'
                    # 
                    # change directory
                    current_dir = os.getcwd()
                    os.chdir(this_imaging_dir)
                    # 
                    # get candidate lines for each ms data
                    #lines_in_ms = cvr.list_lines_in_ms(in_file = in_file,
                    #                                   gal = gal, 
                    #                                   key_handler = self._kh, 
                    #                                   quiet = quiet)
                    # 
                    # extract_continuum
                    cvr.extract_continuum(in_file = in_file, 
                                          out_file = out_file, 
                                          lines_to_flag = lines_to_flag, 
                                          vsys = vsys, 
                                          vwidth = vwidth, 
                                          gal = gal, 
                                          key_handler = self._kh, 
                                          do_statwt = do_statwt, 
                                          do_collapse = do_collapse, 
                                          quiet = quiet, 
                                          overwrite = overwrite, 
                                         )
                    # 
                    # change dir back
                    os.chdir(current_dir)
    
    
    
    ####################
    # concat_continuum #
    ####################
    
    def concat_continuum(   
        self, 
        gal, 
        just_array = None, 
        quiet = False, 
        overwrite = False, 
        ):
        """Concatenate continuum data sets into a single file for one galaxy or part of galaxy.
        """
        
        # 
        # This code is updated from the function concat_phangs_continuum() in phangsPipeline.py / imagingPipeline.py.
        # 
        # dzliu: I updated the just_array behavior. It can be a str, a list, or None. It can be '7m', ['12m','7m'], or '12m+7m'. 
        # dzliu: I renamed concat_continuum_for_gal() to concat_continuum_for_galaxy()
        # dzliu: I added overwrite option, when calling concat_continuum_for_galaxy().
        # 
        
        # 
        # make sure just_array is a list. If it is None, we will process both '7m' and '12m'.
        if just_array is not None:
            if np.isscalar(just_array):
                just_array = [just_array]
        
        # 
        # print starting message
        if not quiet:
            logger.info("")
            logger.info("--------------------------------------------------------")
            logger.info("START: Concatenating continuum from data set.")
            logger.info("--------------------------------------------------------")
            logger.info("Galaxy: "+gal)
        #
        
        ## Unless we just do the 12m, we build a 7m dataset
        #if just_array != '12m':
        # Build 7m dataset
        if (just_array is None) or (np.any([re.match(r'(\b|_)7m(\b|_)', t, re.IGNORECASE) for t in just_array])):
            self.concat_continuum_for_galaxy(
                gal = gal,
                just_array = '7m',
                tag = '7m',
                overwrite = overwrite, 
                )
        
        ## Unless we just do the 7m, we build a 12m dataset
        #if just_array != '7m':
        # Build 12m dataset
        if (just_array is None) or (np.any([re.match(r'(\b|_)12m(\b|_)', t, re.IGNORECASE) for t in just_array])):
            self.concat_continuum_for_galaxy(
                gal = gal,
                just_array = '12m',
                tag = '12m',
                overwrite = overwrite, 
                )
        
        # This can probably be improved, but works for now. Check if
        # we lack either 12m or 7m data, in which case there is no
        # combined data set to make.
        
        has_7m = (len(glob.glob(gal+'*7m*'+'_cont'+'*')) > 0)
        has_12m = (len(glob.glob(gal+'*12m*'+'_cont'+'*')) > 0)
        if (not has_12m) or (not has_7m):
            logger.warning("Missing 12m or 7m for gal "+gal+" continuum ... no combined set made.")
        else:
            # Build combined 12m+7m data
            if (just_array is None) or (np.any([re.match(r'(\b|_)7m(\b|_)', t, re.IGNORECASE) for t in just_array]) and \
                                        np.any([re.match(r'(\b|_)12m(\b|_)', t, re.IGNORECASE) for t in just_array])):
                self.concat_continuum_for_galaxy(
                    gal = gal,
                    just_array = None,
                    tag = '12m+7m',
                    overwrite = overwrite, 
                    )
        
        if not quiet:
            logger.info("--------------------------------------------------------")
            logger.info("END: Concatenating continuum from data set.")
            logger.info("--------------------------------------------------------")
    
    
    
    #################### ###############################
    # concat_continuum # # concat_continuum_for_galaxy #
    #################### ###############################
    
    def concat_continuum_for_galaxy(
        self, 
        gal, 
        just_proj = None,
        just_ms = None,
        just_array = None, 
        tag = '',
        quiet = False, 
        overwrite = False, 
        ):
        """Concatenate the uv data for a given line for all data sets for a galaxy. 
        """
        
        # 
        # This code is updated from the function concat_cont_for_gal() in phangsPipeline.py / imagingPipeline.py.
        # 
        
        # 
        # get imaging dir
        this_gal_imaging_dir = self._kh.get_imaging_dir_for_target(gal)
        logger.info('Imaging dir: '+this_gal_imaging_dir)
        # 
        # prepare the output file name
        if tag != '':
            out_file =  os.path.join(this_gal_imaging_dir, gal+'_'+tag+'_cont'+'.ms')
        else:
            out_file =  os.path.join(this_gal_imaging_dir, gal+'_cont'+'.ms')
        # 
        # check existing output data
        if os.path.isdir(out_file) and not overwrite:
            logger.info('Found existing output data "'+os.path.basename(out_file)+'", will not overwrite it.')
            return
        # 
        # prepare the concatenating file list
        files_to_concat = []
        # 
        # check multipart names for the input galaxy
        this_target_multipart_names = self._kh.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self._kh._ms_dict = '+str(self._kh._ms_dict))
            this_ms_dict = self._kh._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split line data
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get imaging dir
                    this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target_ms_name)
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
                    # get extracted line data as this_in_file
                    this_in_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'_cont'+'.ms'
                    this_in_file = os.path.join(this_imaging_dir, this_in_file)
                    # 
                    # check this_in_file
                    if not os.path.isdir(this_in_file):
                        logger.warning('Warning! The previously extracted continuum uv data measurement set "'+this_in_file+'" was not found!')
                    else:
                        # record this_in_file
                        files_to_concat.append(this_in_file)
        # 
        # end of for each this_target_multipart_names
        if len(files_to_concat) == 0:
            logger.info('No files to concatenate for gal "'+gal+'" continuum. Returning.')
            return
        # 
        # change directory to imaging dir
        current_dir = os.getcwd()
        os.chdir(this_gal_imaging_dir)
        # 
        # Concatenate all of the relevant files
        cvr.concat_ms(in_file_list = files_to_concat, 
                      out_file = out_file, 
                      do_chan0 = False, 
                      quiet = quiet, 
                      overwrite = overwrite, 
                      )
        # 
        # change dir back
        os.chdir(current_dir)
        # 
    
    
    
    
    ###################
    # cleanup_staging #
    ###################
    
    def cleanup_staging(
        self, 
        gal, 
        just_proj = None, 
        just_ms = None, 
        just_array = None, 
        just_line = None, 
        ext = '', 
        append_ext = '', 
        including_cache = False, 
        quiet = False,  
        ):
        
        
        # 
        # get lines
        lines = self._kh.get_line_products(only = just_line) # if just_line is None, then it will return all lines.
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
            logger.info("")
            logger.info("--------------------------------------------------------")
            logger.info("START: Clean up staging.")
            logger.info("--------------------------------------------------------")
        # 
        # check multipart names for the input galaxy
        this_target_multipart_names = self._kh.get_parts_for_linmos(gal)
        if this_target_multipart_names is None:
            this_target_multipart_names = [gal]
        # 
        # loop each multipart name of each galaxy. If this galaxy has no multipart, it is just its galaxy name.
        for this_target_ms_name in this_target_multipart_names:
            # 
            # get ms dict
            #logger.debug('self._kh._ms_dict = '+str(self._kh._ms_dict))
            this_ms_dict = self._kh._ms_dict[this_target_ms_name]
            if len(this_ms_dict) == 0:
                logger.error('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
                raise Exception('The target '+this_target_ms_name+' does not have a valid ms key? Please check your ms_file_key.txt!')
            # 
            # loop proj_tag and array_tag, and split line data
            for this_proj_tag in this_ms_dict.keys():
                for this_array_tag in this_ms_dict[this_proj_tag].keys():
                    # 
                    # get imaging dir
                    this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target_ms_name)
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
                    # change directory
                    current_dir = os.getcwd()
                    os.chdir(this_imaging_dir)
                    # 
                    # check ms
                    temp_ext_list = ['.temp', '.temp.flagversions', '.temp2', '.temp2.flagversions']
                    check_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+ext+'.ms'+append_ext
                    if os.path.isdir(check_file):
                        logger.info('Keeping "'+check_file+'"')
                    else:
                        pass #<TODO># what if data not found?
                    for temp_ext in temp_ext_list:
                        check_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+ext+'.ms'+append_ext+temp_ext
                        if os.path.isdir(check_file):
                            logger.info('Cleaning "'+check_file+'"')
                            shutil.rmtree(check_file)
                    # 
                    # Loop lines for each data set and check line ms
                    for line in lines: 
                        check_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'_'+line+'.ms'
                        if os.path.isdir(check_file):
                            logger.info('Keeping "'+check_file+'"')
                        else:
                            pass #<TODO># what if data not found?
                        for temp_ext in temp_ext_list:
                            check_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'_'+line+'.ms'+temp_ext
                            if os.path.isdir(check_file):
                                logger.info('Cleaning "'+check_file+'"')
                                shutil.rmtree(check_file)
                    # 
                    # check continuum ms
                    check_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+'_cont'+'.ms'
                    if os.path.isdir(check_file):
                        logger.info('Keeping "'+check_file+'"')
                    else:
                        pass #<TODO># what if data not found?
                    for temp_ext in temp_ext_list:
                        check_file = this_target_ms_name+'_'+this_proj_tag+'_'+this_array_tag+ext+'.ms'+append_ext+temp_ext
                        if os.path.isdir(check_file):
                            logger.info('Cleaning "'+check_file+'"')
                            shutil.rmtree(check_file)
                    # 
                    # check chanwidth_for_line cache <TODO>
                    # 
                    # check list_lines_in_ms cache <TODO>
                    # 
                    # change dir back
                    os.chdir(current_dir)
        # 
        # print ending message
        if not quiet:
            logger.info("--------------------------------------------------------")
            logger.info("END: Clean up staging.")
            logger.info("--------------------------------------------------------")

# region Loops

    ########################
    # loop to stage uvdata #
    ########################
    
    # was stage_imaging
    def loop_stage_uvdata(
        self, 
        do_copy = True, 
        do_split = True, 
        do_statwt = False, 
        do_custom_scripts = True,
        do_extract_lines = True, 
        do_extract_cont = True, 
        do_concat_lines = True, 
        do_concat_cont = True, 
        do_cleanup = True, 
        use_symlink = True, 
        make_directories=True,
        ): 
        """
        Loops over the full set of targets, products, and
        configurations to run the uv data staging. Toggle the parts of
        the loop using the do_XXX booleans. Other choices affect the
        algorithms used.
        """
        
        if self._targets_list is None:            
            logger.error("Need a target list.")
            return(None)
 
        if self._all_products is None:            
            logger.error("Need a products list.")
            return(None)

        if make_directories:
            self._kh.make_missing_directories(imaging=True)

        # Loop within each operation. This ensures that any
        # cross-linking between different targets, configurations,
        # etc. is handled.
        
        if do_copy:
            
            for this_target in self._targets_list:

                for this_product in self._all_products():
                    
                    for this_config in self._interf_configs_list:

                        # This needs to be refactored to a recipe,
                        # passing the target, config, product.

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

        if do_custom_scripts:

            for this_target in self._targets_list:

                for this_product in self._all_products():
                    
                    for this_config in self._interf_configs_list:

                        # This needs to be refactored to a recipe,
                        # passing the target, config, product.

                        # Optionally, run custom scripts at this stage. This could, for
                        # example, flag data or carry out uv continuum subtraction. The
                        # line and continuum extensions defined here point the subsequent
                        # programs at the processed data.
                        self.custom_scripts(
                            gal = gal, 
                            )

        if do_extract_lines:

            for this_target in self._targets_list:

                for this_product in self._all_products():
                    
                    for this_config in self._interf_configs_list:

                        # This needs to be refactored to a recipe,
                        # passing the target, config, product.


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

        if do_extract_cont:

            for this_target in self._targets_list:

                for this_product in self._all_products():
                    
                    for this_config in self._interf_configs_list:

                        # This needs to be refactored to a recipe,
                        # passing the target, config, product.

                        # Extract the continuum, avoiding lines and averaging all
                        # frequencies in each SPW together. This step also uses statwt to
                        # empirically weight the data.
                        self.extract_continuum(
                            gal = gal, 
                            just_array = just_array, 
                            do_statwt = True,
                            quiet = quiet, 
                            overwrite = overwrite, 
                            )

        if do_concat_lines:

            for this_target in self._targets_list:

                for this_product in self._all_products():
                    
                    for this_config in self._interf_configs_list:
                        
                        # This needs to be refactored to a recipe,
                        # passing the target, config, product.

                        # Concatenate the extracted lines into the measurement sets that
                        # we will use for imaging. This step also makes a "channel 0"
                        # measurement for each line.
                        self.concat_lines(
                            gal = gal, 
                            just_array = just_array, 
                            just_line = just_line, 
                            quiet = quiet, 
                            overwrite = overwrite, 
                            )

        if do_concat_cont:

            for this_target in self._targets_list:

                for this_product in self._all_products():
                    
                    for this_config in self._interf_configs_list:
                        
                        # This needs to be refactored to a recipe,
                        # passing the target, config, product.

                        self.concat_continuum(
                            gal = gal, 
                            just_array = just_array, 
                            quiet = quiet, 
                            overwrite = overwrite, 
                            )
                        
        if do_cleanup:

            for this_target in self._targets_list:

                for this_product in self._all_products():
                    
                    for this_config in self._interf_configs_list:
                        
                        # This needs to be refactored to a recipe,
                        # passing the target, config, product.

                        # Remove intermediate files. The big space-savers here are the
                        # initial copies of the data. The data after frequency averaging
                        # are smaller by a large factor (~10). For reference, re-copying
                        # all of the PHANGS-ALMA LP takes less than a day on the OSU
                        # system. Full line and continuum exraction takes longer. 

                        self.cleanup_staging(
                            gal = gal, 
                            just_array = just_array, 
                            quiet = quiet, 
                            )
                        
        return()

        # end of loop_stage

    
    




