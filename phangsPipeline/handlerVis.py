"""
UVDataHandler
    
The PHANGS pipeline to handle staging and pre-processing of uv data
before imaging. Works through a single big class (the
UVDataHandler). This needs to be attached to a keyHandler to access
the target, product, and configuration keys and locations of the uv
data.

To run the individual routines, this code needs to be run inside
CASA. See an example application inside stage_7m_co21.py .

Example:

    $ casa
    from phangsPipeline import handlerKeys as kh
    from phangsPipeline import uvdataHandler as uvh
    this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
    this_uvh = uvh.UVDataHandler(key_handler = this_kh, dry_run = False)
    # Set which data to process
    this_uvh.set_line_products(only=['co21'])
    this_uvh.set_interf_configs(only=['12m+7m'])
    this_uvh.set_targets(only=['ngc3621'])
    # Process
    this_uvh.loop_stage_uvdata()

"""

import os, sys, re, shutil
import glob
import numpy as np

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = (sys.argv[0].endswith('start_casa.py')) #<TODO># check whether we are inside CASA environment

if casa_enabled:
    logger.debug('casa_enabled = True')
    import casaVisRoutines as cvr
    reload(cvr) #<TODO><DEBUG># 
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import handlerTemplate

try:
    import utilsFilenames as fnames
except ImportError:
    from phangsPipeline import utilsFilenames as fnames


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
        
    def get_all_input_ms(
        self, 
        target=None, 
        just_for_configs=None
        ):
        """

        Get all the absolute file path(s) to the input measurement
        sets for the input target. Optionally, restrict to array_tags
        that overlap the configurations provided in
        just_for_configs. This is used
        
        Here we need to input target, config and project tag. These should be specified in "ms_file_key.txt". 
        
        For example, `get_ms_filepaths(target='ngc3621_1', config='7m')` returns 
        `filenames = ['ngc3621_1_886_7m_1.ms', 'ngc3621_1_886_7m_2.ms', ...]` and 
        `filepaths = ['/path/to/*886*/sci*/group*/mem*/calibrated/uid*.ms.split.cal', ..., ...]`
        """
        
        ms_filenames = []
        ms_filepaths = []
        for this_target, this_project, this_arraytag, this_obsnum in \
            self.get_ms_projects_and_obsnum_for_target(target=target, config=config):
            # 
            ms_filepath = self.get_input_ms_file(
                target=this_target, 
                project=this_project, 
                arraytag=this_arraytag, 
                obsnum=this_obsnum)

            # append to output lists
            ms_filenames.append(ms_filename) # without '.ms' suffix
            ms_filepaths.append(ms_filepath)
        # output
        return ms_filenames, ms_filepaths


    def _fname_dict(
            self, 
            target=None,
            config=None,
            product=None,
            all_ms_data=False,
            extra_ext='',
            ):
        """
        Internal function to provide file names for the input target, config and product. 
        """
        
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Error checking
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        if target is None:
            logger.error("Need a target.")
            return()

        if config is None:
            logger.error("Need a config.")
            return()

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Initialize
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        fname_dict = {}

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        # Set data file/dir names
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        
        if all_ms_data:
            this_ms_filenames, this_ms_filepaths = self._kh.get_ms_filenames_and_filepaths(target = target, config = config)
            fname_dict['ms_filepaths'] = this_ms_filepaths
            fname_dict['ms_filenames'] = [t+extra_ext+'.ms' for t in this_ms_filenames]
            if product is not None:
                fname_dict['ms_extracted'] = [t+'_'+product+extra_ext+'.ms' for t in this_ms_filenames]
        
        if product is not None:
            fname_dict['ms_concatenated'] = target+'_'+config+'_'+product+extra_ext
        
        return fname_dict
        
    ######################################
    # Loop through all steps and targets #
    ######################################
    
    def loop_stage_uvdata(
        self,
        do_copy = True,
        do_concat_line = True,
        do_concat_cont = True,
        do_custom = False,
        do_contsub = False, 
        do_extract_line = True,
        do_extract_cont = True,
        make_directories = True,
        extra_ext = '',         
        just_projects=None,        
        overwrite = False, 
        ):
        """
        Loops over the full set of targets, products, and configurations
        to run the uv data processing. Toggle the parts of the loop
        using the do_XXX booleans. Other choices affect the algorithms
        used.
        """
        
        if make_directories:
            self._kh.make_missing_directories(imaging=True)
        
                
        target_list = self.get_targets()
        product_list = self.get_products()
        config_list = self.get_interf_configs()

        # This loop goes over target, array, observation number, and
        # project and prepares each observation for imaging and
        # analysis.

        for this_target, this_project, this_array_tag, this_obsnum in \
                self._kh.loop_over_input_ms(targets=target_list,
                                            configs=config_list,
                                            projects=just_projects):

                for this_product in product_list:
                    
                    if do_copy:

                        self.task_copy_and_split(
                            target = this_target, 
                            project = this_project, 
                            array_tag = this_array_tag, 
                            obsnum = this_obsnum, 
                            product = this_product,
                            extra_ext = extra_ext, 
                            # could add algorithm flags here
                            overwrite = overwrite, 
                            )

        # This loop goes over target, product, and configuration and
        # combines all data for each target + product + config
        # combination.

        for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_line=True,just_interf=True):

                if do_concat_line:

                    self.task_split_and_concat(
                        target = this_target, 
                        config = this_config, 
                        product = this_product, 
                        extra_ext = extra_ext,                         
                        overwrite = overwrite, 
                        )
                
                if do_concat_cont:

                    self.task_split_and_concat(
                        target = this_target, 
                        config = this_config, 
                        product = this_product, 
                        extra_ext = extra_ext, 
                        overwrite = overwrite, 
                        )
        
        # This loop processes the concatenated data

        for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_line=True,just_interf=True):

                if do_contsub:
                    
                    self.task_contsub(
                        target = this_target, 
                        project = this_project, 
                        array_tag = this_array_tag, 
                        obsnum = this_obsnum, 
                        extra_ext = extra_ext, 
                        # could add algorithm flags here
                        overwrite = overwrite, 
                        )
        

        if do_extract_line:

            for this_product in self.get_line_products():
            
                common_chan_width = self.task_compute_common_channel(
                    product = this_product, extra_ext = extra_ext
                    )

                        self.task_extract_line(
                            target = this_target, 
                            project = this_project, 
                            array_tag = this_array_tag, 
                            obsnum = this_obsnum, 
                            product = this_product, 
                            extra_ext = extra_ext,
                            # could add algorithm flags here
                            overwrite = overwrite, 
                            )

                if do_extract_cont:

                    for this_product in self.get_cont_products():

                        self.task_extract_continuum(
                            target = this_target, 
                            project = this_project, 
                            array_tag = this_array_tag, 
                            obsnum = this_obsnum, 
                            product = this_product, 
                            extra_ext = extra_ext,
                            # could add algorithm flags here
                            overwrite = overwrite, 
                            )                
                
        return()
        
    ##########################################
    # Tasks - individual operations on data. #
    ##########################################
    
    def task_copy_and_split(
            self, 
            target = None,
            project = None,
            array_tag = None,
            obsnum = None,
            product = None,
            extra_ext_out = '',
            do_split = True, 
            do_statwt = False, 
            use_symlink = True, 
            overwrite = False, 
            ):
        """
        Copy visibility data for one target, project, array_tag,
        obsnum combination from their original location to the imaging
        directory for the target. Then optionally split out only the
        science targets.
        """
        
        if target is None:
            logging.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if project is None:
            logging.error("Please specify a project.")
            raise Exception("Please specify a project.")

        if array_tag is None:
            logging.error("Please specify an array_tag.")
            raise Exception("Please specify an array_tag.")

        if obsnum is None:
            logging.error("Please specify an obsnum.")
            raise Exception("Please specify an obsnum.")
        
        infile = self._kh.get_file_for_input_ms(
            target=target, project=project, array_tag=array_tag, obsnum=obsnum)
        
        field = self._kh.get_field_for_input_ms(
            target=target, project=project, array_tag=array_tag, obsnum=obsnum)
        if (field.lower()).strip() == 'all':
            field = ''

        outfile = fnames.get_staged_msname(
            target=target, project=project, array_tag=array_tag, obsnum=obsnum, 
            product=product, ext=extra_ext_out)

        # If requested, select on SPW for the product

        spw = ''
        if product is not None:

            if product is in self._kh.get_line_products:

                vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target)

                spw = cvr.find_spws_for_line(infile = infile, 
                                             line = this_line, 
                                             vsys = vsys, vwidth = vwidth)

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Copying u-v data for "+outfile)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        # Change to the imaging directory for the target

        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        
        if not self._dry_run and casa_enabled:

            cvr.split_science_targets(
                infile = infile, 
                outfile = outfile,  
                split_field = field,
                split_spw = spw,
                do_split = do_split, 
                do_statwt = do_statwt, 
                use_symlink = use_symlink, 
                overwrite = overwrite, 
                )

        return()

    def task_concat_uvdata(
            self, 
            target = None, 
            product = None, 
            config = None, 
            extra_ext_in = '', 
            extra_ext_out = '', 
            overwrite = False, 
            ):

        """
        Concatenate all measurement sets for the supplied
        target+config+product combination.
        """
        
        if target is None:
            logging.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if product is None:
            logging.error("Please specify a product.")
            raise Exception("Please specify a product.")

        if config is None:
            logging.error("Please specify a config.")
            raise Exception("Please specify a config.")
        
         
        # get the fname dict for the input target and config

        # Come back to this list generation part.

        fname_dict = self._fname_dict(target = target, config = config, product = product, all_ms_data = True, extra_ext = extra_ext)
        ms_concatenated = fname_dict['ms_concatenated']
        ms_list_to_concatenate = fname_dict['ms_extracted']
        ms_list_to_concatenate = sorted(ms_list_to_concatenate)
        str_of_ms_list_to_concatenate = ', '.join(['"'+t+'"' for t in ms_list_to_concatenate]) # for printing
        logger.debug('Current directory: "'+os.getcwd()+'"')
        logger.debug('Concatenating: "'+ms_concatenated+'" <-- ['+str_of_ms_list_to_concatenate+']')

        # Revise here

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Splitting (by spw) then concatenating u-v data for:")
        logger.info("... target: "+target)
        logger.info("... product: "+product)
        logger.info("... config: "+config)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        # Change to the imaging directory for the target

        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        
        if not self._dry_run and casa_enabled:         

            # Add a for loop over ms to make temporary files with only
            # relevant spws (because concat includes no selection)
            
            cvr.concat_ms(in_file_list = ms_list_to_concatenate, 
                          out_file = ms_concatenated, 
                          overwrite = overwrite, 
                          )

            # Clean up those temporary measurement sets.

        return()
    
    def task_contsub(
            self, 
            target = None, 
            project = None, 
            array_tag = None, 
            obsnum = None, 
            product = None,
            extra_ext = '', 
            overwrite = False, 
            ):
        """
        Run u-v plane continuum subtraction on individual input measurement sets.
        """

        if target is None:
            logging.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if project is None:
            logging.error("Please specify a project.")
            raise Exception("Please specify a project.")

        if array_tag is None:
            logging.error("Please specify an array_tag.")
            raise Exception("Please specify an array_tag.")

        if obsnum is None:
            logging.error("Please specify an obsnum.")
            raise Exception("Please specify an obsnum.")
        
        infile = fnames.get_staged_msname(
            target=target, project=project, array_tag=array_tag, obsnum=obsnum, 
            product=product, ext=extra_ext)

        # get target vsys and vwidth
        vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target)

        # get lines to flag as defined in keys
        lines_to_flag = self._kh.get_lines_to_flag_for_continuum_product(product=product)
        
        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("u-v continuum subtraction for "+infile)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        # Change to the imaging directory for the target

        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        
        if not self._dry_run and casa_enabled:

            cvr.contsub(in_file = infile, 
                        lines_to_flag = lines_to_flag, 
                        vsys = vsys, 
                        vwidth = vwidth, 
                        overwrite = overwrite, 
                        )

        return()
    
    def task_run_custom_scripts(
            self, 
            target = None, 
            product = None, 
            config = None, 
            extra_ext = '', 
            ):
        """
        """
        pass


    def task_compute_common_channel(
            self, 
            target = None, 
            product = None, 
            config = None, 
            extra_ext = '', 
            ):
        """
        Compute the coarsest channel width among all ms data for the input target, config and product. 
        """

        # 
        logger.info('START: Computing common channel width among all ms data for target '+target+', config '+config+' and product '+product+'.')

        # 
        # get imaging dir and change directory
        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)

        # 
        # get target vsys and vwidth
        vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target)

        # 
        # get fname dict for the list of ms data for the input target and config
        fname_dict = self._fname_dict(target = target, config = config, product = product, all_ms_data = True, extra_ext = extra_ext)
        this_ms_filenames = fname_dict['ms_filenames']

        # 
        # get line tag for product
        line_tag = self._kh.get_line_tag_for_line_product(product)

        # 
        # loop ms data for the input target and config
        all_chanwidths = []
        logger.debug('Current directory: "'+os.getcwd()+'"')
        for i in range(len(this_ms_filenames)):
            this_ms_filename = this_ms_filenames[i]
            logger.debug('Computing channel width in "'+this_ms_filename+'" for target '+target+', config '+config+', product '+product+', vsys '+str(vsys)+', vwidth '+str(vwidth))
            if not self._dry_run:
                this_chanwidth = cvr.compute_chanwidth_for_line(in_file = this_ms_filename, 
                                                                line = line_tag, 
                                                                vsys = vsys, 
                                                                vwidth = vwidth, 
                                                                )
                if this_chanwidth is None:
                    this_chanwidth = np.nan
                if np.isscalar(this_chanwidth):
                    this_chanwidth = [this_chanwidth]
                else:
                    this_chanwidth = this_chanwidth.flatten()
                all_chanwidths.extend(this_chanwidth)
                logger.debug('Computed channel width '+str(this_chanwidth)+' km/s in "'+this_ms_filename+'" for target '+target+', config '+config+', product '+product+', vsys '+str(vsys)+', vwidth '+str(vwidth))
        # 
        # take the coarsest chanwidth as the common_chanwidth
        if not self._dry_run:
            common_chanwidth = np.nanmax(all_chanwidths)
        else:
            common_chanwidth = 5.0 #<TODO><DEBUG>#
        logger.debug('Common channel width '+str(common_chanwidth)+' km/s for target '+target+', config '+config+', product '+product+', vsys '+str(vsys)+', vwidth '+str(vwidth))
        # 
        logger.info('END: Computing common channel width among all ms data for target '+target+', config '+config+' and product '+product+'.')
        # 
        return common_chanwidth


    def task_extract_line(
            self, 
            target = None, 
            product = None, 
            config = None, 
            extra_ext = '', 
            do_statwt = True, 
            edge_for_statwt = -1, 
            method_for_channel_regridding = 1, 
            overwrite = False, 
            ):
        """
        Extract spectral line data from ms data for the input target, config and product. 
        
        To extract the spectral line data, a common channel width needs to be calculated. 
        """
        # 
        logger.info('START: Extracting spectral line '+product+' for target '+target+' and config '+config+'.')
        # 
        # get user target channel width
        target_chanwidth = self._kh.get_channel_width_for_line_product(product=product)
        # 
        # compute the common channel width for all ms data for the input target, config and product (product is used to select spw in each ms data)
        common_chanwidth = self.task_compute_common_channel_width(target=target, product=product, config=config, extra_ext=extra_ext)
        # 
        # compute the rebinning factor to make the rebinned chanwidth as close to the target channel width as possible (but not exceeding it)
        one_plus_eps = 1.0+1e-3 #<TODO># documentation
        interpolate_cw = common_chanwidth * one_plus_eps #<TODO># documentation
        rat = target_chanwidth / interpolate_cw
        rebin_fac = int(np.round(rat))
        if rebin_fac < 1:
            rebin_fac = 1
        # 
        logger.info('Computed common channel width '+str(common_chanwidth)+' km/s for target '+target+', config '+config+' and product '+product+'.')
        logger.info('Will rebin by a factor of '+str(rebin_fac)+' to a final channel width of '+str(rebin_fac*interpolate_cw)+' km/s, given the target channel width of '+str(target_chanwidth)+' km/s as defined in "config_definitions.txt".')
        # 
        # get imaging dir and change directory
        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        # 
        # get target vsys and vwidth
        vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target)
        # 
        # get fname dict for the list of ms data for the input target and config
        fname_dict = self._fname_dict(target = target, config = config, product = product, all_ms_data = True, extra_ext = extra_ext)
        this_ms_filenames = fname_dict['ms_filenames']
        # 
        # set edge_for_statwt <TODO>
        edge_for_statwt = 25
        if product == 'co21': edge_for_statwt = 25
        if product == '13co21': edge_for_statwt = 25
        if product == 'c18o21': edge_for_statwt = 20
        # 
        # get line tag for product
        line_tag = self._kh.get_line_tag_for_line_product(product)
        # 
        # extract line for each ms data
        logger.debug('Current directory: "'+os.getcwd()+'"')
        for i in range(len(this_ms_filenames)):
            this_ms_filename = this_ms_filenames[i]
            this_ms_extracted = fname_dict['ms_extracted'][i]
            logger.debug('Extracting spectral line: "'+this_ms_extracted+'" <-- "'+this_ms_filename+'"')
            if not self._dry_run:
                # 
                # extract line channels data for each ms data
                # this includes regridding and rebinning the velocity axis
                if method_for_channel_regridding == 1:
                    # first regridding then rebinning
                    cvr.extract_line(in_file = this_ms_filename, 
                                     out_file = this_ms_extracted, 
                                     line = line_tag, 
                                     vsys = vsys, 
                                     vwidth = vwidth, 
                                     chan_fine = interpolate_cw, 
                                     rebin_factor = rebin_fac, 
                                     do_regrid_only = False, 
                                     do_regrid_first = True, 
                                     do_statwt = do_statwt, 
                                     edge_for_statwt = edge_for_statwt, 
                                     overwrite = overwrite,
                                    )
                # 
                elif method_for_channel_regridding == 2:
                    # first rebinning then regridding
                    cvr.extract_line(in_file = this_ms_filename, 
                                     out_file = this_ms_extracted, 
                                     line = line_tag, 
                                     vsys = vsys, 
                                     vwidth = vwidth, 
                                     chan_fine = target_chanwidth, 
                                     rebin_factor = rebin_fac, 
                                     do_regrid_only = False, 
                                     do_regrid_first = False, 
                                     do_statwt = do_statwt, 
                                     edge_for_statwt = edge_for_statwt, 
                                     overwrite = overwrite,
                                    )
                # 
                elif method_for_channel_regridding == 3:
                    # do only one regridding step
                    # this is not recommended because this creates non-uniform noise due to non-integer binning, 
                    # i.e., the spectral sawtooth issue
                    cvr.extract_line(in_file = this_ms_filename, 
                                     out_file = this_ms_extracted, 
                                     line = line_tag, 
                                     vsys = vsys, 
                                     vwidth = vwidth, 
                                     chan_fine = target_chanwidth, 
                                     rebin_factor = 0, 
                                     do_regrid_only = True, 
                                     do_regrid_first = False, 
                                     do_statwt = do_statwt, 
                                     edge_for_statwt = edge_for_statwt, 
                                     overwrite = overwrite,
                                    )
                else:
                    logger.error('Wrong value for method_for_channel_regridding! It is '+str(method_for_channel_regridding)+' but only 1, 2, or 3 is accepted. 1 is the default.')
                    raise ValueError('Wrong value for method_for_channel_regridding!')
        # 
        logger.info('END: Extracting spectral line '+product+' for target '+target+' and config '+config+'.')


    def task_extract_continuum(
            self, 
            target = None, 
            product = None, 
            config = None, 
            extra_ext = '', 
            do_statwt = True, 
            do_collapse = True, 
            overwrite = False, 
            ):
        """
        Extract continuum data from ms data for the input target, config and product. 
        
        """
        # 
        logger.info('START: Extracting continuum '+product+' for target '+target+' and config '+config+'.')
        # 
        # get imaging dir and change directory
        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        # 
        # get target vsys and vwidth
        vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target)
        # 
        # get lines to flag as defined in keys
        lines_to_flag = self._kh.get_lines_to_flag_for_continuum_product(product=product)
        # 
        # get fname dict for the list of ms data for the input target and config
        fname_dict = self._fname_dict(target = target, config = config, product = product, all_ms_data = True, extra_ext = extra_ext)
        this_ms_filenames = fname_dict['ms_filenames']
        # 
        # extract continuum for each ms data
        logger.debug('Current directory: "'+os.getcwd()+'"')
        for i in range(len(this_ms_filenames)):
            this_ms_filename = this_ms_filenames[i]
            this_ms_extracted = fname_dict['ms_extracted'][i]
            logger.debug('Extracting continuum: "'+this_ms_extracted+'" <-- "'+this_ms_filename+'"')
            if not self._dry_run:
                # extract_continuum
                cvr.extract_continuum(in_file = this_ms_filename, 
                                      out_file = this_ms_extracted, 
                                      lines_to_flag = lines_to_flag, 
                                      vsys = vsys, 
                                      vwidth = vwidth, 
                                      do_statwt = do_statwt, 
                                      do_collapse = do_collapse, 
                                      overwrite = overwrite, 
                                     )
        # 
        # 
        logger.info('END: Extracting continuum '+product+' for target '+target+' and config '+config+'.')
    
    ###################################
    # Recipes - combinations of tasks #
    ###################################

    ########################################
    # 20200226: DEPRECATED FUNCTIONS BELOW #
    ########################################
    
    
    
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
    
    








