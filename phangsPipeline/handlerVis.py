"""
handlerVis (VisHandler object)
    
The PHANGS pipeline to handle staging and pre-processing of uv data
before imaging. Works through a single big class (the
VisHandler). This needs to be attached to a keyHandler to access the
target, product, and configuration keys and locations of the uv data.

To run the individual routines, this code needs to be run inside
CASA. See an example application inside stage_7m_co21.py .

Example:

    $ casa
    from phangsPipeline import handlerKeys as kh
    from phangsPipeline import uvdataHandler as uvh
    this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
    this_uvh = uvh.VisHandler(key_handler = this_kh, dry_run = False)
    # Set which data to process
    this_uvh.set_line_products(only=['co21'])
    this_uvh.set_interf_configs(only=['12m+7m'])
    this_uvh.set_targets(only=['ngc3621'])
    # Process
    this_uvh.loop_stage_uvdata()

"""

# 20200303: "get_ms_filenames_and_filepaths" --> "get_all_input_ms"

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


class VisHandler(handlerTemplate.HandlerTemplate):
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
        
    ######################################
    # Loop through all steps and targets #
    ######################################
    
    def loop_stage_uvdata(
        self,
        do_copy = True,
        do_concat = True,
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
        product_list = self.get_all_products()
        config_list = self.get_interf_configs()

        # Our first step uses CASA's split to extract the relevant
        # fields and spectral windows from each input data set.

        for this_target, this_project, this_array_tag, this_obsnum in \
                self._kh.loop_over_input_ms(target=target_list,
                                            config=config_list,
                                            project=just_projects):

                for this_product in product_list:
                    
                    if do_copy:

                        self.task_copy_and_split(
                            target = this_target, 
                            project = this_project, 
                            array_tag = this_array_tag, 
                            obsnum = this_obsnum, 
                            product = this_product,
                            # could add algorithm flags here
                            overwrite = overwrite, 
                            )

        # Next we concatenate linked data sets (sharing a
        # target+config+product ) into unified measurement sets.

        for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_line=True,just_interf=True):

                if do_concat:

                    self.task_concat_uvdata(
                        target = this_target, 
                        config = this_config, 
                        product = this_product, 
                        extra_ext_out = "noregrid",                         
                        overwrite = overwrite, 
                        )
        
        # Next we apply uv continuum subtraction. We may later offer
        # an algorithm choice to do this before or after
        # regridding. The correct choice depends on some details of
        # the observation setup. This one assumes one SPW maps to one
        # line and so may be more common.

        for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_line=True,just_interf=True):

                if do_contsub:
                    
                    self.task_contsub(
                        target = this_target, 
                        config = this_config, 
                        product = this_product, 
                        extra_ext_in = "noregrid", 
                        # could add algorithm flags here
                        overwrite = overwrite, 
                        )
        
        # Now we reprocess the data to have the desired spectral
        # setup. This involves rebinning and regridding for line
        # products and flagging and integration for continuum
        # products.

        for this_target, this_product, this_config in \
                self.looper(do_targets=True,do_products=True,do_configs=True,
                            just_line=True,just_interf=True):

                if this_product in self._kh.get_line_products():
                    if do_extract_line:

                        self.task_extract_line(
                            target = this_target, 
                            project = this_project, 
                            array_tag = this_array_tag, 
                            obsnum = this_obsnum, 
                            product = this_product, 
                            extra_ext_in = "noregrid",
                            # could add algorithm flags here
                            overwrite = overwrite, 
                            )

                if this_product in self._kh.get_continuum_products():
                    if do_extract_cont:

                        self.task_extract_continuum(
                            target = this_target, 
                            project = this_project, 
                            array_tag = this_array_tag, 
                            obsnum = this_obsnum, 
                            product = this_product, 
                            extra_ext_in = "noregrid",
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

            if product in self._kh.get_line_products():

                this_line = self._kh.get_line_tag_for_line_product(product)
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
            just_projects = None,
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
        
        # Generate the list of staged measurement sets to combine
        
        staged_ms_list = []        
        for this_target, this_project, this_array_tag, this_obsnum in \
                self._kh.loop_over_input_ms(target=target, config=config,
                                            project=just_projects):
                
                this_staged_ms = fnames.get_staged_msname(
                    target=this_target, project=this_project, array_tag=this_array_tag, 
                    this_obsnum=obsnum, product=product, ext=extra_ext_in)
                staged_ms_list.append(this_staged_ms)

        # Generate the outfile name

        outfile = fnames.get_vis_filename(
            target=target, config=config, product=product, 
            ext=extra_ext_out, suffix=None)
                
        # Revise here

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Splitting (by spw) then concatenating u-v data for:")
        logger.info("... target: "+target)
        logger.info("... product: "+product)
        logger.info("... config: "+config)
        logger.info("... files: "+str(staged_ms_list))
        logger.info("... output: "+str(outfile))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        # Change to the imaging directory for the target
        
        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        
        # Concatenate the measurement sets

        if not self._dry_run and casa_enabled:         
            
            cvr.concat_ms(infile_list = staged_ms_list, 
                          outfile = outfile, 
                          overwrite = overwrite, 
                          )

        return()
    
    def task_contsub(
            self, 
            target = None, 
            config = None, 
            product = None, 
            extra_ext_in = '', 
            overwrite = False, 
            ):
        """
        Run u-v plane continuum subtraction on an individual input
        measurement set.
        """

        if target is None:
            logging.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if product is None:
            logging.error("Please specify a product.")
            raise Exception("Please specify a product.")

        if config is None:
            logging.error("Please specify an config.")
            raise Exception("Please specify an config.")
        
        infile = fnames.get_vis_filename(
            target=target, config=config, product=product, 
            ext=extra_ext_in, suffix=None)
        
        # get target vsys and vwidth
        vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target)

        # get lines to exclude from the continuum fit ... this is
        # tricky because those lines are defined for the
        # continuum. Need to examine conventions

        # lines_to_exclude = self._kh.get_lines_to_flag_for_continuum_product(product=product)
        
        lines_to_exclude = self._kh.get_line_tag_for_line_product(product)

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("u-v continuum subtraction for "+infile)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        # Change to the imaging directory for the target

        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        
        if not self._dry_run and casa_enabled:

            cvr.contsub(infile = infile, 
                        lines_to_flag = lines_to_exclude, 
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
