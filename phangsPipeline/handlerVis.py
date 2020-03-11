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

# Spectral lines
import utilsLines as lines

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

        
#region Loops
        
    ######################################
    # Loop through all steps and targets #
    ######################################
    
    def loop_stage_uvdata(
        self,
        do_copy = True,
        do_concat = True,
        do_remove_staging = False,
        do_custom = False,
        do_contsub = False, 
        do_extract_line = True,
        do_extract_cont = True,
        extra_ext = '',       
        make_directories = True,
        statwt_cont = True,
        collapse_cont = True,
        timebin = '0s',
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

                        self.task_split(
                            target = this_target, 
                            project = this_project, 
                            array_tag = this_array_tag, 
                            obsnum = this_obsnum, 
                            product = this_product,
                            timebin = timebin,
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
        
        # Clean up the staged measurement sets. They cost time to
        # re-split but have a huge disk imprint and are totally
        # redundant with the concatenated data and original
        # data. Saves 33% on disk space.

        for this_target, this_project, this_array_tag, this_obsnum in \
                self._kh.loop_over_input_ms(target=target_list,
                                            config=config_list,
                                            project=just_projects):

                for this_product in product_list:
                    
                    if do_remove_staging:

                        self.remove_staged_products(
                            target = this_target, 
                            project = this_project, 
                            array_tag = this_array_tag, 
                            obsnum = this_obsnum, 
                            product = this_product,
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
                            config = this_config, 
                            product = this_product, 
                            exact = False,
                            extra_ext_in = "noregrid",
                            contsub = "prefer",
                            # could add algorithm flags here                            
                            overwrite = overwrite, 
                            )

                if this_product in self._kh.get_continuum_products():
                    if do_extract_cont:

                        self.task_extract_continuum(
                            target = this_target, 
                            product = this_project, 
                            config = this_config, 
                            extra_ext_in = "noregrid", 
                            do_statwt = statwt_cont, 
                            do_collapse = collapse_cont, 
                            overwrite = overwrite, 
                            )                
                
        return()

#endregion

#region Tasks
        
    ##########################################
    # Tasks - individual operations on data. #
    ##########################################
    
    def task_split(
            self, 
            target = None,
            project = None,
            array_tag = None,
            obsnum = None,
            product = None,
            extra_ext_out = '',
            do_statwt = False, 
            timebin = '0s',
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
            logger.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if project is None:
            logger.error("Please specify a project.")
            raise Exception("Please specify a project.")

        if array_tag is None:
            logger.error("Please specify an array_tag.")
            raise Exception("Please specify an array_tag.")

        if obsnum is None:
            logger.error("Please specify an obsnum.")
            raise Exception("Please specify an obsnum.")
        
        infile = self._kh.get_file_for_input_ms(
            target=target, project=project, array_tag=array_tag, obsnum=obsnum)

        if infile is None:
            logger.error("Infile not found. Returning.")
            return()

        if not os.path.isdir(infile):
            logger.error("Infile not found. Returning.")
            return()
        
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
                max_chanwidth_kms = self._kh.get_channel_width_for_line_product(product)

                if not self._dry_run and casa_enabled:
                    spw = cvr.find_spws_for_line(infile = infile, 
                                                 line = this_line, 
                                                 max_chanwidth_kms = max_chanwidth_kms,
                                                 vsys_kms = vsys, vwidth_kms = vwidth)
                    if spw is None:
                        logger.warning("No SPWs meet the selection criteria. Skipping.")
                        return()

        # TBD - Work out time binning using some logic once we
        # generate some context from the keyhandler. For now, just
        # pass along the parameter from the user.

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
                field = field,
                spw = spw,             
                timebin = timebin,   
                do_statwt = do_statwt, 
                overwrite = overwrite, 
                )

        return()

    def remove_staged_products(
        self,
        target = None,
        project = None,
        array_tag = None,
        obsnum = None,
        product = None,
        extra_ext = '',
        ):
        """
        Remove 'staged' visibility products, which are intermediate
        between the calibrated data and the concated measurement sets
        that we begin processing on. Run this step after concat to
        reduct the disk footprint of the pipeline.
        """
        
        if target is None:
            logger.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if project is None:
            logger.error("Please specify a project.")
            raise Exception("Please specify a project.")

        if array_tag is None:
            logger.error("Please specify an array_tag.")
            raise Exception("Please specify an array_tag.")

        if obsnum is None:
            logger.error("Please specify an obsnum.")
            raise Exception("Please specify an obsnum.")
        
        infile = fnames.get_staged_msname(
            target=target, project=project, array_tag=array_tag, obsnum=obsnum, 
            product=product, ext=extra_ext)
            
        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Clearing intermediate staged u-v data for "+infile)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
       
        # Change to the imaging directory for the target

        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)

        if not self._dry_run:
            
            os.system('rm -rf '+infile)

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
            logger.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if product is None:
            logger.error("Please specify a product.")
            raise Exception("Please specify a product.")

        if config is None:
            logger.error("Please specify a config.")
            raise Exception("Please specify a config.")

        # Change to the imaging directory for the target
        
        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        
        # Generate the list of staged measurement sets to combine
        
        staged_ms_list = []        
        for this_target, this_project, this_array_tag, this_obsnum in \
                self._kh.loop_over_input_ms(target=target, config=config,
                                            project=just_projects):
                
                this_staged_ms = fnames.get_staged_msname(
                    target=this_target, project=this_project, array_tag=this_array_tag, 
                    obsnum=this_obsnum, product=product, ext=extra_ext_in)
                if os.path.isdir(this_staged_ms):
                    staged_ms_list.append(this_staged_ms)
                else:
                    logger.warning("MS not found and will be dropped from concat: "+str(this_staged_ms))
                    logger.warning("This might or might not be a problem.")


        if len(staged_ms_list) == 0:
            logger.warning("No measurement sets to concatenate, returning.")
            return()

        # Generate the outfile name

        outfile = fnames.get_vis_filename(
            target=target, config=config, product=product, 
            ext=extra_ext_out, suffix=None)
                
        # Revise here

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Concatenating staged and split u-v data for:")
        logger.info("... target: "+target)
        logger.info("... product: "+product)
        logger.info("... config: "+config)
        logger.info("... files: "+str(staged_ms_list))
        logger.info("... output: "+str(outfile))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
                   
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
            logger.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if product is None:
            logger.error("Please specify a product.")
            raise Exception("Please specify a product.")

        if config is None:
            logger.error("Please specify an config.")
            raise Exception("Please specify an config.")
        
        infile = fnames.get_vis_filename(
            target=target, config=config, product=product, 
            ext=extra_ext_in, suffix=None)

        # Get lines to exclude. Previously had used the continuum fit,
        # but that mixes product definitions ... this is tricky
        # because those lines are defined for the continuum. Need to
        # examine conventions, but probably make this a part of
        # product definition.
        
        # get target vsys and vwidth
        vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target)

        # lines_to_exclude = self._kh.get_lines_to_flag_for_continuum_product(product=product)        
        lines_to_exclude = self._kh.get_line_tag_for_line_product(product)

        ranges_to_exclude = lines.get_ghz_range_for_list(
            line_list=lines_to_exclude, vsys_kms=vsys, vwidth_kms=vwidth)
            
        # Need to look up fit order and other tuning parameters eventually
        fitorder = 0

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("u-v continuum subtraction for "+infile)
        logger.info("... excluding frequency ranges:"+str(ranges_to_exclude))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        # Change to the imaging directory for the target

        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)
        
        if not self._dry_run and casa_enabled:

            cvr.contsub(infile = infile, 
                        # outfile is TBD
                        ranges_to_exclude = ranges_to_exclude,              
                        overwrite = overwrite, 
                        fitorder = fitorder,
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

    def task_extract_line(
            self, 
            target = None, 
            product = None, 
            config = None, 
            exact = False,
            extra_ext_in = '',
            contsub = "prefer",
            extra_ext_out = '',
            do_statwt = True, 
            edge_for_statwt = -1,
            method = "regrid_then_rebin",
            overwrite = False, 
            ):
        """
        Extract spectral line data from ms data for the input target,
        config and product.        
        """
        
        # Error checking

        if target is None:
            logger.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if product is None:
            logger.error("Please specify a product.")
            raise Exception("Please specify a product.")

        if config is None:
            logger.error("Please specify a config.")
            raise Exception("Please specify a config.")

        # Option re: continuum subtraction

        valid_contsub_options = ['prefer','require','none']
        if contsub.lower().strip() not in valid_contsub_options:
            logger.error("Please choose a valid contsub option.")
            logger.error("Valid options are:"+str(valid_contsub_options))
            raise Exception("Please choose a valid contsub option.")

            
        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)

        if contsub == 'require':

            infile = fnames.get_vis_filename(
                target=target, config=config, product=product, 
                ext=extra_ext_in, suffix='contsub')

        if contsub == 'prefer':

            infile = fnames.get_vis_filename(
                target=target, config=config, product=product, 
                ext=extra_ext_in, suffix='contsub')

            if not os.path.isdir(infile):
                infile = fnames.get_vis_filename(
                    target=target, config=config, product=product, 
                    ext=extra_ext_in, suffix=None)

        if contsub == 'none':
            infile = fnames.get_vis_filename(
                target=target, config=config, product=product, 
                ext=extra_ext_in, suffix=None)

        if not os.path.isdir(infile):
            logger.warning("Input file not found. Returning.")
            return()
                
        outfile = fnames.get_vis_filename(
            target=target, config=config, product=product, 
            ext=extra_ext_out, suffix=None)

        # Extract the spectral information needed for the regrid

        vsys_kms, vwidth_kms = self._kh.get_system_velocity_and_velocity_width_for_target(target)
        line_to_extract = self._kh.get_line_tag_for_line_product(product)

        valid_methods = ['regrid_then_rebin','rebin_then_regrid','just_regrid','just_rebin']
        if method.lower().strip() not in valid_methods:
            logger.error("Not a valid line extraction medod - "+str(method))
            raise Exception("Please specify a valid line extraction method.")

        target_chanwidth = self._kh.get_channel_width_for_line_product(product)

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Extracting spectral product:")
        logger.info("Line: "+str(line_to_extract))
        logger.info("Vsys [km/s]: "+str(vsys_kms))
        logger.info("Vwidth [km/s]: "+str(vwidth_kms))
        logger.info("Method: "+str(method))
        logger.info("From file: "+str(infile))
        logger.info("To file: "+str(outfile))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")

        if not self._dry_run and casa_enabled:

            cvr.batch_extract_line(
                infile_list = [infile],
                outfile = outfile,
                target_chan_kms = target_chanwidth,
                line = line_to_extract,
                vsys_kms=vsys_kms, vwidth_kms=vwidth_kms,
                method = 'regrid_then_rebin',
                exact = exact,
                overwrite = True,
                )

        return()

    def task_extract_continuum(
            self, 
            target = None, 
            product = None, 
            config = None, 
            extra_ext_in = '', 
            do_statwt = True, 
            do_collapse = True, 
            overwrite = False, 
            ):
        """
        Extract continuum data from ms data for the input target, config and product.         
        """
        
        # Error checking

        if target is None:
            logger.error("Please specify a target.")
            raise Exception("Please specify a target.")

        if product is None:
            logger.error("Please specify a product.")
            raise Exception("Please specify a product.")

        if config is None:
            logger.error("Please specify a config.")
            raise Exception("Please specify a config.")

        infile = fnames.get_vis_filename(
            target=target, config=config, product=product, 
            ext=extra_ext_in, suffix=None)

        # get target vsys and vwidth
        vsys, vwidth = self._kh.get_system_velocity_and_velocity_width_for_target(target)

        # get lines to flag as defined in keys
        lines_to_flag = self._kh.get_lines_to_flag_for_continuum_product(product=product)
                
        if len(lines_to_flag) > 0:
            ranges_to_exclude = lines.get_ghz_range_for_list(
                line_list=lines_to_flag, vsys_kms=vsys, vwidth_kms=vwidth)
        else:
            ranges_to_exclude = []

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Extracting continuum from "+infile)
        logger.info("... flagging ranges: "+str(ranges_to_exclude))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
       
        this_imaging_dir = self._kh.get_imaging_dir_for_target(target, changeto=True)

        if not self._dry_run and casa_enabled:

            cvr.extract_continuum(
                infile = infile, 
                outfile = outfile, 
                ranges_to_exclude = ranges_to_exclude,
                do_statwt = do_statwt,
                do_collapse = do_collapse, 
                overwrite = overwrite, 
                )

        return()

#endregion
