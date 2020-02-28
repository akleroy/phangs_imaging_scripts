"""imagingHandler

This module makes images and cubes out of the uv data products created by uvdataHandler. 

This code needs to be run inside CASA. 

Example:
    $ casa
    from phangsPipeline import keyHandler as kh
    from phangsPipeline import imagingHandler as imh
    this_kh = kh.KeyHandler(master_key = 'config_keys/master_key.txt')
    this_imh = uvh.ImagingHandler(key_handler = this_kh)
    this_imh.set_targets(only = ['ngc3627'])
    this_imh.loop_imaging()

Args:
    make_dirty_image (boolean): If `True` then make dirty image cube from ms data. Default is `True`.
    revert_to_dirty (boolean): If `True` then reset current image cube by the dirty image cube (mainly for debugging). Default is `False`.
    read_in_clean_mask (boolean): If `True` then read in clean mask as defined in "cleanmask_keys.txt". Default is `True`.
    run_multiscale_clean (boolean): If `True` then run multiscale clean. Default is `True`.
    revert_to_multiscale (boolean): If `True` then reset current image cube by the multiscale cleaned image cube (mainly for debugging). Default is `False`.
    make_singlescale_mask (boolean): If `True` then make singlescale clean mask based on current image cube. Default is `True`.
    run_singlescale_clean (boolean): If `True` then run singlescale clean. Default is `True`.
    do_export_to_fits (boolean): If `True` then export ms data folders into fits-format image cube files. Default is `True`.

Notes:
    This code calls following functions:
        import casaImagingRoutines as imr
        import casaMaskingRoutines as msr
        if make_dirty_image: imr.make_dirty_map             
        if revert_to_dirty: imr.replace_cube_with_copy      
        if read_in_clean_mask: msr.import_and_align_mask    
        if run_multiscale_clean: imr.multiscale_loop        
        if revert_to_multiscale: imr.replace_cube_with_copy 
        if make_singlescale_mask: msr.signal_mask           
        if run_singlescale_clean: imr.singlescale_loop      
        if do_export_to_fits: imr.export_to_fits            

"""

#<DONE># 20200210 dzliu: self._kh._cleanmask_dict is always None. It is not yet implemented in "keyHandler.py"!
#<TODO># 20200214 dzliu: will users want to do imaging for individual project instead of concatenated ms?
#<TODO># 20200214 dzliu: needing KeyHandler API: 
#<DONE># 20200214 dzliu:     self._kh._cleanmask_dict   --> self._kh.get_cleanmask() # input target name, output clean mask file
#<TODO># 20200214 dzliu:     self._kh._target_dict      --> self._kh.get_target_dict() # for rastring, decstring
#<TODO># 20200214 dzliu:     self._kh._override_dict    --> self._kh.get_overrides()
#<TODO># 20200214 dzliu:     self._kh._dir_keys         --> self._kh.get_target_name_for_multipart_name()
#<TODO># 20200214 dzliu:     self._kh._config_dict['interf_config']['clean_scales_arcsec'] # angular scales
#<TODO># 20200218 dzliu: CASA 5.4.0 works, but CASA 5.6.0 does not work!!
#<TODO># 20200218 dzliu: revert does not work!
#<TODO># 20200218 dzliu: need to test 'cont'

import os, sys, re, shutil
import glob
import numpy as np

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

casa_enabled = (sys.argv[0].endswith('start_casa.py')) #<TODO># check whether we are inside CASA environment

if casa_enabled:
    logger.debug('casa_enabled = True')
    import casaImagingRoutines as imr
    import casaMaskingRoutines as msr
    reload(imr) #<TODO><DEBUG># 
    reload(msr) #<TODO><DEBUG># #<TODO># we still need imr.cleanCall for dry_run?
else:
    logger.debug('casa_enabled = False')
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from clean_call import CleanCall

import utils
import line_list
import handlerTemplate


class ImagingHandler(handlerTemplate.HandlerTemplate):
    """
    Class to makes image cubes out of uv data from each spectral line and continuum of each galaxy. 
    """
    
    ############
    # __init__ #
    ############
    
    def __init__(
        self, 
        key_handler = None, 
        dry_run = False, 
        ):
        
        # inherit template class
        handlerTemplate.HandlerTemplate.__init__(self,key_handler = key_handler, dry_run = dry_run)
        
    ################
    # loop_imaging #
    ################
    
    def loop_imaging(
        self, 
        do_dirty_image = True,
        do_revert_to_dirty = True,
        do_read_clean_mask = True, 
        do_multiscale_clean = True
        do_revert_to_multiscale = True,
        do_singlscale_mask = True,
        do_singlescale_clean = True        
        do_export_to_fits = True, 
        extra_ext_in = None,
        suffix_in = None,
        extra_ext_out = None,        
        recipe = 'phangsalma',
        make_directories = True, 
        dynamic_sizing = True,
        force_square = False,
        overwrite = False, 
        ): 
        """
        Loops over the full set of targets, products, and
        configurations to do the imaging. Toggle the parts of the loop
        using the do_XXX booleans. Other choices affect algorithms
        used.
        """

        if len(self.get_targets()) == 0:            
            logger.error("Need a target list.")
            return(None)
        
        if len(self.get_all_products()) == 0:            
            logger.error("Need a products list.")
            return(None)
        
        known_recipes = ['phangsalma']
        if recipe not in known_recipes:
            logger.error("Recipe not known ", recipe)
            return(None)

        if make_directories:
            self._kh.make_missing_directories(imaging=True)
        
        # Loop over targets, products, and interferometric configs

        for this_target, this_product, this_config in \
                self.looper(do_targets=True, do_products=True, do_configs=True, just_interf=True):

            # Change to the relevant directory

            this_imaging_dir = self._kh.get_imaging_dir_for_target(this_target, changeto=True)

            # print starting message
            logger.info("")
            logger.info("--------------------------------------------------------")
            logger.info("START: Imaging ", target, config, product)
            logger.info("--------------------------------------------------------")
            logger.info('Imaging recipe: '+recipe)

            if recipe == 'phangsalma':
                self.recipe_phangsalma_imaging(
                    target=this_target,
                    product=this_product,
                    config=this_config,
                    extra_ext_in = extra_ext_in,
                    suffix_in = suffix_in,
                    extra_ext_out = extra_ext_out,
                    do_dirty_image = True,
                    do_revert_to_dirty = True,
                    do_read_clean_mask = True, 
                    do_multiscale_clean = True
                    do_revert_to_multiscale = True,
                    do_singlscale_mask = True,
                    do_singlescale_clean = True        
                    do_export_to_fits = True, 
                    dynamic_sizing = dynamic_sizing,
                    force_square = force_square, 
                    overwrite = overwrite, 
                    )
                
            # print ending message
            logger.info("--------------------------------------------------------")
            logger.info("END: Imaging ", target, config, product)
            logger.info("--------------------------------------------------------")
        # 
        # end of for looper
    # 
    # end of loop_imaging()

    ####################
    # File name lookup #
    ####################
    
    ##################################################################
    # Tasks - discrete steps on target, product, config combinations #
    ##################################################################
    
    def task_initialize_clean_call(
        self,
        target = None, 
        config = None, 
        product = None, 
        extra_ext_in = None,
        suffix_in = None,
        extra_ext_out = None,
        stage = 'dirty',
        ):
        """
        Initialize a clean call object for a target, config, product
        combination and an imaging stage.
        """

        # check user input. 
        if target is None or config is None or product is None:
            logger.error('Please input target, config and product!')
            logger.error('e.g., target = "ngc3627_1", config = "12m+7m", product = "co21"')
            raise Exception('Please input target, config and product!')

        # Look up the recipes for this case
        recipe_list = self._kh.get_imaging_recipes(config=config, product=product)

        # Initialize the clean call with the appropriate recipe list
        clean_call = CleanCall(recipe_list)

        # Get the visibility name
        vis_file = self._kh.get_vis_filename(
            target=target, product=product, config=config,
            ext=extra_ext_in, suffix_in=None)

        # Test existence
        full_vis_file = self._kh.get_imaging_dir_for_target(target=target)+vis_file
        if not os.pwd.isdir(full_vis_file):
            logger.error('Visibility file not found: ', vis_file)
            raise Exception('Missing visibility needed for imaging.')            

        # Set the visibility file (note that we assume we are in the working directory)

        clean_call.set_param('vis', vis_file, nowarning=True)

        # Set the output image file name (note no suffix for imaging root)
         
        image_root = self._kh.get_cube_filename(
            target=target, product=product, config=config,
            ext=extra_ext_out, casa=True, casaext='')

        clean_call.set_param('imagename', image_root, nowarning=True)

        # Get the phase center associated with the target
        rastring, decstring = self._kh.get_phasecenter_for_target(target=target)
        phasecenter = 'J2000 '+rastring+' '+decstring

        clean_call.set_param('phasecenter', phasecenter)

        return(clean_call)

    def task_pick_cell_and_imsize(
        self, 
        clean_call = None,
        oversamp=5,
        force_square=False,
        check_files=True,
        ):
        """
        Pick a cell size and imsize for a clean call by inspecting the
        visibility data, then attach these values to the clean call.
        
        This will call casaImagingRoutines.estimate_cell_and_imsize()
        first, then apply custom overrides if they exist.
        """
                
        if clean_call is None:
            logger.warning("Require a clean_call object. Returning.")
            return()

        # Look up the file name
        vis = clean_call.get_param('vis')

        # Check visibility existence
        if check_files:
            if not (os.path.isdir(vis_file)):
                logger.warning("Missing visibility data "+vis_file)
                return()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Calculating target cell size and image size for:")
        logger.info(vis_file)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
        
        # Call the estimating routine
        if (not self._dry_run) and casa_enabled:
            cell, imsize = \
                imr.estimate_cell_and_imsize(vis_file, 
                                             oversamp,
                                             force_square=force_square)
        else:
            cell, imsize = '0.1arcsec', [1000,1000]
            logger.info('DRY RUN skips calling imr.estimate_cell_and_imsize()')
          
        logger.info('cell='+cell_size_str+'arcsec, imsize=['+x_size_string+','+y_size_string+']')

        clean_call.set_param('cell',cell,nowarning=True)
        clean_call.set_param('imsize',imsize,nowarning=True)

        # Return

        return(cell, imsize)
    
    def task_assign_multiscales(
        self, 
        clean_call = None,
        config = None, 
        ):
        """
        Look up the appropriate scales for multi-scale clean
        associated with this config.
        """
        
        if clean_call is None:
            logger.warning("Require a clean_call object. Returning.")
            return()

        scales_to_clean = self._kh.get_clean_scales_for_config(config=config)
        clean_call.set_multiscale_arcsec(scales=scales_to_clean)

        return()

    def task_make_dirty_image(
        clean_call = None,
        backup = True,
        ):
        """
        Execute the attached clean call with zero iterations to create
        a dirty image.
        """

        if clean_call is None:
            logger.warning("Require a clean_call object. Returning.")
            return()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Making dirty image:")
        logger.info(str(clean_call.get_param('imagename')))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        if (not self._dry_run) and casa_enabled:
            imr.make_dirty_image(clean_call)
            if backup:
                imr.copy_imaging(
                    input_root=clean_call.get_param('imagename'),
                    output_root=clean_call.get_param('imagename')+'_dirty',
                    wipe_first=True)
                
        return()

    def task_revert_to_imaging(
        clean_call = None,
        tag = 'dirty'
        ):
        """
        Reset the current imaging stack to an earlier imaging
        state. Requires a tag that will be used to define the old
        imaging moved into place. Wipes the current imaging.
        """

        if clean_call is None:
            logger.warning("Require a clean_call object. Returning.")
            return()

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Resetting to "+tag+" imaging:")
        logger.info(str(clean_call.get_param('imagename')))
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
            
        if (not self._dry_run) and casa_enabled:
            imr.copy_imaging(
                input_root=clean_call.get_param('imagename')+'_'+tag,
                output_root=clean_call.get_param('imagename'),
                wipe_first=True)

        return()

    def task_read_clean_mask(
        clean_call = None,
        target = None,
        product = None,
        config = None,
        ):
        """
        """

        if read_in_clean_mask:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("READING IN THE CLEAN MASK.")
                logger.info("")
                if clean_call.clean_mask_file is not None:
                    msr.import_and_align_mask(
                        in_file=clean_call.clean_mask_file,
                        out_file=clean_call.image_root+'.mask',
                        template=clean_call.image_root+'.image',
                        )
                    clean_call.usemask = 'user'
                else:
                    logger.info("No clean mask defined.")
                    clean_call.usemask = 'pb'
            else:
                logger.info("DRY RUN skips setting the clean mask")
                clean_call.usemask = 'pb'
        pass

    def task_multiscale_clean(
        clean_call = None,
        ):
        """
        """
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("RUNNING THE MULTISCALE CLEAN.")
                logger.info("")
                imr.multiscale_loop(
                    clean_call = clean_call,
                    record_file = clean_call.image_root+'_multiscale_record.txt',
                    delta_flux_threshold=0.01,
                    absolute_threshold=None,
                    snr_threshold=4.0,
                    stop_at_negative=True,
                    max_loop = 20
                    )
            else:
                logger.info("DRY RUN skips running imr.multiscale_loop() and outputting to "+clean_call.image_root+'_multiscale'+'.*')

 multiscale_loop(
    clean_call = None,
    record_file=None,
    delta_flux_threshold=0.02,
    absolute_delta=True,
    absolute_threshold=None,
    snr_threshold=4.0,
    stop_at_negative=True,
    max_loop = 20
    ):
    """
    Carry out an iterative multiscale clean loop.
    """
    
    # Check that we have a vile clean call

    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")
    
    # Figure out the scales to use in pixel units

    cell_as_num = float((clean_call.cell_size.split('arcsec'))[0])
    scales_as_pix = []
    for scale in clean_call.scales_as_angle:
        scales_as_pix.append(int(scale/cell_as_num))
        
    clean_call.deconvolver = 'multiscale'
    clean_call.scales_as_pix = scales_as_pix
    clean_call.calcres = False
    clean_call.calcpsf = False

    logger.info("I will use the following scales: ")
    logger.info("... as pixels: " + str(clean_call.scales_as_pix))
    logger.info("... as arcseconds: " + str(clean_call.scales_as_angle))

    # Call the loop

    clean_loop(
        clean_call=clean_call,
        record_file=record_file,
        delta_flux_threshold=delta_flux_threshold,
        absolute_delta=absolute_delta,
        absolute_threshold=absolute_threshold,
        snr_threshold=snr_threshold,
        stop_at_negative=stop_at_negative,
        max_loop = max_loop, 
        log_ext = 'multiscale', #<TODO># set log_ext to output log files
        )

        pass

    def task_singlescale_mask(
        clean_call = None
        ):
        """
        """
        pass

        if make_singlescale_mask:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("MAKING THE MASK FOR SINGLE SCALE CLEAN.")
                logger.info("")
                msr.signal_mask(
                    cube_root=clean_call.image_root,
                    out_file=clean_call.image_root+'.mask',
                    operation='AND',
                    high_snr=4.0,
                    low_snr=2.0,
                    absolute=False)
                clean_call.usemask='user'
            else:
                logger.info("DRY RUN skips running msr.signal_mask() and setting the clean mask")
   
    def task_singlescale_clean(        
        clean_call = None,
        ):
        """
        """
        pass
    
        if run_singlescale_clean:
            if (not self._dry_run) and casa_enabled:
                logger.info("")
                logger.info("RUNNING THE SINGLE SCALE CLEAN.")
                logger.info("")
                imr.singlescale_loop(
                    clean_call = clean_call,
                    record_file = clean_call.image_root+'_singlescale_record.txt',
                    delta_flux_threshold=0.01,
                    absolute_delta=True,
                    absolute_threshold=None,
                    snr_threshold=1.0,
                    stop_at_negative=False,
                    remask=False,
                    max_loop = 20
                    )
            else:
                logger.info("DRY RUN skips running imr.singlescale_loop() and outputting to "+clean_call.image_root+'_singlescale'+'.*')

                singlescale_loop(
                    clean_call = None,
                    scales_as_angle=[],
                    record_file=None,
                    delta_flux_threshold=0.02,
                    absolute_delta=True,
                    absolute_threshold=None,
                    snr_threshold=4.0,
                    stop_at_negative=True,
                    remask=False,
                    max_loop = 20
                    ):
    """
    Carry out an iterative multiscale clean loop.
    """
    
    # Check that we have a vile clean call

    if not isinstance(clean_call, CleanCall):
        logger.error("Please input a valid clean call!")
        raise Exception("Please input a valid clean call!")
        
    clean_call.deconvolver = 'hogbom'
    clean_call.calcres = False
    clean_call.calcpsf = False

    # Call the loop

    clean_loop(
        clean_call=clean_call,
        record_file=record_file,
        delta_flux_threshold=delta_flux_threshold,
        absolute_delta=absolute_delta,
        absolute_threshold=absolute_threshold,
        snr_threshold=snr_threshold,
        stop_at_negative=stop_at_negative,
        remask=remask,
        max_loop = max_loop, 
        log_ext = 'singlescale', #<TODO># set log_ext to output log files
        )


    def task_export_to_fits(
        clean_call = None,
        tag = None,
        ):
        """ 
        Export the results of a clean_call to FITS files. Optionally
        append _tag to the file.
        """

        if clean_call is None:
            logger.warning("Require a clean_call object. Returning.")
            return()

        image_root = str(clean_call.get_param('imagename'))
        if tag is not None:
            image_root += '_'+tag

        logger.info("")
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("Exporting to FITS:"+image_root)
        logger.info("&%&%&%&%&%&%&%&%&%&%&%&%&%")
        logger.info("")
                    
        if (not self._dry_run) and casa_enabled:
            imr.export_to_fits(image_root)
        
        return()

    ###########################
    # recipe_build_clean_call #
    ###########################
    
    def recipe_build_clean_call(
        self, 
        target = None, 
        config = None, 
        product = None, 
        tag = '', 
        force_square = False, 
        overwrite = False, 
        dry_run = False, 
        ):
        """Build a clean call before running the clean task with the function phangsImagingRecipe().
        """
                

        # 
        # get cell size and imsize
        cell_size, x_size, y_size = \
            self.task_pick_cell_and_imsize(\
                clean_call.vis, 
                force_square=force_square) #<TODO># dzliu: this can be improved, and this already included checking overrides by ms file name.
        image_size = [int(x_size), int(y_size)]
        
        clean_call.cell_size = cell_size
        clean_call.image_size = image_size
        
        # Look up the line and data product
        # TODO - this should generic, just get freq for line given the line name
        
        if product == 'co21':
            clean_call.specmode = 'cube'
            clean_call.restfreq_ghz = line_list.line_list['co21']

        if product == 'c18o21':
            clean_call.specmode = 'cube'
            clean_call.restfreq_ghz = line_list.line_list['c18o21']

        # Continuum case

        if product == 'cont':
            clean_call.specmode = 'mfs'
            clean_call.restfreq_ghz = -1.0
        
        # Set angular scales to be used in multiscale clean
        #scales_for_clean = self._config_dict['interf_config'][this_config]['clean_scales_arcsec']
        
        if config == '7m':
            clean_call.pblimit = 0.25
            clean_call.smallscalebias = 0.6
            clean_call.scales_as_angle = [0, 5, 10]
        elif config == '12m':
            clean_call.smallscalebias = 0.6
            clean_call.scales_as_angle = [0, 1, 2.5, 5]
        elif config == '12m+7m':
            clean_call.smallscalebias = 0.8
            clean_call.scales_as_angle = [0, 1, 2.5, 5, 10]
        
        # Look up overrides in the imaging parameters
        #logger.debug(str(clean_call))
        #logger.debug('Checking overrides for "'+clean_call.image_root+'"')
        if self._kh.has_overrides_for_key(clean_call.image_root):
            clean_call.smallscalebias = float(self._kh.get_overrides(key = clean_call.image_root, param = 'smallscalebias', default = clean_call.smallscalebias))
            clean_call.image_size[0] = int(self._kh.get_overrides(key = clean_call.image_root, param = 'x_size', default = clean_call.image_size[0]))
            clean_call.image_size[1] = int(self._kh.get_overrides(key = clean_call.image_root, param = 'y_size', default = clean_call.image_size[1]))
            clean_call.pblimit = float(self._kh.get_overrides(key = clean_call.image_root, param = 'pblimit', default = clean_call.pblimit))
            clean_call.scales_as_angle = self._kh.get_overrides(key = clean_call.image_root, param = 'scales_as_angle', default = clean_call.scales_as_angle)
            if type(clean_call.scales_as_angle) is str:
                clean_call.scales_as_angle = np.array(list(filter(None, clean_call.scales_as_angle.split(',')))).astype(float).tolist()
        #logger.debug(str(clean_call))
        
        # Define the clean mask (note one mask per galaxy, so we need to convert target multipart name to target name)
        dir_key = self._kh._dir_keys # read_dir_key() #<TODO># Need KeyHandler function get_target_name_by_multipart_name()
        if target in dir_key:
            target_name = dir_key[target]
        else:
            target_name = target
        
        this_cleanmask = self._kh.get_cleanmask_filename(target = target_name, product = product)
        if this_cleanmask is None:
            #logger.warning('Warning! Clean mask is not defined for target "'+target+'" in "cleanmask_key.txt"! cleanmask_dict: '+str(cleanmask_dict.keys()))
            this_cleanmask = ''
        
        if os.path.isfile(this_cleanmask):
            clean_call.clean_mask_file = this_cleanmask
        else:
            clean_call.clean_mask_file = None
            #logger.warning('Warning! Clean mask file for target "'+target_name+'" and product "'+product+'" was not found: "'+this_cleanmask+'"')
        
        # 
        # Return
        return clean_call
    
    # end of recipe_build_clean_call()
      
    #############################
    # recipe_imaging_one_target #
    #############################
    
    def recipe_phangsalma_imaging(
        self, 
        target=None,
        product=None,
        config=None,
        extra_ext_in = None,
        suffix_in = None,
        extra_ext_out = None,
        do_dirty_image = True,
        do_revert_to_dirty = True,
        do_read_clean_mask = True, 
        do_multiscale_clean = True
        do_revert_to_multiscale = True,
        do_singlscale_mask = True,
        do_singlescale_clean = True        
        do_export_to_fits = True, 
        dynamic_sizing = True,
        force_square = False,
        export_dirty = False,
        export_multiscale = False,
        overwrite = False, 
        ): 
        """
        PHANGS-ALMA basic imaging recipe. 
        
        Major steps:
       
        (1) Dirty imaging
        (2) Align a broad clean mask (if present)
        (3) Lightly masked multiscale clean to S/N ~ 4
        (4) Heavily masked single scale clean until convergence
        (5) Export to FITS

        Operates by passing a "clean_call" object between steps. Can
        restart from any individual step.
        """

        cell = None
        imsize = None

        # These calls instantiate a clean call object with default
        # clean parameters. For the next few steps, we use the clean
        # call appropriate for stage "dirty" imaging.

        clean_call = self.task_initialize_clean_call(
            target = target, config = config, product = product, 
            extra_ext_in = extra_ext_in, suffix_in = suffix_in,
            extra_ext_out = extra_ext_out,
            stage = 'dirty')
            
        if dynamic_sizing:
            if cell is None or imsize is None:
                cell, imsize = self.task_pick_cell_and_imsize(
                    clean_call=clean_call, force_square=force_square)
            else:                    
                clean_call.set_param('cell',cell,nowarning=True)
                clean_call.set_param('imsize',imsize,nowarning=True)                    

        # Make a dirty image (niter=0)

        if do_dirty_image:

            self.task_make_dirty_image(clean_call=clean_call)

        if do_export_to_fits and export_dirty:

            self.task_export_to_fits(clean_call=clean_call, tag='dirty')

        # Reset the current imaging to the dirty image.

        if do_revert_to_dirty:            

            self.task_revert_to_imaging(clean_call=clean_call, tag='dirty')

        # Read and align the clean mask to the astrometry of the image.

        if do_read_clean_mask:            

            self.task_read_clean_mask(
                clean_call = clean_call,
                target = target, config = config, product = product)

        # For the next few steps, we use the clean call appropriate
        # for stage "multiscale" imaging.

        clean_call = self.task_initialize_clean_call(
            target = target, config = config, product = product, 
            extra_ext_in = extra_ext_in, suffix_in = suffix_in,
            extra_ext_out = extra_ext_out,
            stage = 'multiscale')

        if dynamic_sizing:
            clean_call.set_param('cell',cell,nowarning=True)
            clean_call.set_param('imsize',imsize,nowarning=True)                    

        # Look up the angular scales to clean for this config
        
        self.task_assign_multiscales(config=config, clean_call=clean_call)
            
        # Run a multiscale clean until it converges.

        if do_multiscale_clean:
            
            self.task_multiscale_clean(clean_call=clean_call)

        if do_export_to_fits and export_multiscale:

            self.task_export_to_fits(clean_call=clean_call, tag='multiscale')

        # Reset the current imaging to the results of the multiscale clean.

        if do_revert_to_multiscale:
                        
            self.task_revert_to_imaging(clean_call=clean_call, tag='multiscale')

        # For the next few steps, we use the clean call appropriate
        # for stage "singlescale" imaging.

        clean_call = self.task_initialize_clean_call(
            target = target, config = config, product = product, 
            extra_ext_in = extra_ext_in, suffix_in = suffix_in,
            extra_ext_out = extra_ext_out,
            stage = 'singlescale')

        if dynamic_sizing:
            clean_call.set_param('cell',cell,nowarning=True)
            clean_call.set_param('imsize',imsize,nowarning=True)                    

        # Make a signal-to-noise based mask for use in singlescale clean.

        if do_singlescale_mask:
            
            self.task_singlescale_mask(clean_call=clean_call)

        # Run a singlescale clean until it converges.

        if do_singlescale_clean:

            self.task_singlescale_clean(clean_call=clean_call)

        # Export the products of the current clean to FITS files.

        if do_export_to_fits:

            self.task_export_to_fits(clean_call=clean_call)
    
    # end of recipe_phangsalma_imaging()
    
    


    
    
    
    
    
    







